#include <GCore/algorithms/intersection.h>

#include "GCore/Components/MeshComponent.h"
#include "GCore/Components/PointsComponent.h"
#include "GCore/Components/XformComponent.h"

#ifdef GPU_GEOM_ALGORITHM
#include "GPUContext/compute_context.hpp"
#include "GPUContext/program_vars.hpp"
#include "GPUContext/raytracing_context.hpp"
#include "RHI/ResourceManager/resource_allocator.hpp"
#include "RHI/internal/resources.hpp"
#include "RHI/rhi.hpp"
#include "Scene/SceneTypes.slang"
#include "glm/ext/matrix_transform.hpp"
#include "glm/gtc/type_ptr.hpp"
#include "nvrhi/nvrhi.h"

RUZINO_NAMESPACE_OPEN_SCOPE
ResourceAllocator resource_allocator_;
std::shared_ptr<ShaderFactory> shader_factory;

ResourceAllocator& get_resource_allocator()
{
    init_gpu_geometry_algorithms();
    return resource_allocator_;
}

void init_gpu_geometry_algorithms()
{
    if (!shader_factory) {
        resource_allocator_.set_device(RHI::get_device());
        shader_factory = std::make_shared<ShaderFactory>();
        shader_factory->add_search_path(RENDERER_SHADER_DIR "shaders");
        shader_factory->add_search_path(GEOM_COMPUTE_SHADER_DIR);
        resource_allocator_.shader_factory = shader_factory.get();
    }
}

void deinit_gpu_geometry_algorithms()
{
    resource_allocator_.terminate();
    shader_factory.reset();
}

nvrhi::rt::AccelStructHandle get_geomtry_tlas(
    const Geometry& geometry,
    MeshDesc* out_mesh_desc,
    nvrhi::BufferHandle* out_vertex_buffer)
{
    init_gpu_geometry_algorithms();
    auto mesh_component = geometry.get_component<MeshComponent>();
    if (!mesh_component) {
        return nullptr;
    }

    auto transform = glm::identity<glm::mat4>();
    auto xform_component = geometry.get_component<XformComponent>();
    if (xform_component) {
        transform = xform_component->get_transform();
    }

    // First build the BLAS for the mesh
    MeshDesc mesh_desc;
    auto device = RHI::get_device();

    size_t total_buffer_size = 0;
    size_t index_buffer_offset = 0;
    size_t normal_buffer_offset = 0;
    size_t texcoord_buffer_offset = 0;

    auto vertices = mesh_component->get_vertices();
    auto indices = mesh_component->get_face_vertex_indices();
    auto normals = mesh_component->get_normals();
    auto uvs = mesh_component->get_texcoords_array();

    // Calculate buffer offsets and total size
    total_buffer_size = vertices.size() * 3 * sizeof(float);
    index_buffer_offset = total_buffer_size;
    total_buffer_size += indices.size() * sizeof(int);

    if (!normals.empty()) {
        normal_buffer_offset = total_buffer_size;
        total_buffer_size += normals.size() * 3 * sizeof(float);
    }

    if (!uvs.empty()) {
        texcoord_buffer_offset = total_buffer_size;
        total_buffer_size += uvs.size() * 2 * sizeof(float);
    }

    // Create vertex buffer
    nvrhi::BufferDesc desc =
        nvrhi::BufferDesc{}
            .setCanHaveRawViews(true)
            .setByteSize(total_buffer_size)
            .setIsVertexBuffer(true)
            .setInitialState(nvrhi::ResourceStates::ShaderResource)
            .setCpuAccess(nvrhi::CpuAccessMode::None)
            .setIsAccelStructBuildInput(true)
            .setKeepInitialState(true)
            .setDebugName("meshVertexBuffer");

    auto vertexBuffer = device->createBuffer(desc);

    // Create single command list for all operations
    auto commandlist = resource_allocator_.create(CommandListDesc{});
    commandlist->open();

    // Copy vertices
    commandlist->writeBuffer(
        vertexBuffer, vertices.data(), vertices.size() * 3 * sizeof(float), 0);

    // Copy indices
    commandlist->writeBuffer(
        vertexBuffer,
        indices.data(),
        indices.size() * sizeof(int),
        index_buffer_offset);

    // Copy normals if available
    if (!normals.empty()) {
        commandlist->writeBuffer(
            vertexBuffer,
            normals.data(),
            normals.size() * 3 * sizeof(float),
            normal_buffer_offset);
    }

    // Copy UVs if available
    if (!uvs.empty()) {
        commandlist->writeBuffer(
            vertexBuffer,
            uvs.data(),
            uvs.size() * 2 * sizeof(float),
            texcoord_buffer_offset);
    }

    commandlist->close();
    device->executeCommandList(commandlist);
    device->waitForIdle();

    // Set up mesh_desc
    mesh_desc.vbOffset = 0;
    mesh_desc.ibOffset = index_buffer_offset;
    mesh_desc.normalOffset = normal_buffer_offset;
    mesh_desc.texCrdOffset = texcoord_buffer_offset;
    mesh_desc.normalInterpolation = InterpolationType::Vertex;
    mesh_desc.texCrdInterpolation = InterpolationType::Vertex;

    // Create BLAS
    nvrhi::rt::AccelStructDesc blas_desc;
    nvrhi::rt::GeometryDesc geometry_desc;
    geometry_desc.geometryType = nvrhi::rt::GeometryType::Triangles;
    nvrhi::rt::GeometryTriangles triangles;
    triangles.setVertexBuffer(vertexBuffer)
        .setVertexOffset(0)
        .setIndexBuffer(vertexBuffer)
        .setIndexOffset(index_buffer_offset)
        .setIndexCount(indices.size())
        .setVertexCount(vertices.size())
        .setVertexStride(3 * sizeof(float))
        .setVertexFormat(nvrhi::Format::RGB32_FLOAT)
        .setIndexFormat(nvrhi::Format::R32_UINT);
    geometry_desc.setTriangles(triangles);
    blas_desc.addBottomLevelGeometry(geometry_desc);
    blas_desc.isTopLevel = false;
    auto BLAS = device->createAccelStruct(blas_desc);

    commandlist->open();
    nvrhi::utils::BuildBottomLevelAccelStruct(
        commandlist, BLAS, blas_desc);
    commandlist->close();
    device->executeCommandList(commandlist);
    device->waitForIdle();

    // Now create TLAS
    nvrhi::rt::AccelStructDesc tlas_desc;
    tlas_desc.isTopLevel = true;
    tlas_desc.topLevelMaxInstances = 1;
    nvrhi::rt::InstanceDesc instance_desc;
    instance_desc.setBLAS(BLAS);
    instance_desc.setInstanceID(0);
    instance_desc.setInstanceMask(0xFF);
    // set the transform
    nvrhi::rt::AffineTransform affine_transform;

    memcpy(
        affine_transform,
        glm::value_ptr(transpose(transform)),
        sizeof(nvrhi::rt::AffineTransform));

    instance_desc.setTransform(affine_transform);

    auto TLAS = device->createAccelStruct(tlas_desc);

    commandlist->open();
    commandlist->buildTopLevelAccelStruct(
        TLAS, std::vector{ instance_desc }.data(), 1);
    commandlist->close();
    device->executeCommandList(commandlist);
    device->waitForIdle();

    // Clean up commandlist
    resource_allocator_.destroy(commandlist);

    // Fill output parameters if provided
    if (out_mesh_desc) {
        *out_mesh_desc = mesh_desc;
    }

    if (out_vertex_buffer) {
        *out_vertex_buffer = vertexBuffer;
    }

    assert(TLAS);
    return TLAS;
}

nvrhi::BufferHandle IntersectToBuffer(
    const nvrhi::BufferHandle& ray_buffer,
    size_t ray_count,
    const Geometry& BaseMesh)
{
    init_gpu_geometry_algorithms();
    auto mesh_component = BaseMesh.get_component<MeshComponent>();

    if (!mesh_component) {
        return nullptr;
    }

    std::vector<glm::vec3> vertices = mesh_component->get_vertices();

    if (vertices.empty()) {
        return nullptr;
    }

    std::vector<glm::vec3> normals = mesh_component->get_normals();
    std::vector<int> indices = mesh_component->get_face_vertex_indices();
    std::vector<glm::vec2> uvs = mesh_component->get_texcoords_array();

    MeshDesc mesh_desc;

    auto& resource_allocator = get_resource_allocator();
    auto device = RHI::get_device();

    nvrhi::BufferHandle vertex_buffer;
    auto accel =
        get_geomtry_tlas(BaseMesh, &mesh_desc, std::addressof(vertex_buffer));

    auto commandlist = resource_allocator.create(CommandListDesc{});

    ProgramDesc desc;
    desc.shaderType = nvrhi::ShaderType::AllRayTracing;
    desc.set_path(GEOM_COMPUTE_SHADER_DIR "intersection.slang");
    auto program = resource_allocator.create(desc);

    // Create output buffer for intersection results
    auto result_buffer = resource_allocator.create(
        nvrhi::BufferDesc{}
            .setByteSize(ray_count * sizeof(PointSample))
            .setStructStride(sizeof(PointSample))
            .setInitialState(nvrhi::ResourceStates::ShaderResource)
            .setKeepInitialState(true)
            .setCanHaveUAVs(true)
            .setDebugName("resultBuffer"));

    auto mesh_cb = resource_allocator.create(
        nvrhi::BufferDesc{}
            .setByteSize(sizeof(MeshDesc))
            .setStructStride(sizeof(MeshDesc))
            .setInitialState(nvrhi::ResourceStates::ConstantBuffer)
            .setKeepInitialState(true)
            .setCpuAccess(nvrhi::CpuAccessMode::Write)
            .setIsConstantBuffer(true)
            .setDebugName("meshDescBuffer"));

    void* data = device->mapBuffer(mesh_cb, nvrhi::CpuAccessMode::Write);
    memcpy(data, &mesh_desc, sizeof(MeshDesc));
    device->unmapBuffer(mesh_cb);

    ProgramVars program_vars(resource_allocator, program);
    // Bind resources to the shader
    program_vars["mesh"] = mesh_cb;
    program_vars["g_vertexBuffer"] = vertex_buffer;
    program_vars["g_rays"] = ray_buffer;
    program_vars["g_Result"] = result_buffer;
    program_vars["SceneBVH"] = accel;

    program_vars.finish_setting_vars();

    RaytracingContext raytracing_context(resource_allocator, program_vars);

    // Set up shader names
    raytracing_context.announce_raygeneration("RayGen");
    raytracing_context.announce_hitgroup("ClosestHit");
    raytracing_context.announce_miss("Miss");
    raytracing_context.finish_announcing_shader_names();

    // Execute ray tracing
    raytracing_context.begin();
    raytracing_context.trace_rays({}, program_vars, ray_count, 1, 1);
    raytracing_context.finish();

    // Clean up resources except for result_buffer
    resource_allocator.destroy(vertex_buffer);
    resource_allocator.destroy(mesh_cb);
    resource_allocator.destroy(program);
    resource_allocator.destroy(commandlist);

    // Return the GPU buffer with results
    return result_buffer;
}

std::vector<PointSample> IntersectWithBuffer(
    const nvrhi::BufferHandle& ray_buffer,
    size_t ray_count,
    const Geometry& BaseMesh)
{
    auto& resource_allocator = get_resource_allocator();
    auto device = RHI::get_device();

    // Get the GPU buffer with intersection results
    auto result_buffer = IntersectToBuffer(ray_buffer, ray_count, BaseMesh);

    // If we got a null buffer, return empty results
    if (!result_buffer) {
        return std::vector<PointSample>(ray_count);
    }

    std::vector<PointSample> result;
    result.resize(ray_count);

    // Create readback buffer
    auto readback_buffer = resource_allocator.create(
        nvrhi::BufferDesc{}
            .setByteSize(ray_count * sizeof(PointSample))
            .setCpuAccess(nvrhi::CpuAccessMode::Read)
            .setDebugName("resultReadbackBuffer"));

    // Create command list to copy data
    auto commandlist = resource_allocator.create(CommandListDesc{});
    commandlist->open();
    commandlist->copyBuffer(
        readback_buffer, 0, result_buffer, 0, ray_count * sizeof(PointSample));
    commandlist->close();
    device->executeCommandList(commandlist);
    device->waitForIdle();

    // Map and read the results
    void* mapped_data =
        device->mapBuffer(readback_buffer, nvrhi::CpuAccessMode::Read);
    memcpy(result.data(), mapped_data, ray_count * sizeof(PointSample));
    device->unmapBuffer(readback_buffer);

    // Clean up resources
    resource_allocator.destroy(result_buffer);
    resource_allocator.destroy(readback_buffer);
    resource_allocator.destroy(commandlist);

    return result;
}

std::vector<PointSample> Intersect(
    const std::vector<glm::ray>& rays,
    const Geometry& BaseMesh)
{
    init_gpu_geometry_algorithms();
    auto& resource_allocator = get_resource_allocator();
    auto device = RHI::get_device();

    // Create ray buffer with input rays
    auto ray_buffer = resource_allocator.create(
        nvrhi::BufferDesc{}
            .setByteSize(rays.size() * sizeof(glm::ray))
            .setStructStride(sizeof(glm::ray))
            .setInitialState(nvrhi::ResourceStates::ShaderResource)
            .setKeepInitialState(true)
            .setCpuAccess(nvrhi::CpuAccessMode::Write)
            .setDebugName("rayBuffer"));

    // Copy rays to the ray buffer
    void* data = device->mapBuffer(ray_buffer, nvrhi::CpuAccessMode::Write);
    memcpy(data, rays.data(), rays.size() * sizeof(glm::ray));
    device->unmapBuffer(ray_buffer);

    // Call the implementation with the buffer
    auto result = IntersectWithBuffer(ray_buffer, rays.size(), BaseMesh);

    // Clean up the buffer we created
    resource_allocator.destroy(ray_buffer);

    return result;
}

std::vector<PointSample> Intersect(
    const std::vector<glm::vec3>& start_point,
    const std::vector<glm::vec3>& next_point,
    const Geometry& BaseMesh)
{
    std::vector<glm::ray> rays;
    rays.reserve(start_point.size());
    for (size_t i = 0; i < start_point.size(); ++i) {
        rays.push_back(
            glm::ray(start_point[i], next_point[i] - start_point[i]));
    }

    return Intersect(rays, BaseMesh);
}

std::vector<PointSample> IntersectInterweaved(
    const std::vector<glm::vec3>& start_point,
    const std::vector<glm::vec3>& next_point,
    const Geometry& BaseMesh)
{
    std::vector<glm::ray> rays;
    rays.reserve(start_point.size() * next_point.size());
    for (size_t i = 0; i < start_point.size(); ++i) {
        for (size_t j = 0; j < next_point.size(); ++j) {
            rays.push_back(
                glm::ray(start_point[i], next_point[j] - start_point[i]));
        }
    }
    return Intersect(rays, BaseMesh);
}

nvrhi::BufferHandle FindNeighborsFromPositionBuffer(
    const nvrhi::BufferHandle& position_buffer,
    size_t point_count,
    float radius,
    unsigned& out_pair_count)
{
    init_gpu_geometry_algorithms();
    
    if (!position_buffer || point_count == 0) {
        out_pair_count = 0;
        return nullptr;
    }

    auto& resource_allocator = get_resource_allocator();
    auto device = RHI::get_device();

    // Create single command list for all operations
    auto commandlist = resource_allocator.create(CommandListDesc{});

    // Step 1: Create AABB buffer and compute AABBs using toAABB.slang
    auto aabb_buffer = resource_allocator.create(
        nvrhi::BufferDesc{}
            .setByteSize(point_count * 2 * sizeof(glm::vec3))
            .setInitialState(nvrhi::ResourceStates::UnorderedAccess)
            .setKeepInitialState(true)
            .setCanHaveRawViews(true)
            .setCanHaveUAVs(true)
            .setIsAccelStructBuildInput(true)
            .setDebugName("aabbBuffer"));

    // Compile and run toAABB compute shader
    ProgramDesc aabb_desc;
    aabb_desc.set_path(GEOM_COMPUTE_SHADER_DIR "Points/toAABB.slang")
        .set_entry_name("main")
        .set_shader_type(nvrhi::ShaderType::Compute);

    auto aabb_program = resource_allocator.create(aabb_desc);

    if (!aabb_program->get_error_string().empty()) {
        spdlog::error(
            "toAABB shader compilation failed: {}",
            aabb_program->get_error_string());
        resource_allocator.destroy(aabb_program);
        resource_allocator.destroy(aabb_buffer);
        resource_allocator.destroy(commandlist);
        out_pair_count = 0;
        return nullptr;
    }

    ProgramVars aabb_vars(resource_allocator, aabb_program);

    // Create constant buffer for radius
    auto radius_buffer = resource_allocator.create(
        nvrhi::BufferDesc{}
            .setByteSize(16)
            .setIsConstantBuffer(true)
            .setInitialState(nvrhi::ResourceStates::ConstantBuffer)
            .setKeepInitialState(true)
            .setDebugName("radiusBuffer"));

    commandlist->open();
    float radius_data[4] = { radius, 0, 0, 0 };
    commandlist->writeBuffer(radius_buffer, radius_data, sizeof(radius_data));
    commandlist->close();
    device->executeCommandList(commandlist);
    device->waitForIdle();

    aabb_vars["SearchParams"] = radius_buffer;
    aabb_vars["positions"] = position_buffer;
    aabb_vars["AABBs"] = aabb_buffer;
    aabb_vars.finish_setting_vars();

    ComputeContext compute_context(resource_allocator, aabb_vars);
    compute_context.finish_setting_pso();
    compute_context.begin();
    compute_context.dispatch({}, aabb_vars, point_count, 32);
    compute_context.finish();

    resource_allocator.destroy(aabb_program);
    resource_allocator.destroy(radius_buffer);

    // Step 2: Build BLAS with AABBs
    nvrhi::rt::AccelStructDesc blas_desc;
    nvrhi::rt::GeometryDesc geometry_desc;
    geometry_desc.geometryType = nvrhi::rt::GeometryType::AABBs;
    geometry_desc.flags = nvrhi::rt::GeometryFlags::NoDuplicateAnyHitInvocation;
    nvrhi::rt::GeometryAABBs aabbs;
    aabbs.setBuffer(aabb_buffer)
        .setOffset(0)
        .setCount(point_count)
        .setStride(2 * sizeof(glm::vec3));
    geometry_desc.setAABBs(aabbs);
    blas_desc.addBottomLevelGeometry(geometry_desc);
    blas_desc.isTopLevel = false;
    auto BLAS = resource_allocator.create(blas_desc);

    commandlist->open();
    nvrhi::utils::BuildBottomLevelAccelStruct(commandlist, BLAS, blas_desc);
    commandlist->close();
    device->executeCommandList(commandlist);
    device->waitForIdle();

    // Step 3: Build TLAS
    nvrhi::rt::AccelStructDesc tlas_desc;
    tlas_desc.isTopLevel = true;
    tlas_desc.topLevelMaxInstances = 1;

    nvrhi::rt::InstanceDesc instance_desc;
    instance_desc.setBLAS(BLAS);
    instance_desc.setInstanceID(0);
    instance_desc.setInstanceMask(0xFF);

    nvrhi::rt::AffineTransform affine_transform;
    memset(affine_transform, 0, sizeof(nvrhi::rt::AffineTransform));
    affine_transform[0] = 1.0f;
    affine_transform[5] = 1.0f;
    affine_transform[10] = 1.0f;
    instance_desc.setTransform(affine_transform);

    auto TLAS = resource_allocator.create(tlas_desc);

    commandlist->open();
    commandlist->buildTopLevelAccelStruct(
        TLAS, std::vector{ instance_desc }.data(), 1);
    commandlist->close();
    device->executeCommandList(commandlist);
    device->waitForIdle();

    // Step 4: Create output buffers for pairs
    size_t max_pairs = point_count * 30;

    auto pairs_buffer = resource_allocator.create(
        nvrhi::BufferDesc{}
            .setByteSize(max_pairs * sizeof(PointPairs))
            .setStructStride(sizeof(PointPairs))
            .setInitialState(nvrhi::ResourceStates::UnorderedAccess)
            .setKeepInitialState(true)
            .setCanHaveUAVs(true)
            .setDebugName("pairsBuffer"));

    auto pairs_count_buffer = resource_allocator.create(
        nvrhi::BufferDesc{}
            .setByteSize(sizeof(unsigned))
            .setStructStride(sizeof(unsigned))
            .setInitialState(nvrhi::ResourceStates::UnorderedAccess)
            .setKeepInitialState(true)
            .setCanHaveUAVs(true)
            .setDebugName("pairsCountBuffer"));

    commandlist->open();
    uint32_t zero = 0;
    commandlist->writeBuffer(pairs_count_buffer, &zero, sizeof(uint32_t));
    commandlist->close();
    device->executeCommandList(commandlist);
    device->waitForIdle();

    // Step 5: Run contact.slang ray tracing shader
    ProgramDesc contact_desc;
    contact_desc.shaderType = nvrhi::ShaderType::AllRayTracing;
    contact_desc.set_path(GEOM_COMPUTE_SHADER_DIR "Points/contact.slang");
    auto contact_program = resource_allocator.create(contact_desc);

    if (!contact_program || !contact_program->getBufferPointer() ||
        !contact_program->get_error_string().empty()) {
        if (!contact_program->get_error_string().empty()) {
            spdlog::error(
                "contact shader compilation failed: {}",
                contact_program->get_error_string());
        }
        resource_allocator.destroy(contact_program);
        resource_allocator.destroy(pairs_count_buffer);
        resource_allocator.destroy(pairs_buffer);
        resource_allocator.destroy(TLAS);
        resource_allocator.destroy(BLAS);
        resource_allocator.destroy(aabb_buffer);
        resource_allocator.destroy(commandlist);
        out_pair_count = 0;
        return nullptr;
    }

    ProgramVars contact_vars(resource_allocator, contact_program);

    auto contact_radius_buffer = resource_allocator.create(
        nvrhi::BufferDesc{}
            .setByteSize(16)
            .setIsConstantBuffer(true)
            .setInitialState(nvrhi::ResourceStates::ConstantBuffer)
            .setKeepInitialState(true)
            .setDebugName("contactRadiusBuffer"));

    commandlist->open();
    float contact_radius_data[4] = { radius, 0, 0, 0 };
    commandlist->writeBuffer(
        contact_radius_buffer,
        contact_radius_data,
        sizeof(contact_radius_data));
    commandlist->close();
    device->executeCommandList(commandlist);
    device->waitForIdle();

    contact_vars["SearchParams"] = contact_radius_buffer;
    contact_vars["SceneBVH"] = TLAS;
    contact_vars["positions"] = position_buffer;
    contact_vars["Pairs"] = pairs_buffer;
    contact_vars["PairsCount"] = pairs_count_buffer;
    contact_vars.finish_setting_vars();

    RaytracingContext raytracing_context(resource_allocator, contact_vars);
    raytracing_context.announce_raygeneration("RayGen");
    raytracing_context.announce_hitgroup(
        "ClosestHit", "AnyHit", "Intersection");
    raytracing_context.announce_miss("Miss");
    raytracing_context.finish_announcing_shader_names();

    raytracing_context.begin();
    raytracing_context.trace_rays({}, contact_vars, point_count, 1, 1);
    raytracing_context.finish();

    // Step 6: Read back the pair count
    auto count_readback_buffer = resource_allocator.create(
        nvrhi::BufferDesc{}
            .setByteSize(sizeof(unsigned))
            .setCpuAccess(nvrhi::CpuAccessMode::Read)
            .setInitialState(nvrhi::ResourceStates::CopyDest)
            .setKeepInitialState(true)
            .setDebugName("pairsCountReadbackBuffer"));

    commandlist->open();
    commandlist->copyBuffer(
        count_readback_buffer, 0, pairs_count_buffer, 0, sizeof(unsigned));
    commandlist->close();
    device->executeCommandList(commandlist);
    device->waitForIdle();

    void* count_data =
        device->mapBuffer(count_readback_buffer, nvrhi::CpuAccessMode::Read);
    memcpy(&out_pair_count, count_data, sizeof(unsigned));
    device->unmapBuffer(count_readback_buffer);

    // Clean up temporary resources
    resource_allocator.destroy(aabb_buffer);
    resource_allocator.destroy(BLAS);
    resource_allocator.destroy(TLAS);
    resource_allocator.destroy(pairs_count_buffer);
    resource_allocator.destroy(contact_program);
    resource_allocator.destroy(contact_radius_buffer);
    resource_allocator.destroy(count_readback_buffer);
    resource_allocator.destroy(commandlist);

    return pairs_buffer;
}

nvrhi::BufferHandle FindNeighborsToBuffer(
    const Geometry& point_cloud,
    float radius,
    unsigned& out_pair_count)
{
    init_gpu_geometry_algorithms();
    auto points_component = point_cloud.get_component<PointsComponent>();

    if (!points_component) {
        out_pair_count = 0;
        return nullptr;
    }

    std::vector<glm::vec3> positions = points_component->get_vertices();

    if (positions.empty()) {
        out_pair_count = 0;
        return nullptr;
    }

    auto& resource_allocator = get_resource_allocator();
    auto device = RHI::get_device();

    size_t point_count = positions.size();

    // Create position buffer from geometry
    auto position_buffer = resource_allocator.create(
        nvrhi::BufferDesc{}
            .setByteSize(point_count * sizeof(glm::vec3))
            .setStructStride(sizeof(glm::vec3))
            .setInitialState(nvrhi::ResourceStates::ShaderResource)
            .setKeepInitialState(true)
            .setDebugName("positionBuffer"));

    auto commandlist = resource_allocator.create(CommandListDesc{});
    commandlist->open();
    commandlist->writeBuffer(
        position_buffer, positions.data(), point_count * sizeof(glm::vec3));
    commandlist->close();
    device->executeCommandList(commandlist);
    device->waitForIdle();
    resource_allocator.destroy(commandlist);

    // Call the core function with the position buffer
    auto result = FindNeighborsFromPositionBuffer(
        position_buffer, point_count, radius, out_pair_count);

    // Clean up the position buffer we created
    resource_allocator.destroy(position_buffer);

    return result;
}

std::vector<PointPairs> FindNeighbors(
    const Geometry& point_cloud,
    float radius,
    unsigned& out_pair_count)
{
    auto& resource_allocator = get_resource_allocator();
    auto device = RHI::get_device();

    // Get the GPU buffer with neighbor pairs
    auto pairs_buffer =
        FindNeighborsToBuffer(point_cloud, radius, out_pair_count);

    // If we got a null buffer or no pairs found, return empty results
    if (!pairs_buffer || out_pair_count == 0) {
        if (pairs_buffer) {
            resource_allocator.destroy(pairs_buffer);
        }
        return std::vector<PointPairs>();
    }

    std::vector<PointPairs> result;
    result.resize(out_pair_count);

    // Create readback buffer
    auto readback_buffer = resource_allocator.create(
        nvrhi::BufferDesc{}
            .setByteSize(out_pair_count * sizeof(PointPairs))
            .setCpuAccess(nvrhi::CpuAccessMode::Read)
            .setInitialState(nvrhi::ResourceStates::CopyDest)
            .setKeepInitialState(true)
            .setDebugName("pairsReadbackBuffer"));

    // Create command list to copy data
    auto commandlist = resource_allocator.create(CommandListDesc{});
    commandlist->open();
    commandlist->copyBuffer(
        readback_buffer,
        0,
        pairs_buffer,
        0,
        out_pair_count * sizeof(PointPairs));
    commandlist->close();
    device->executeCommandList(commandlist);
    device->waitForIdle();

    // Map and read the results
    void* mapped_data =
        device->mapBuffer(readback_buffer, nvrhi::CpuAccessMode::Read);
    memcpy(result.data(), mapped_data, out_pair_count * sizeof(PointPairs));
    device->unmapBuffer(readback_buffer);

    // Clean up resources
    resource_allocator.destroy(pairs_buffer);
    resource_allocator.destroy(readback_buffer);
    resource_allocator.destroy(commandlist);

    return result;
}

RUZINO_NAMESPACE_CLOSE_SCOPE
#endif
