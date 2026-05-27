#include <GCore/algorithms/gpu_geometry.h>

#ifdef GPU_GEOM_ALGORITHM

#include <RHI/rhi.hpp>
#include <spdlog/spdlog.h>

#include <memory>

#include "RHI/shaderCompiler.h"
#include "RHI/ShaderFactory/shader.hpp"

RUZINO_NAMESPACE_OPEN_SCOPE

static ResourceAllocator resource_allocator_;
static std::shared_ptr<ShaderFactory> shader_factory;

ResourceAllocator& get_resource_allocator()
{
    init_gpu_geometry_algorithms();
    return resource_allocator_;
}

void init_gpu_geometry_algorithms()
{
    if (shader_factory)
        return;

    // Hold an RHI reference so the device outlives all geometry GPU resources.
    // This reference is released in deinit after the cache is cleared.
    RHI::init(false);

    resource_allocator_.set_device(RHI::get_device());
    shader_factory = std::make_shared<ShaderFactory>();
    shader_factory->add_search_path(
        SlangShaderCompiler::get_shader_dir(ShaderDirType::Renderer)
            .string() +
        "/shaders");
    shader_factory->add_search_path(
        SlangShaderCompiler::get_shader_dir(ShaderDirType::GeomNodes)
            .string());
    shader_factory->add_search_path(
        SlangShaderCompiler::get_shader_dir(ShaderDirType::GeomCompute)
            .string());
    resource_allocator_.shader_factory = shader_factory.get();
}

void deinit_gpu_geometry_algorithms()
{
    if (!shader_factory)
        return;

    // Release all cached GPU resources while the device is still alive
    // (kept alive by our RHI reference).
    resource_allocator_.terminate();
    resource_allocator_.device = nullptr;
    resource_allocator_.shader_factory = nullptr;
    shader_factory.reset();

    // Release our RHI reference — the device may be destroyed here.
    RHI::shutdown();
}

RUZINO_NAMESPACE_CLOSE_SCOPE

#endif // GPU_GEOM_ALGORITHM
