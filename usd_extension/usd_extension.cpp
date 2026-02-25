#include <GCore/usd_extension.h>
#include <pxr/base/gf/matrix4d.h>
#include <pxr/usd/usdGeom/capsule.h>
#include <pxr/usd/usdGeom/cone.h>
#include <pxr/usd/usdGeom/cube.h>
#include <pxr/usd/usdGeom/curves.h>
#include <pxr/usd/usdGeom/cylinder.h>
#include <pxr/usd/usdGeom/mesh.h>
#include <pxr/usd/usdGeom/pointInstancer.h>
#include <pxr/usd/usdGeom/points.h>
#include <pxr/usd/usdGeom/primvarsAPI.h>
#include <pxr/usd/usdGeom/sphere.h>
#include <pxr/usd/usdShade/materialBindingAPI.h>
#include <pxr/usd/usdSkel/bindingAPI.h>
#include <pxr/usd/usdSkel/cache.h>
#include <pxr/usd/usdSkel/skeleton.h>
#include <pxr/usd/usdSkel/skeletonQuery.h>
#include <pxr/usd/usdVol/openVDBAsset.h>
#include <pxr/usd/usdVol/volume.h>
#include <spdlog/spdlog.h>

#include "GCore/Components/CurveComponent.h"
#include "GCore/Components/InstancerComponent.h"
#include "GCore/Components/MaterialComponent.h"
#include "GCore/Components/MeshComponent.h"
#include "GCore/Components/MeshUSDView.h"  // Only need USD view
#include "GCore/Components/PointsComponent.h"
#include "GCore/Components/SkelComponent.h"
#include "GCore/Components/VolumeComponent.h"
#include "GCore/Components/XformComponent.h"
#include "glm/glm.hpp"
#include "glm/gtc/type_ptr.hpp"

RUZINO_NAMESPACE_OPEN_SCOPE

// Helper function to check if prim is a basic shape type that can be converted
// to mesh
static bool is_convertible_geom_shape(const pxr::UsdPrim& prim)
{
    using namespace pxr;
    return prim.IsA<UsdGeomMesh>() || prim.IsA<UsdGeomSphere>() ||
           prim.IsA<UsdGeomCylinder>() || prim.IsA<UsdGeomCone>() ||
           prim.IsA<UsdGeomCube>() || prim.IsA<UsdGeomCapsule>();
}

bool generate_parametric_mesh(
    const std::string& type_name,
    const ParametricShapeParams& params,
    pxr::VtArray<pxr::GfVec3f>& out_points,
    pxr::VtArray<int>& out_face_counts,
    pxr::VtArray<int>& out_face_indices,
    pxr::VtArray<pxr::GfVec3f>& out_normals)
{
    using namespace pxr;
    float radius = static_cast<float>(params.radius);
    float height = static_cast<float>(params.height);
    float size = static_cast<float>(params.size);

    if (type_name == "Sphere" || type_name == "UsdGeomSphere") {
        int segments = 16;
        int rings = 8;

        for (int ring = 0; ring <= rings; ++ring) {
            float phi = M_PI * float(ring) / float(rings);
            float y = cosf(phi) * radius;
            float ring_radius = sinf(phi) * radius;

            for (int seg = 0; seg <= segments; ++seg) {
                float theta = 2.0f * M_PI * float(seg) / float(segments);
                float x = cosf(theta) * ring_radius;
                float z = sinf(theta) * ring_radius;
                out_points.push_back(GfVec3f(x, y, z));
                out_normals.push_back(
                    GfVec3f(x / radius, y / radius, z / radius));
            }
        }

        for (int ring = 0; ring < rings; ++ring) {
            for (int seg = 0; seg < segments; ++seg) {
                int first = ring * (segments + 1) + seg;
                int second = first + segments + 1;

                out_face_counts.push_back(3);
                out_face_indices.push_back(first);
                out_face_indices.push_back(second);
                out_face_indices.push_back(first + 1);

                out_face_counts.push_back(3);
                out_face_indices.push_back(first + 1);
                out_face_indices.push_back(second);
                out_face_indices.push_back(second + 1);
            }
        }
        return true;
    }

    if (type_name == "Cylinder" || type_name == "UsdGeomCylinder") {
        float half_height = height / 2.0f;
        int segments = 16;

        int top_center_idx = 0;
        int bottom_center_idx = 1;
        out_points.push_back(GfVec3f(0, half_height, 0));
        out_points.push_back(GfVec3f(0, -half_height, 0));
        out_normals.push_back(GfVec3f(0, 1, 0));
        out_normals.push_back(GfVec3f(0, -1, 0));

        for (int seg = 0; seg <= segments; ++seg) {
            float theta = 2.0f * M_PI * float(seg) / float(segments);
            float x = cosf(theta) * radius;
            float z = sinf(theta) * radius;
            out_points.push_back(GfVec3f(x, half_height, z));
            out_normals.push_back(GfVec3f(0, 1, 0));
        }

        int bottom_ring_start = out_points.size();
        for (int seg = 0; seg <= segments; ++seg) {
            float theta = 2.0f * M_PI * float(seg) / float(segments);
            float x = cosf(theta) * radius;
            float z = sinf(theta) * radius;
            out_points.push_back(GfVec3f(x, -half_height, z));
            out_normals.push_back(GfVec3f(0, -1, 0));
        }

        for (int seg = 0; seg < segments; ++seg) {
            out_face_counts.push_back(3);
            out_face_indices.push_back(top_center_idx);
            out_face_indices.push_back(2 + seg);
            out_face_indices.push_back(2 + seg + 1);
        }

        for (int seg = 0; seg < segments; ++seg) {
            out_face_counts.push_back(3);
            out_face_indices.push_back(bottom_center_idx);
            out_face_indices.push_back(bottom_ring_start + seg + 1);
            out_face_indices.push_back(bottom_ring_start + seg);
        }

        int side_top_start = out_points.size();
        for (int seg = 0; seg <= segments; ++seg) {
            float theta = 2.0f * M_PI * float(seg) / float(segments);
            float x = cosf(theta) * radius;
            float z = sinf(theta) * radius;
            out_points.push_back(GfVec3f(x, half_height, z));
            out_normals.push_back(GfVec3f(x / radius, 0, z / radius));
        }
        int side_bottom_start = out_points.size();
        for (int seg = 0; seg <= segments; ++seg) {
            float theta = 2.0f * M_PI * float(seg) / float(segments);
            float x = cosf(theta) * radius;
            float z = sinf(theta) * radius;
            out_points.push_back(GfVec3f(x, -half_height, z));
            out_normals.push_back(GfVec3f(x / radius, 0, z / radius));
        }

        for (int seg = 0; seg < segments; ++seg) {
            out_face_counts.push_back(4);
            out_face_indices.push_back(side_top_start + seg);
            out_face_indices.push_back(side_bottom_start + seg);
            out_face_indices.push_back(side_bottom_start + seg + 1);
            out_face_indices.push_back(side_top_start + seg + 1);
        }

        return true;
    }

    if (type_name == "Cube" || type_name == "UsdGeomCube") {
        float half = size / 2.0f;

        out_points = { { -half, -half, -half }, { half, -half, -half },
                       { half, half, -half },   { -half, half, -half },
                       { -half, -half, half },  { half, -half, half },
                       { half, half, half },    { -half, half, half } };

        out_face_counts = { 4, 4, 4, 4, 4, 4 };
        out_face_indices = { 0, 1, 2, 3, 5, 4, 7, 6, 4, 0, 3, 7,
                             1, 5, 6, 2, 3, 2, 6, 7, 4, 5, 1, 0 };

        out_normals = { { 0, 0, -1 }, { 0, 0, -1 }, { 0, 0, -1 }, { 0, 0, -1 },
                        { 0, 0, 1 },  { 0, 0, 1 },  { 0, 0, 1 },  { 0, 0, 1 } };

        return true;
    }

    if (type_name == "Cone" || type_name == "UsdGeomCone") {
        float half_height = height / 2.0f;
        int segments = 16;

        int apex_idx = 0;
        out_points.push_back(GfVec3f(0, half_height, 0));
        out_normals.push_back(GfVec3f(0, 1, 0));

        int bottom_center_idx = 1;
        out_points.push_back(GfVec3f(0, -half_height, 0));
        out_normals.push_back(GfVec3f(0, -1, 0));

        int bottom_ring_start = out_points.size();
        for (int seg = 0; seg <= segments; ++seg) {
            float theta = 2.0f * M_PI * float(seg) / float(segments);
            float x = cosf(theta) * radius;
            float z = sinf(theta) * radius;
            out_points.push_back(GfVec3f(x, -half_height, z));
            out_normals.push_back(GfVec3f(0, -1, 0));
        }

        for (int seg = 0; seg < segments; ++seg) {
            out_face_counts.push_back(3);
            out_face_indices.push_back(bottom_center_idx);
            out_face_indices.push_back(bottom_ring_start + seg + 1);
            out_face_indices.push_back(bottom_ring_start + seg);
        }

        for (int seg = 0; seg < segments; ++seg) {
            float theta1 = 2.0f * M_PI * float(seg) / float(segments);
            float theta2 = 2.0f * M_PI * float(seg + 1) / float(segments);

            GfVec3f p1(
                cosf(theta1) * radius, -half_height, sinf(theta1) * radius);
            GfVec3f p2(
                cosf(theta2) * radius, -half_height, sinf(theta2) * radius);
            GfVec3f apex(0, half_height, 0);

            GfVec3f edge1 = p1 - apex;
            GfVec3f edge2 = p2 - apex;
            GfVec3f normal = GfCross(edge1, edge2);
            normal.Normalize();

            int side_p1 = out_points.size();
            out_points.push_back(GfVec3f(p1));
            out_normals.push_back(normal);

            int side_p2 = side_p1 + 1;
            out_points.push_back(GfVec3f(p2));
            out_normals.push_back(normal);

            int side_apex = side_p2 + 1;
            out_points.push_back(GfVec3f(apex));
            out_normals.push_back(normal);

            out_face_counts.push_back(3);
            out_face_indices.push_back(side_p1);
            out_face_indices.push_back(side_p2);
            out_face_indices.push_back(side_apex);
        }

        return true;
    }

    if (type_name == "Capsule" || type_name == "UsdGeomCapsule") {
        float half_height = height / 2.0f;
        int segments = 16;

        int top_center_idx = 0;
        int bottom_center_idx = 1;
        out_points.push_back(GfVec3f(0, half_height + radius, 0));
        out_points.push_back(GfVec3f(0, -half_height - radius, 0));
        out_normals.push_back(GfVec3f(0, 1, 0));
        out_normals.push_back(GfVec3f(0, -1, 0));

        for (int seg = 0; seg <= segments; ++seg) {
            float theta = 2.0f * M_PI * float(seg) / float(segments);
            float x = cosf(theta) * radius;
            float z = sinf(theta) * radius;
            out_points.push_back(GfVec3f(x, half_height + radius, z));
            out_normals.push_back(GfVec3f(0, 1, 0));
        }

        int bottom_ring_start = out_points.size();
        for (int seg = 0; seg <= segments; ++seg) {
            float theta = 2.0f * M_PI * float(seg) / float(segments);
            float x = cosf(theta) * radius;
            float z = sinf(theta) * radius;
            out_points.push_back(GfVec3f(x, -half_height - radius, z));
            out_normals.push_back(GfVec3f(0, -1, 0));
        }

        for (int seg = 0; seg < segments; ++seg) {
            out_face_counts.push_back(3);
            out_face_indices.push_back(top_center_idx);
            out_face_indices.push_back(2 + seg);
            out_face_indices.push_back(2 + seg + 1);
        }

        for (int seg = 0; seg < segments; ++seg) {
            out_face_counts.push_back(3);
            out_face_indices.push_back(bottom_center_idx);
            out_face_indices.push_back(bottom_ring_start + seg + 1);
            out_face_indices.push_back(bottom_ring_start + seg);
        }

        int side_top_start = out_points.size();
        for (int seg = 0; seg <= segments; ++seg) {
            float theta = 2.0f * M_PI * float(seg) / float(segments);
            float x = cosf(theta) * radius;
            float z = sinf(theta) * radius;
            out_points.push_back(GfVec3f(x, half_height + radius, z));
            out_normals.push_back(GfVec3f(x / radius, 0, z / radius));
        }
        int side_bottom_start = out_points.size();
        for (int seg = 0; seg <= segments; ++seg) {
            float theta = 2.0f * M_PI * float(seg) / float(segments);
            float x = cosf(theta) * radius;
            float z = sinf(theta) * radius;
            out_points.push_back(GfVec3f(x, -half_height - radius, z));
            out_normals.push_back(GfVec3f(x / radius, 0, z / radius));
        }

        for (int seg = 0; seg < segments; ++seg) {
            out_face_counts.push_back(4);
            out_face_indices.push_back(side_top_start + seg);
            out_face_indices.push_back(side_bottom_start + seg);
            out_face_indices.push_back(side_bottom_start + seg + 1);
            out_face_indices.push_back(side_top_start + seg + 1);
        }

        return true;
    }

    return false;
}

static bool convert_shape_to_mesh(
    const pxr::UsdPrim& prim,
    pxr::UsdTimeCode time,
    pxr::VtArray<pxr::GfVec3f>& out_points,
    pxr::VtArray<int>& out_face_counts,
    pxr::VtArray<int>& out_face_indices,
    pxr::VtArray<pxr::GfVec3f>& out_normals)
{
    using namespace pxr;

    if (prim.IsA<UsdGeomMesh>()) {
        UsdGeomMesh mesh(prim);
        mesh.GetPointsAttr().Get(&out_points, time);
        mesh.GetFaceVertexCountsAttr().Get(&out_face_counts, time);
        mesh.GetFaceVertexIndicesAttr().Get(&out_face_indices, time);
        if (mesh.GetNormalsAttr()) {
            mesh.GetNormalsAttr().Get(&out_normals, time);
        }
        return true;
    }

    ParametricShapeParams params;
    std::string type_name;

    if (prim.IsA<UsdGeomSphere>()) {
        type_name = "Sphere";
        UsdGeomSphere sphere(prim);
        sphere.GetRadiusAttr().Get(&params.radius, time);
    }
    else if (prim.IsA<UsdGeomCylinder>()) {
        type_name = "Cylinder";
        UsdGeomCylinder cyl(prim);
        cyl.GetRadiusAttr().Get(&params.radius, time);
        cyl.GetHeightAttr().Get(&params.height, time);
    }
    else if (prim.IsA<UsdGeomCube>()) {
        type_name = "Cube";
        UsdGeomCube cube(prim);
        cube.GetSizeAttr().Get(&params.size, time);
    }
    else if (prim.IsA<UsdGeomCone>()) {
        type_name = "Cone";
        UsdGeomCone cone(prim);
        cone.GetRadiusAttr().Get(&params.radius, time);
        cone.GetHeightAttr().Get(&params.height, time);
    }
    else if (prim.IsA<UsdGeomCapsule>()) {
        type_name = "Capsule";
        UsdGeomCapsule capsule(prim);
        capsule.GetRadiusAttr().Get(&params.radius, time);
        capsule.GetHeightAttr().Get(&params.height, time);
    }
    else {
        return false;
    }

    return generate_parametric_mesh(
        type_name,
        params,
        out_points,
        out_face_counts,
        out_face_indices,
        out_normals);
}

// ============================================================================
// read_geometry_from_usd - 从 USD prim 读取几何数据到 Geometry
// ============================================================================
bool read_geometry_from_usd(
    Geometry& geometry,
    const pxr::UsdPrim& prim,
    pxr::UsdTimeCode time)
{
    using namespace pxr;

    if (!prim) {
        spdlog::error("[read_geometry_from_usd] Invalid prim");
        return false;
    }

    // Check if prim is a convertible geometry type
    if (!is_convertible_geom_shape(prim)) {
        spdlog::error(
            "[read_geometry_from_usd] Prim '{}' is not a convertible geometry "
            "type",
            prim.GetPath().GetString());
        return false;
    }

    // Convert shape to mesh data
    VtArray<GfVec3f> points;
    VtArray<int> face_counts;
    VtArray<int> face_indices;
    VtArray<GfVec3f> normals;

    if (!convert_shape_to_mesh(
            prim, time, points, face_counts, face_indices, normals)) {
        spdlog::error(
            "[read_geometry_from_usd] Failed to convert prim to mesh");
        return false;
    }

    // 创建或获取 MeshComponent
    auto mesh_comp = geometry.get_component<MeshComponent>();
    if (!mesh_comp) {
        mesh_comp = std::make_shared<MeshComponent>(&geometry);
        geometry.attach_component(mesh_comp);
    }

    auto mesh_view = get_usd_view(*mesh_comp);

    // Set vertices
    if (!points.empty()) {
        mesh_view.set_vertices(points);
    }

    // Set topology
    if (!face_counts.empty() && !face_indices.empty()) {
        mesh_view.set_face_topology(face_counts, face_indices);
    }

    // Set normals
    if (!normals.empty()) {
        mesh_view.set_normals(normals);
    }

    // 读取颜色 (for mesh types)
    VtArray<GfVec3f> colors;
    if (prim.IsA<UsdGeomMesh>()) {
        UsdGeomMesh usd_mesh(prim);
        if (usd_mesh.GetDisplayColorAttr()) {
            usd_mesh.GetDisplayColorAttr().Get(&colors, time);
            if (!colors.empty()) {
                mesh_view.set_display_colors(colors);
            }
        }

        // 读取 UV 坐标
        UsdGeomPrimvarsAPI primvar_api(usd_mesh);
        auto uv_primvar = primvar_api.GetPrimvar(TfToken("UVMap"));
        if (!uv_primvar) {
            uv_primvar = primvar_api.GetPrimvar(TfToken("st"));
        }
        if (uv_primvar) {
            VtArray<GfVec2f> uvs;
            uv_primvar.Get(&uvs, time);
            mesh_view.set_uv_coordinates(uvs);
        }
    }

    // 读取 Transform (作为 Geometry 的 XformComponent)
    UsdGeomXformable xformable(prim);
    if (xformable) {
        GfMatrix4d xform_matrix = xformable.ComputeLocalToWorldTransform(time);
        if (xform_matrix != GfMatrix4d().SetIdentity()) {
            auto xform_comp = geometry.get_component<XformComponent>();
            if (!xform_comp) {
                xform_comp = std::make_shared<XformComponent>(&geometry);
                geometry.attach_component(xform_comp);
            }

            auto translation = xform_matrix.ExtractTranslation();
            xform_comp->translation.clear();
            xform_comp->translation.push_back(
                glm::vec3(translation[0], translation[1], translation[2]));
            xform_comp->rotation.push_back(glm::vec3(0.0f));  // TODO: 提取旋转
            xform_comp->scale.push_back(glm::vec3(1.0f));     // TODO: 提取缩放
        }
    }

    // 读取 Skeleton 数据 (for mesh types)
    if (prim.IsA<UsdGeomMesh>()) {
        UsdGeomMesh usd_mesh(prim);
        UsdSkelBindingAPI binding = UsdSkelBindingAPI(usd_mesh);
        SdfPathVector targets;
        binding.GetSkeletonRel().GetTargets(&targets);

        if (targets.size() == 1) {
            auto stage = prim.GetStage();
            auto skel_prim = stage->GetPrimAtPath(targets[0]);
            UsdSkelSkeleton skeleton(skel_prim);

            if (skeleton) {
                UsdSkelCache skelCache;
                UsdSkelSkeletonQuery skelQuery =
                    skelCache.GetSkelQuery(skeleton);

                auto skel_component = geometry.get_component<SkelComponent>();
                if (!skel_component) {
                    skel_component = std::make_shared<SkelComponent>(&geometry);
                    geometry.attach_component(skel_component);
                }

                VtArray<GfMatrix4f> xforms;
                skelQuery.ComputeJointLocalTransforms(&xforms, time);

                skel_component->localTransforms = xforms;
                skel_component->jointOrder = skelQuery.GetJointOrder();
                skel_component->topology = skelQuery.GetTopology();

                VtArray<float> jointWeight;
                binding.GetJointWeightsAttr().Get(&jointWeight, time);

                VtArray<GfMatrix4d> bindTransforms;
                skeleton.GetBindTransformsAttr().Get(&bindTransforms, time);
                skel_component->bindTransforms = bindTransforms;

                VtArray<int> jointIndices;
                binding.GetJointIndicesAttr().Get(&jointIndices, time);
                skel_component->jointWeight = jointWeight;
                skel_component->jointIndices = jointIndices;
            }
        }
    }

    return true;
}

// Utility functions for ConstMeshUSDView
static pxr::VtArray<pxr::GfVec3f> vec3f_array_to_vt_array(
    const std::vector<glm::vec3>& array)
{
    return pxr::VtArray<pxr::GfVec3f>(
        reinterpret_cast<const pxr::GfVec3f*>(array.data()),
        reinterpret_cast<const pxr::GfVec3f*>(array.data() + array.size()));
}

static pxr::VtArray<pxr::GfVec2f> vec2f_array_to_vt_array(
    const std::vector<glm::vec2>& array)
{
    return pxr::VtArray<pxr::GfVec2f>(
        reinterpret_cast<const pxr::GfVec2f*>(array.data()),
        reinterpret_cast<const pxr::GfVec2f*>(array.data() + array.size()));
}

static pxr::VtArray<float> float_array_to_vt_array(
    const std::vector<float>& array)
{
    return pxr::VtArray<float>(array.begin(), array.end());
}

static pxr::VtArray<int> int_array_to_vt_array(const std::vector<int>& array)
{
    return pxr::VtArray<int>(array.begin(), array.end());
}

bool legal(const std::string& string)
{
    if (string.empty()) {
        return false;
    }
    if (std::find_if(string.begin(), string.end(), [](char val) {
            return val == '(' || val == ')' || val == ',';
        }) == string.end()) {
        return true;
    }
    return false;
}

bool write_geometry_to_usd(
    const Geometry& geometry,
    pxr::UsdStageRefPtr stage,
    const pxr::SdfPath& sdf_path,
    pxr::UsdTimeCode time)
{
    Geometry geom_copy = geometry;
    geom_copy.apply_transform();

    auto mesh = geom_copy.get_component<MeshComponent>();
    auto points = geom_copy.get_component<PointsComponent>();
    auto curve = geom_copy.get_component<CurveComponent>();
    auto volume = geom_copy.get_component<VolumeComponent>();
    auto instancer = geom_copy.get_component<InstancerComponent>();

    assert(!(points && mesh));

    pxr::SdfPath actual_path = sdf_path;
    if (instancer) {
        actual_path = sdf_path.AppendPath(pxr::SdfPath("Prototype"));
    }

    if (mesh) {
        auto mesh_usdview = get_usd_view(*mesh);
        pxr::UsdGeomMesh usdgeom = pxr::UsdGeomMesh::Define(stage, actual_path);

        if (usdgeom) {
            usdgeom.CreatePointsAttr().Set(mesh_usdview.get_vertices(), time);
            usdgeom.CreateFaceVertexCountsAttr().Set(
                mesh_usdview.get_face_vertex_counts(), time);
            usdgeom.CreateFaceVertexIndicesAttr().Set(
                mesh_usdview.get_face_vertex_indices(), time);

            auto primVarAPI = pxr::UsdGeomPrimvarsAPI(usdgeom);

            if (!mesh_usdview.get_normals().empty()) {
                usdgeom.CreateNormalsAttr().Set(
                    mesh_usdview.get_normals(), time);

                // Check if normals are per-vertex or per-face-vertex
                // (face-varying)
                size_t num_normals = mesh_usdview.get_normals().size();
                size_t num_vertices = mesh_usdview.get_vertices().size();
                size_t num_face_vertices =
                    mesh_usdview.get_face_vertex_indices().size();

                if (num_normals == num_vertices) {
                    // One normal per vertex - vertex interpolation
                    usdgeom.SetNormalsInterpolation(pxr::UsdGeomTokens->vertex);
                }
                else if (num_normals == num_face_vertices) {
                    // One normal per face-vertex - faceVarying
                    // interpolation
                    usdgeom.SetNormalsInterpolation(
                        pxr::UsdGeomTokens->faceVarying);
                }
                else {
                    // Mismatch - log error but set faceVarying as fallback
                    usdgeom.SetNormalsInterpolation(
                        pxr::UsdGeomTokens->faceVarying);
                }

                usdgeom.CreateSubdivisionSchemeAttr().Set(
                    pxr::UsdGeomTokens->none);
            }
            else {
                auto normals_attr = usdgeom.GetNormalsAttr();
                if (normals_attr) {
                    normals_attr.Block();
                }
            }

            if (!mesh_usdview.get_display_colors().empty()) {
                auto colorPrimvar = primVarAPI.CreatePrimvar(
                    pxr::TfToken("displayColor"),
                    pxr::SdfValueTypeNames->Color3fArray);

                // Determine interpolation based on color count
                size_t num_colors = mesh_usdview.get_display_colors().size();
                size_t num_vertices = mesh_usdview.get_vertices().size();
                size_t num_faces = mesh_usdview.get_face_vertex_counts().size();

                if (num_colors == num_vertices) {
                    colorPrimvar.SetInterpolation(pxr::UsdGeomTokens->vertex);
                }
                else if (num_colors == num_faces) {
                    colorPrimvar.SetInterpolation(pxr::UsdGeomTokens->uniform);
                }
                else {
                    // Default to vertex interpolation
                    colorPrimvar.SetInterpolation(pxr::UsdGeomTokens->vertex);
                }

                colorPrimvar.Set(mesh_usdview.get_display_colors(), time);
            }
            else {
                auto colorPrimvar =
                    primVarAPI.GetPrimvar(pxr::TfToken("displayColor"));
                if (colorPrimvar) {
                    colorPrimvar.GetAttr().Block();
                }
            }

            if (!mesh_usdview.get_uv_coordinates().empty()) {
                auto primvar = primVarAPI.CreatePrimvar(
                    pxr::TfToken("UVMap"),
                    pxr::SdfValueTypeNames->TexCoord2fArray);
                primvar.Set(mesh_usdview.get_uv_coordinates(), time);
                if (mesh_usdview.get_uv_coordinates().size() ==
                    mesh_usdview.get_vertices().size()) {
                    primvar.SetInterpolation(pxr::UsdGeomTokens->vertex);
                }
                else {
                    primvar.SetInterpolation(pxr::UsdGeomTokens->faceVarying);
                }
            }
            else {
                auto uvPrimvar = primVarAPI.GetPrimvar(pxr::TfToken("UVMap"));
                if (uvPrimvar) {
                    uvPrimvar.GetAttr().Block();
                }
            }

            usdgeom.CreateDoubleSidedAttr().Set(true);

            // Write all mesh attributes
            for (const std::string& name :
                 mesh->get_vertex_scalar_quantity_names()) {
                auto values = mesh_usdview.get_vertex_scalar_quantity(name);
                const std::string primvar_name = "vertex:scalar:" + name;
                auto primvar = primVarAPI.CreatePrimvar(
                    pxr::TfToken(primvar_name.c_str()),
                    pxr::SdfValueTypeNames->FloatArray);
                primvar.SetInterpolation(pxr::UsdGeomTokens->vertex);
                primvar.Set(values, time);
            }

            for (const std::string& name :
                 mesh->get_face_scalar_quantity_names()) {
                auto values = mesh_usdview.get_face_scalar_quantity(name);
                const std::string primvar_name = "face:scalar:" + name;
                auto primvar = primVarAPI.CreatePrimvar(
                    pxr::TfToken(primvar_name.c_str()),
                    pxr::SdfValueTypeNames->FloatArray);
                primvar.SetInterpolation(pxr::UsdGeomTokens->uniform);
                primvar.Set(values, time);
            }

            for (const std::string& name :
                 mesh->get_vertex_vector_quantity_names()) {
                auto values = mesh_usdview.get_vertex_vector_quantity(name);
                const std::string primvar_name = "vertex:vector:" + name;
                auto primvar = primVarAPI.CreatePrimvar(
                    pxr::TfToken(primvar_name.c_str()),
                    pxr::SdfValueTypeNames->Vector3fArray);
                primvar.SetInterpolation(pxr::UsdGeomTokens->vertex);
                primvar.Set(values, time);
            }

            for (const std::string& name :
                 mesh->get_face_vector_quantity_names()) {
                auto values = mesh_usdview.get_face_vector_quantity(name);
                const std::string primvar_name = "face:vector:" + name;
                auto primvar = primVarAPI.CreatePrimvar(
                    pxr::TfToken(primvar_name.c_str()),
                    pxr::SdfValueTypeNames->Vector3fArray);
                primvar.SetInterpolation(pxr::UsdGeomTokens->uniform);
                primvar.Set(values, time);
            }

            for (const std::string& name :
                 mesh->get_face_corner_parameterization_quantity_names()) {
                auto values =
                    mesh_usdview.get_face_corner_parameterization_quantity(
                        name);
                const std::string primvar_name =
                    "face_corner:parameterization:" + name;
                auto primvar = primVarAPI.CreatePrimvar(
                    pxr::TfToken(primvar_name.c_str()),
                    pxr::SdfValueTypeNames->TexCoord2fArray);
                primvar.SetInterpolation(pxr::UsdGeomTokens->faceVarying);
                primvar.Set(values, time);
            }

            for (const std::string& name :
                 mesh->get_vertex_parameterization_quantity_names()) {
                auto values =
                    mesh_usdview.get_vertex_parameterization_quantity(name);
                const std::string primvar_name =
                    "vertex:parameterization:" + name;
                auto primvar = primVarAPI.CreatePrimvar(
                    pxr::TfToken(primvar_name.c_str()),
                    pxr::SdfValueTypeNames->TexCoord2fArray);
                primvar.SetInterpolation(pxr::UsdGeomTokens->vertex);
                primvar.Set(values, time);
            }
        }
    }
    else if (points) {
        pxr::UsdGeomPoints usdpoints =
            pxr::UsdGeomPoints::Define(stage, actual_path);

        usdpoints.CreatePointsAttr().Set(
            vec3f_array_to_vt_array(points->get_vertices()), time);

        if (points->get_width().size() > 0) {
            usdpoints.CreateWidthsAttr().Set(
                float_array_to_vt_array(points->get_width()), time);
        }

        auto PrimVarAPI = pxr::UsdGeomPrimvarsAPI(usdpoints);
        if (points->get_display_color().size() > 0) {
            pxr::UsdGeomPrimvar colorPrimvar = PrimVarAPI.CreatePrimvar(
                pxr::TfToken("displayColor"),
                pxr::SdfValueTypeNames->Color3fArray);
            colorPrimvar.SetInterpolation(pxr::UsdGeomTokens->vertex);
            colorPrimvar.Set(
                vec3f_array_to_vt_array(points->get_display_color()), time);
        }
    }
    else if (curve) {
        pxr::UsdGeomBasisCurves usd_curve =
            pxr::UsdGeomBasisCurves::Define(stage, actual_path);
        if (usd_curve) {
            // Set curve type and basis first
            usd_curve.CreateTypeAttr().Set(
                curve->get_type() == CurveComponent::CurveType::Linear
                    ? pxr::UsdGeomTokens->linear
                    : pxr::UsdGeomTokens->cubic,
                time);

            // Set basis for cubic curves
            if (curve->get_type() != CurveComponent::CurveType::Linear) {
                usd_curve.CreateBasisAttr().Set(
                    pxr::UsdGeomTokens->bspline, time);
            }

            // Set wrap mode
            if (curve->get_periodic())
                usd_curve.CreateWrapAttr().Set(pxr::UsdGeomTokens->periodic);
            else
                usd_curve.CreateWrapAttr().Set(pxr::UsdGeomTokens->nonperiodic);

            // Set geometry data
            usd_curve.CreatePointsAttr().Set(
                vec3f_array_to_vt_array(curve->get_vertices()), time);
            usd_curve.CreateWidthsAttr().Set(
                float_array_to_vt_array(curve->get_width()), time);
            usd_curve.CreateCurveVertexCountsAttr().Set(
                int_array_to_vt_array(curve->get_vert_count()), time);
            usd_curve.CreateNormalsAttr().Set(
                vec3f_array_to_vt_array(curve->get_curve_normals()), time);
            usd_curve.CreateDisplayColorAttr().Set(
                vec3f_array_to_vt_array(curve->get_display_color()), time);
        }
    }
    else if (volume) {
        pxr::UsdVolVolume usd_volume =
            pxr::UsdVolVolume::Define(stage, actual_path);

        if (!usd_volume) {
            return false;
        }

        auto vdb_asset_path = actual_path.AppendChild(pxr::TfToken("field"));
        auto openvdb_asset =
            pxr::UsdVolOpenVDBAsset::Define(stage, vdb_asset_path);

        if (!openvdb_asset) {
            return false;
        }

        auto file_name = "volume" + actual_path.GetName() +
                         std::to_string(time.GetValue()) + ".vdb";
        volume->write_disk(file_name);

        openvdb_asset.CreateFilePathAttr().Set(
            pxr::SdfAssetPath(file_name), time);

        if (!usd_volume.CreateFieldRelationship(
                pxr::TfToken("field"), vdb_asset_path)) {
            return false;
        }

        usd_volume.MakeVisible(time);
    }
    else {
        return false;
    }

    // Handle instancer
    if (instancer) {
        auto instancer_component =
            pxr::UsdGeomPointInstancer::Define(stage, sdf_path);
        instancer_component.CreatePrototypesRel().SetTargets({ actual_path });

        const auto& positions_glm = instancer->get_positions();
        const auto& orientations_glm = instancer->get_orientations();
        const auto& scales_glm = instancer->get_scales();

        size_t instance_count = instancer->get_instance_count();

        // Convert positions directly
        pxr::VtVec3fArray positions(instance_count);
        for (size_t i = 0; i < instance_count; ++i) {
            positions[i] = pxr::GfVec3f(
                positions_glm[i].x, positions_glm[i].y, positions_glm[i].z);
        }

        // Handle orientations if rotations are enabled
        if (instancer->has_rotations_enabled()) {
            // Only write orientations if array is not empty
            if (!orientations_glm.empty()) {
                pxr::VtQuathArray orientations(instance_count);
                for (size_t i = 0; i < instance_count; ++i) {
                    const auto& q = orientations_glm[i];
                    orientations[i] = pxr::GfQuath(q.w, q.x, q.y, q.z);
                }
                instancer_component.CreateOrientationsAttr().Set(
                    orientations, time);
            }
            else {
                // No orientations data, clear the attribute
                auto orientations_attr =
                    instancer_component.GetOrientationsAttr();
                if (orientations_attr) {
                    orientations_attr.Block();
                }
            }

            // Only write scales if array is not empty
            if (!scales_glm.empty()) {
                pxr::VtVec3fArray scales(instance_count);
                for (size_t i = 0; i < instance_count; ++i) {
                    scales[i] = pxr::GfVec3f(
                        scales_glm[i].x, scales_glm[i].y, scales_glm[i].z);
                }
                instancer_component.CreateScalesAttr().Set(scales, time);
            }
            else {
                // No scales data, clear the attribute
                auto scales_attr = instancer_component.GetScalesAttr();
                if (scales_attr) {
                    scales_attr.Block();
                }
            }
        }
        else {
            // Clear orientation attribute when rotations are disabled
            auto orientations_attr = instancer_component.GetOrientationsAttr();
            if (orientations_attr) {
                orientations_attr.Block();
            }
        }

        instancer_component.CreateProtoIndicesAttr().Set(
            pxr::VtIntArray(instancer->get_proto_indices()), time);
        instancer_component.CreatePositionsAttr().Set(positions, time);
    }

    // Handle materials
    auto material_component = geometry.get_component<MaterialComponent>();
    if (material_component) {
        auto usdgeom = pxr::UsdGeomXformable::Get(stage, actual_path);
        if (legal(std::string(material_component->textures[0].c_str()))) {
            auto material_path = material_component->get_material_path();
            auto material =
                material_component->define_material(stage, material_path);
            usdgeom.GetPrim().ApplyAPI(pxr::UsdShadeTokens->MaterialBindingAPI);
            pxr::UsdShadeMaterialBindingAPI(usdgeom).Bind(material);
        }
    }

    pxr::UsdGeomImageable(stage->GetPrimAtPath(actual_path)).MakeVisible();
    return true;
}

// ============================================================================
// Modifier Mode Implementation
// ============================================================================

pxr::SdfPath get_modifier_output_path(
    const pxr::SdfPath& prim_path,
    int modifier_index)
{
    std::string modifier_name = "modifier_" + std::to_string(modifier_index);
    return prim_path.AppendPath(pxr::SdfPath("modifiers"))
        .AppendPath(pxr::SdfPath(modifier_name));
}

bool write_geometry_as_over_spec(
    const Geometry& geometry,
    pxr::UsdStageRefPtr stage,
    const pxr::SdfPath& sdf_path,
    pxr::UsdTimeCode time,
    pxr::SdfLayerHandle modifier_layer)
{
    using namespace pxr;

    if (!modifier_layer) {
        spdlog::error("[write_geometry_as_over_spec] Invalid modifier layer");
        return false;
    }

    if (!stage) {
        spdlog::error("[write_geometry_as_over_spec] Invalid stage");
        return false;
    }

    if (sdf_path.IsEmpty()) {
        spdlog::error("[write_geometry_as_over_spec] Empty sdf_path");
        return false;
    }

    Geometry geom_copy = geometry;
    geom_copy.apply_transform();

    auto mesh = geom_copy.get_component<MeshComponent>();
    if (!mesh) {
        spdlog::error("[write_geometry_as_over_spec] No mesh component found");
        return false;
    }

    auto mesh_usdview = get_usd_view(*mesh);

    spdlog::debug(
        "[MODIFIER] write_geometry_as_over_spec: Creating over spec at "
        "'{}' in layer '{}'",
        sdf_path.GetString(),
        modifier_layer->GetIdentifier());

    // Create the prim spec in the modifier layer (not root layer!)
    // Use SdfJustCreatePrimInLayer to create prim hierarchy
    if (!SdfJustCreatePrimInLayer(modifier_layer, sdf_path)) {
        spdlog::error(
            "[write_geometry_as_over_spec] Failed to create prim hierarchy at: "
            "{}",
            sdf_path.GetString());
        return false;
    }

    // Now get the prim spec
    SdfPrimSpecHandle prim_spec = modifier_layer->GetPrimAtPath(sdf_path);
    if (!prim_spec) {
        spdlog::error(
            "[write_geometry_as_over_spec] Failed to get prim spec at: {}",
            sdf_path.GetString());
        return false;
    }

    spdlog::debug(
        "[MODIFIER] Prim spec created, vertices count: {}, faces count: "
        "{}",
        mesh_usdview.get_vertices().size(),
        mesh_usdview.get_face_vertex_counts().size());

    // Set specifier to "over" (not "def") - this is the key for
    // non-destructive editing
    prim_spec->SetSpecifier(SdfSpecifierOver);

    // Set the type name
    prim_spec->SetTypeName(TfToken("Mesh"));

    // Helper to set or update an attribute value
    // Important: We must check if attribute already exists before creating
    auto set_attribute_value = [&prim_spec, &modifier_layer, &time, &sdf_path](
                                   const TfToken& attr_name,
                                   const SdfValueTypeName& type_name,
                                   const VtValue& value) -> bool {
        // Build correct attribute path: /prim_path.attr_name
        SdfPath attr_path = sdf_path.AppendProperty(attr_name);

        // Check if attribute spec already exists in the layer
        SdfAttributeSpecHandle attr_spec =
            modifier_layer->GetAttributeAtPath(attr_path);

        if (!attr_spec) {
            // Create new attribute spec only if it doesn't exist
            attr_spec = SdfAttributeSpec::New(
                prim_spec, attr_name, type_name, SdfVariabilityVarying);
            if (!attr_spec) {
                spdlog::debug(
                    "[MODIFIER] Failed to create attribute: {}",
                    attr_name.GetString());
                return false;
            }
        }

        // Set the value (this updates existing spec)
        attr_spec->SetDefaultValue(value);

        // Also set time sample if not default time
        if (time != UsdTimeCode::Default()) {
            modifier_layer->SetTimeSample(
                attr_spec->GetPath(), time.GetValue(), value);
        }

        return true;
    };

    // Write points
    if (!mesh_usdview.get_vertices().empty()) {
        set_attribute_value(
            TfToken("points"),
            SdfValueTypeNames->Point3fArray,
            VtValue(mesh_usdview.get_vertices()));
    }

    // Write face topology
    if (!mesh_usdview.get_face_vertex_counts().empty()) {
        set_attribute_value(
            TfToken("faceVertexCounts"),
            SdfValueTypeNames->IntArray,
            VtValue(mesh_usdview.get_face_vertex_counts()));
    }

    if (!mesh_usdview.get_face_vertex_indices().empty()) {
        set_attribute_value(
            TfToken("faceVertexIndices"),
            SdfValueTypeNames->IntArray,
            VtValue(mesh_usdview.get_face_vertex_indices()));
    }

    // Write normals
    if (!mesh_usdview.get_normals().empty()) {
        set_attribute_value(
            TfToken("normals"),
            SdfValueTypeNames->Vector3fArray,
            VtValue(mesh_usdview.get_normals()));
    }

    // Note: primvars (UV, colors) are not written here because SdfAttributeSpec
    // cannot handle namespaced names like "primvars:UVMap" directly.
    // For full primvar support, we would need to use SdfPrimvarSpec API.
    // For now, the over spec only overrides basic geometry (points, topology,
    // normals). Primvars from the original prim will still be visible through
    // composition.

    spdlog::debug(
        "[MODIFIER] Successfully wrote over spec to '{}', points: {}, "
        "faces: {}",
        sdf_path.GetString(),
        mesh_usdview.get_vertices().size(),
        mesh_usdview.get_face_vertex_counts().size());

    // Debug: print first few points
    const auto& verts = mesh_usdview.get_vertices();
    if (!verts.empty()) {
        spdlog::debug(
            "[MODIFIER] First 3 points: ({}, {}, {}), ({}, {}, {}), ({}, "
            "{}, {})",
            verts[0][0],
            verts[0][1],
            verts[0][2],
            verts.size() > 1 ? verts[1][0] : 0,
            verts.size() > 1 ? verts[1][1] : 0,
            verts.size() > 1 ? verts[1][2] : 0,
            verts.size() > 2 ? verts[2][0] : 0,
            verts.size() > 2 ? verts[2][1] : 0,
            verts.size() > 2 ? verts[2][2] : 0);
    }

    return true;
}

RUZINO_NAMESPACE_CLOSE_SCOPE
