#include <GCore/usd_extension.h>
#include <pxr/base/gf/matrix4d.h>
#include <pxr/usd/usdGeom/curves.h>
#include <pxr/usd/usdGeom/mesh.h>
#include <pxr/usd/usdGeom/pointInstancer.h>
#include <pxr/usd/usdGeom/points.h>
#include <pxr/usd/usdGeom/primvarsAPI.h>
#include <pxr/usd/usdShade/materialBindingAPI.h>
#include <pxr/usd/usdVol/openVDBAsset.h>
#include <pxr/usd/usdVol/volume.h>

#include "GCore/Components/CurveComponent.h"
#include "GCore/Components/InstancerComponent.h"
#include "GCore/Components/MaterialComponent.h"
#include "GCore/Components/MeshComponent.h"
#include "GCore/Components/PointsComponent.h"
#include "GCore/Components/VolumeComponent.h"
#include "GCore/Components/XformComponent.h"
#include "glm/gtc/type_ptr.hpp"
#include "glm/gtx/matrix_decompose.hpp"

USTC_CG_NAMESPACE_OPEN_SCOPE

// Utility functions for ConstMeshUSDView
static pxr::VtArray<pxr::GfVec3f> vec3f_array_to_vt_array(
    const std::vector<glm::vec3>& array)
{
    pxr::VtArray<pxr::GfVec3f> vt_array(array.size());
    for (size_t i = 0; i < array.size(); ++i) {
        vt_array[i] = pxr::GfVec3f(array[i].x, array[i].y, array[i].z);
    }
    return vt_array;
}

static pxr::VtArray<pxr::GfVec2f> vec2f_array_to_vt_array(
    const std::vector<glm::vec2>& array)
{
    pxr::VtArray<pxr::GfVec2f> vt_array(array.size());
    for (size_t i = 0; i < array.size(); ++i) {
        vt_array[i] = pxr::GfVec2f(array[i].x, array[i].y);
    }
    return vt_array;
}
static pxr::VtArray<float> float_array_to_vt_array(
    const std::vector<float>& array)
{
    pxr::VtArray<float> vt_array(array.size());
    std::memcpy(vt_array.data(), array.data(), array.size() * sizeof(float));
    return vt_array;
}

// int
static pxr::VtArray<int> int_array_to_vt_array(const std::vector<int>& array)
{
    pxr::VtArray<int> vt_array(array.size());
    std::memcpy(vt_array.data(), array.data(), array.size() * sizeof(int));
    return vt_array;
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
    auto mesh = geometry.get_component<MeshComponent>();
    auto points = geometry.get_component<PointsComponent>();
    auto curve = geometry.get_component<CurveComponent>();
    auto volume = geometry.get_component<VolumeComponent>();
    auto instancer = geometry.get_component<InstancerComponent>();

    assert(!(points && mesh));

    pxr::SdfPath actual_path = sdf_path;
    if (instancer) {
        actual_path = sdf_path.AppendPath(pxr::SdfPath("Prototype"));
    }

    if (mesh) {
        auto mesh_usdview = mesh->get_usd_view();
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
                    // One normal per face-vertex - faceVarying interpolation
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
                usdgeom.CreateNormalsAttr().Block();
            }
            if (!mesh_usdview.get_display_colors().empty()) {
                auto colorPrimvar = primVarAPI.CreatePrimvar(
                    pxr::TfToken("displayColor"),
                    pxr::SdfValueTypeNames->Color3fArray);
                colorPrimvar.SetInterpolation(pxr::UsdGeomTokens->vertex);
                colorPrimvar.Set(mesh_usdview.get_display_colors(), time);
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

            for (std::string& name : mesh->get_vertex_color_quantity_names()) {
                auto values = mesh_usdview.get_vertex_color_quantity(name);
                const std::string primvar_name = "vertex:color:" + name;
                auto primvar = primVarAPI.CreatePrimvar(
                    pxr::TfToken(primvar_name.c_str()),
                    pxr::SdfValueTypeNames->Color3fArray);
                primvar.SetInterpolation(pxr::UsdGeomTokens->vertex);
                primvar.Set(values, time);
            }

            for (const std::string& name :
                 mesh->get_face_color_quantity_names()) {
                auto values = mesh_usdview.get_face_color_quantity(name);
                const std::string primvar_name = "face:color:" + name;
                auto primvar = primVarAPI.CreatePrimvar(
                    pxr::TfToken(primvar_name.c_str()),
                    pxr::SdfValueTypeNames->Color3fArray);
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

        auto transforms = instancer->get_instances();

        pxr::VtVec3fArray positions = pxr::VtVec3fArray(transforms.size());

        // Fast path: if no rotations, directly extract positions
        if (!instancer->has_rotations_enabled()) {
            for (size_t i = 0; i < transforms.size(); ++i) {
                // Directly extract translation from the last column of the
                // matrix
                positions[i] = pxr::GfVec3f(
                    transforms[i][3][0],
                    transforms[i][3][1],
                    transforms[i][3][2]);
            }
            // Clear orientation attribute when rotations are disabled
            auto orientations_attr = instancer_component.GetOrientationsAttr();
            if (orientations_attr) {
                orientations_attr.Block();
            }
        }
        else {
            // Full path: decompose matrices for positions, orientations, and
            // scales
            pxr::VtQuathArray orientations =
                pxr::VtQuathArray(transforms.size());
            pxr::VtVec3fArray scales = pxr::VtVec3fArray(transforms.size());

            for (size_t i = 0; i < transforms.size(); ++i) {
                glm::vec3 translation;
                glm::quat rotation;
                glm::vec3 scale;

                // Decompose GLM matrix
                glm::vec3 skew;
                glm::vec4 perspective;
                glm::decompose(
                    transforms[i],
                    scale,
                    rotation,
                    translation,
                    skew,
                    perspective);

                positions[i] =
                    pxr::GfVec3f(translation.x, translation.y, translation.z);
                orientations[i] = pxr::GfQuath(
                    rotation.w, rotation.x, rotation.y, rotation.z);
                scales[i] = pxr::GfVec3f(scale.x, scale.y, scale.z);
            }
            instancer_component.CreateOrientationsAttr().Set(
                orientations, time);
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

    // Handle transforms
    auto xform_component = geometry.get_component<XformComponent>();
    auto usdgeom = pxr::UsdGeomXformable::Get(stage, actual_path);
    auto xform_op = usdgeom.GetTransformOp();
    if (!xform_op) {
        xform_op = usdgeom.AddTransformOp();
    }

    if (xform_component) {
        assert(
            xform_component->translation.size() ==
            xform_component->rotation.size());
        glm::mat4 final_transform = xform_component->get_transform();
        pxr::GfMatrix4d usd_transform;
        const float* src = glm::value_ptr(final_transform);
        double* dst = usd_transform.GetArray();
        for (int i = 0; i < 16; ++i) {
            dst[i] = static_cast<double>(src[i]);
        }
        xform_op.Set(usd_transform, time);
    }
    else {
        xform_op.Set(pxr::GfMatrix4d().SetIdentity(), time);
    }

    pxr::UsdGeomImageable(stage->GetPrimAtPath(actual_path)).MakeVisible();
    return true;
}

USTC_CG_NAMESPACE_CLOSE_SCOPE
