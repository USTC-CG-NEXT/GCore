#pragma once

#include <pxr/usd/sdf/layer.h>
#include <pxr/usd/usd/stage.h>

#include "GCore/GOP.h"
#include "GCore/api.h"

RUZINO_NAMESPACE_OPEN_SCOPE class Stage;

void GEOMETRY_API init(Stage* stage);

// 从 USD prim 读取几何数据到 Geometry 对象
bool GEOMETRY_API read_geometry_from_usd(
    Geometry& geometry,
    const pxr::UsdPrim& prim,
    pxr::UsdTimeCode time);

// 从 SdfLayer 直接读取几何数据（跳过 composition，用于 modifier mode）
bool GEOMETRY_API read_geometry_from_layer(
    Geometry& geometry,
    pxr::SdfLayerHandle layer,
    const pxr::SdfPath& prim_path,
    pxr::UsdTimeCode time);

// 将 Geometry 对象写入到 USD stage
bool GEOMETRY_API write_geometry_to_usd(
    const Geometry& geometry,
    pxr::UsdStageRefPtr stage,
    const pxr::SdfPath& sdf_path,
    pxr::UsdTimeCode time);

// ============================================================================
// Modifier Mode Functions (Non-destructive layer composition)
// ============================================================================

// 将 Geometry 作为 over spec 写入（非破坏性修改器模式）
// modifier_layer: 指定要写入的 layer（通常为 session layer 或专用 modifier
// layer）
bool GEOMETRY_API write_geometry_as_over_spec(
    const Geometry& geometry,
    pxr::UsdStageRefPtr stage,
    const pxr::SdfPath& sdf_path,
    pxr::UsdTimeCode time,
    pxr::SdfLayerHandle modifier_layer);

// 获取修改器输出路径
// 格式: /path/to/prim/modifiers/modifier_N
pxr::SdfPath GEOMETRY_API
get_modifier_output_path(const pxr::SdfPath& prim_path, int modifier_index);

// ============================================================================
// Parametric Shape Mesh Generation
// ============================================================================

struct ParametricShapeParams {
    double radius = 1.0;
    double height = 2.0;
    double size = 2.0;
};

bool GEOMETRY_API generate_parametric_mesh(
    const std::string& type_name,
    const ParametricShapeParams& params,
    pxr::VtArray<pxr::GfVec3f>& out_points,
    pxr::VtArray<int>& out_face_counts,
    pxr::VtArray<int>& out_face_indices,
    pxr::VtArray<pxr::GfVec3f>& out_normals);

RUZINO_NAMESPACE_CLOSE_SCOPE
