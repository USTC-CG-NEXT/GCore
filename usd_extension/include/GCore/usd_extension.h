#pragma once

#include "GCore/GOP.h"
#include "GCore/api.h"

#include <pxr/usd/usd/stage.h>

USTC_CG_NAMESPACE_OPEN_SCOPE class Stage;

void GEOMETRY_API init(Stage* stage);

bool GEOMETRY_API write_geometry_to_usd(
    const Geometry& geometry,
    pxr::UsdStageRefPtr stage,
    const pxr::SdfPath& sdf_path,
    pxr::UsdTimeCode time);

USTC_CG_NAMESPACE_CLOSE_SCOPE
