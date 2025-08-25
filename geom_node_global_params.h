#pragma once

#include "GCore/api.h"


USTC_CG_NAMESPACE_OPEN_SCOPE
struct PickEvent;
struct GeomNodeGlobalParams {
    pxr::SdfPath prim_path;
    PickEvent* pick;
};
USTC_CG_NAMESPACE_CLOSE_SCOPE