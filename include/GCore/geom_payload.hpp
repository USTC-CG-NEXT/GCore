#pragma once

#include <glm/glm.hpp>
#include <memory>
#include <string>
#include <vector>
#ifdef GEOM_USD_EXTENSION
#include <pxr/usd/sdf/layer.h>
#include <pxr/usd/sdf/path.h>
#include <pxr/usd/usd/common.h>
#endif

struct PickEvent {
    glm::vec3 point;
    glm::vec3 normal;
#ifdef GEOM_USD_EXTENSION
    pxr::SdfPath prim_path;
    pxr::SdfPath instancer_path;
#endif

    PickEvent(
        const glm::vec3& pt,
        const glm::vec3& norm
#ifdef GEOM_USD_EXTENSION
        ,
        const pxr::SdfPath& path,
        const pxr::SdfPath& inst
#endif
        )
        : point(pt),
          normal(norm)
#ifdef GEOM_USD_EXTENSION
          ,
          prim_path(path),
          instancer_path(inst)
#endif

    {
    }
};

struct GeomPayload {
#ifdef GEOM_USD_EXTENSION
    pxr::UsdStageRefPtr stage;
    pxr::UsdTimeCode current_time = pxr::UsdTimeCode(0);
    pxr::SdfPath prim_path;

    // Modifier mode support (non-destructive layer composition)
    bool is_modifier_mode = false;
    pxr::SdfLayerHandle modifier_layer;
    pxr::SdfPath modifier_output_path;
    pxr::SdfPath modifier_input_path;  // Where input_geometry should read from
    int current_modifier_index =
        -1;  // Current modifier in stack (-1 = not in modifier mode)
#endif

    std::shared_ptr<PickEvent> pick;
    float delta_time = 0.0f;
    bool has_simulation = false;
    bool is_simulating = false;

    std::string stage_filepath_;
};
