#pragma once

#include <pxr/usd/usdShade/material.h>

#include <string>

#include "GCore/Components.h"
#include "GCore/GOP.h"
#include "GCore/MaterialObject.h"

RUZINO_NAMESPACE_OPEN_SCOPE
struct GEOMETRY_API MaterialComponent : public GeometryComponent {
    explicit MaterialComponent(Geometry* attached_operand)
        : GeometryComponent(attached_operand)
    {
    }

    void apply_transform(const glm::mat4& transform) override
    {
    }

    size_t hash() const override
    {
        size_t h = 0;
        for (const auto& texture : material_object.textures) {
            h ^= std::hash<std::string>{}(texture) + 0x9e3779b9 + (h << 6) +
                 (h >> 2);
        }
        h ^= std::hash<bool>{}(material_object.textured) + 0x9e3779b9 +
             (h << 6) + (h >> 2);
        auto color_hash = std::hash<float>{};
        h ^= (color_hash(material_object.base_color.r) ^
              (color_hash(material_object.base_color.g) << 1) ^
              (color_hash(material_object.base_color.b) << 2)) +
             0x9e3779b9 + (h << 6) + (h >> 2);
        h ^= std::hash<bool>{}(material_object.alpha_cutout) + 0x9e3779b9 +
             (h << 6) + (h >> 2);
        h ^= std::hash<float>{}(material_object.opacity_threshold) +
             0x9e3779b9 + (h << 6) + (h >> 2);
        h ^= std::hash<float>{}(material_object.roughness) + 0x9e3779b9 +
             (h << 6) + (h >> 2);
        h ^= std::hash<std::string>{}(material_object.wrap_mode) + 0x9e3779b9 +
             (h << 6) + (h >> 2);
        return h;
    }

    GeometryComponentHandle copy(Geometry* operand) const override
    {
        auto ret = std::make_shared<MaterialComponent>(operand);
        ret->material_object = this->material_object;
        return ret;
    }

    std::string to_string() const override
    {
        return {};
    }

    pxr::UsdShadeMaterial define_material(
        pxr::UsdStageRefPtr stage,
        pxr::SdfPath path) const;

    void set_materialx_path(pxr::SdfPath path)
    {
        mtlx_material_path = path;
    }

    pxr::SdfPath get_material_path() const;

    // The authored material data. The write path (usd_extension) materializes
    // deferred in-graph textures into image files (filling
    // material_object.textures) and authors the USD material network.
    MaterialObject material_object;

    pxr::SdfPath mtlx_material_path;
};

RUZINO_NAMESPACE_CLOSE_SCOPE
