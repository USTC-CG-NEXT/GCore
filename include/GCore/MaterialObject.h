#pragma once

#include <GCore/Texture/TextureObject.h>
#include <GCore/api.h>

#include <glm/glm.hpp>
#include <string>
#include <vector>

RUZINO_NAMESPACE_OPEN_SCOPE

// Socket payload produced by create_material and attached by set_material.
// Pure data — no geometry attachment, no USD authoring: the write path
// (usd_extension) materializes deferred textures and authors the
// UsdPreviewSurface network from these fields.
struct MaterialObject {
    // Texture sources. Priority: texture_object (in-graph computed, storage
    // deferred to write time) > textures[0] (file path). With textured=false
    // both are ignored and the material is flat base_color.
    TextureObject texture_object;
    std::vector<std::string> textures;
    bool textured = true;
    glm::vec3 base_color{ 0.5f, 0.5f, 0.5f };

    // Alpha-cutout (mask) transparency: binds the texture's alpha channel as
    // UsdPreviewSurface opacity with opacityMode=mask.
    bool alpha_cutout = false;
    float opacity_threshold = 0.5f;
    float roughness = 0.5f;
    // UsdUVTexture wrapS/wrapT token: mirror | clamp | periodic | black.
    std::string wrap_mode = "clamp";

    bool has_texture_source() const
    {
        return (textured && !textures.empty()) || texture_object.has_cpu() ||
               texture_object.has_gpu();
    }
};

RUZINO_NAMESPACE_CLOSE_SCOPE
