#pragma once

#include <GCore/Texture/Texture.h>
#include <GCore/api.h>

#include <cstdint>
#include <memory>

RUZINO_NAMESPACE_OPEN_SCOPE

// Aggregate texture socket payload.
//
// Carries a texture that was computed in the node graph without forcing it
// to disk (storage is deferred to USD-write time) and without committing to
// CPU or GPU residency up front:
//   - cpu: an evaluable Texture (DataTexture2D pixels, FunctionTexture2D
//     procedural closure, ...). Materializable: write paths bake it to an
//     image file next to the stage when authoring the USD material.
//   - gpu_handle: opaque renderer-side texture handle (e.g. descriptor /
//     bindless id). Interpretation is deferred to the consumer; geometry
//     code never dereferences it.
// Both optionals may be set (GPU consumers prefer the handle, CPU
// consumers prefer the pixels). An empty TextureObject means "no texture".
struct TextureObject {
    std::shared_ptr<Texture> cpu;
    std::uint64_t gpu_handle = 0;

    bool has_cpu() const
    {
        return cpu != nullptr;
    }
    bool has_gpu() const
    {
        return gpu_handle != 0;
    }
    bool empty() const
    {
        return !has_cpu() && !has_gpu();
    }
};

RUZINO_NAMESPACE_CLOSE_SCOPE
