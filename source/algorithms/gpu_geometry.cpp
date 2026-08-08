#include <GCore/algorithms/gpu_geometry.h>

#ifdef GPU_GEOM_ALGORITHM

#include <RHI/rhi.hpp>
#include <cstdlib>
#include <spdlog/spdlog.h>

#include <memory>

#include "RHI/shaderCompiler.h"
#include "RHI/ShaderFactory/shader.hpp"

RUZINO_NAMESPACE_OPEN_SCOPE

static ResourceAllocator resource_allocator_;
static std::shared_ptr<ShaderFactory> shader_factory;
static bool gpu_alive_ = false;
static bool atexit_registered_ = false;

ResourceAllocator& get_resource_allocator()
{
    if (gpu_alive_)
        return resource_allocator_;
    init_gpu_geometry_algorithms();
    return resource_allocator_;
}

bool is_gpu_alive() { return gpu_alive_; }

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
            .string());
    shader_factory->add_search_path(
        SlangShaderCompiler::get_shader_dir(ShaderDirType::GeomNodes)
            .string());
    shader_factory->add_search_path(
        SlangShaderCompiler::get_shader_dir(ShaderDirType::GeomCompute)
            .string());
    resource_allocator_.shader_factory = shader_factory.get();
    gpu_alive_ = true;

    // Register cleanup via atexit so it runs before static destructors.
    // This ensures all GPU resources are returned to the allocator cache,
    // then the allocator cache is cleared, then RHI::shutdown() releases
    // the device — all while the relevant DLLs are still loaded.
    if (!atexit_registered_) {
        atexit_registered_ = true;
        std::atexit(deinit_gpu_geometry_algorithms);
    }
}

void deinit_gpu_geometry_algorithms()
{
    if (!gpu_alive_)
        return;

    gpu_alive_ = false;

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
