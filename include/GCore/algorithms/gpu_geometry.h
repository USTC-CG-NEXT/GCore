#pragma once

#include "GCore/api.h"

#ifdef GPU_GEOM_ALGORITHM
#include "RHI/ResourceManager/resource_allocator.hpp"
#endif

RUZINO_NAMESPACE_OPEN_SCOPE

#ifdef GPU_GEOM_ALGORITHM

// Initialize GPU geometry algorithms: set up the resource allocator,
// shader factory, and hold an RHI reference to keep the device alive.
GEOMETRY_API void init_gpu_geometry_algorithms();

// Tear down GPU geometry algorithms: release all cached GPU resources,
// then release the RHI reference so the device can be destroyed.
GEOMETRY_API void deinit_gpu_geometry_algorithms();

// Returns true between init and deinit — safe to use the allocator.
GEOMETRY_API bool is_gpu_alive();

GEOMETRY_API ResourceAllocator& get_resource_allocator();

#endif

RUZINO_NAMESPACE_CLOSE_SCOPE
