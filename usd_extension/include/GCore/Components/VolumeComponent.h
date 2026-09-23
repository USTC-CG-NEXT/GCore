#pragma once

#include <functional>

#include "GCore/Components.h"
#include "GCore/GOP.h"
#include "GCore/api.h"
#include "GCore/geom_payload.hpp"
#include "openvdb/openvdb.h"

RUZINO_NAMESPACE_OPEN_SCOPE
class GEOMETRY_API VolumeComponent : public GeometryComponent {
   public:
    explicit VolumeComponent(Geometry* attached_operand)
        : GeometryComponent(attached_operand)
    {
    }

    GeometryComponentHandle copy(Geometry* operand) const override;
    std::string to_string() const override;

    size_t hash() const override
    {
        if (!grid) {
            return 0;
        }

        // Hashing every voxel would be far too expensive, and the interior of
        // a fog volume is fully active, so the grid is summarised by a few
        // cheap invariants plus a bounded sample of its values.
        size_t h = 0;
        const auto mix = [&h](size_t value) {
            h ^= value + 0x9e3779b9 + (h << 6) + (h >> 2);
        };
        const auto mix_coord = [&mix](const openvdb::Coord& coord) {
            mix(std::hash<int>{}(static_cast<int>(coord.x())));
            mix(std::hash<int>{}(static_cast<int>(coord.y())));
            mix(std::hash<int>{}(static_cast<int>(coord.z())));
        };

        mix(std::hash<int>{}(static_cast<int>(grid->getGridClass())));
        mix(std::hash<float>{}(grid->background()));
        mix(std::hash<size_t>{}(static_cast<size_t>(grid->activeVoxelCount())));

        const openvdb::Vec3d voxel_size = grid->voxelSize();
        for (int i = 0; i < 3; ++i) {
            mix(std::hash<double>{}(voxel_size[i]));
        }

        const openvdb::CoordBBox bbox = grid->evalActiveVoxelBoundingBox();
        mix_coord(bbox.min());
        mix_coord(bbox.max());

        constexpr size_t kMaxSampledValues = 1024;
        size_t sampled = 0;
        for (auto iter = grid->cbeginValueOn();
             iter && sampled < kMaxSampledValues;
             ++iter, ++sampled) {
            mix_coord(iter.getCoord());
            mix(std::hash<float>{}(*iter));
        }

        return h;
    }

    void add_grid(openvdb::FloatGrid::Ptr grid);

    void apply_transform(const glm::mat4& transform) override;
    void write_disk(const std::string& file_name) const;
    openvdb::FloatGrid::ConstPtr get_grid() const;

   private:
    openvdb::FloatGrid::Ptr grid;
};
RUZINO_NAMESPACE_CLOSE_SCOPE
