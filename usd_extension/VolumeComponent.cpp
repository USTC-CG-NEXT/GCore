
#include "GCore/Components/VolumeComponent.h"

RUZINO_NAMESPACE_OPEN_SCOPE
GeometryComponentHandle VolumeComponent::copy(Geometry* operand) const
{
    auto new_component = std::make_shared<VolumeComponent>(operand);
    // Geometry copies share their components and rely on copy() to detach
    // them, so handing out the same grid would let a write through one
    // Geometry show up in the other.
    new_component->add_grid(grid ? grid->deepCopy() : nullptr);
    return new_component;
}

std::string VolumeComponent::to_string() const
{
    return std::string("VolumeComponent");
}

void VolumeComponent::add_grid(openvdb::FloatGrid::Ptr grid)
{
    this->grid = grid;
}

void VolumeComponent::apply_transform(const glm::mat4& transform)
{
    if (!grid) {
        return;
    }

    // OpenVDB maps index space to world space with row vectors
    // (world = index * W), whereas glm is column-vector based
    // (world = T * index), so the transform has to be transposed on the way
    // in and post-multiplied onto the existing grid transform.
    // glm is stored column-major (transform[column][row]) while Mat4d is
    // row-major, so reading both with the same index pair already transposes.
    openvdb::math::Mat4d vdb_transform;
    for (int row = 0; row < 4; ++row) {
        for (int column = 0; column < 4; ++column) {
            vdb_transform(row, column) =
                static_cast<double>(transform[row][column]);
        }
    }

    auto grid_transform = grid->transform().copy();
    grid_transform->postMult(vdb_transform);
    grid->setTransform(grid_transform);
}

void VolumeComponent::write_disk(const std::string& file_name) const
{
    openvdb::io::File file(file_name);
    openvdb::GridPtrVec grids;
    grids.push_back(grid);
    file.write(grids);
    file.close();
}

openvdb::FloatGrid::ConstPtr VolumeComponent::get_grid() const
{
    return grid;
}

RUZINO_NAMESPACE_CLOSE_SCOPE
