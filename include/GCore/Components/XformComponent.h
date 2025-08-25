#pragma once

#include "GCore/Components.h"
#include "GCore/api.h"


USTC_CG_NAMESPACE_OPEN_SCOPE
// Stores the chain of transformation

class GEOMETRY_API XformComponent : public GeometryComponent {
   public:
    GeometryComponentHandle copy(Geometry* operand) const override;
    std::string to_string() const override;

    explicit XformComponent(Geometry* attached_operand)
        : GeometryComponent(attached_operand)
    {
    }

    void apply_transform(const glm::mat4& transform) override;

    glm::mat4 get_transform() const;

    std::vector<glm::vec3> translation;
    std::vector<glm::vec3> scale;
    std::vector<glm::vec3> rotation;
};

USTC_CG_NAMESPACE_CLOSE_SCOPE
