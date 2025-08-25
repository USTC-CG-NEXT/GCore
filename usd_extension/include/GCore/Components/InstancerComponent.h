#pragma once
#include "GCore/Components.h"
#include "GCore/api.h"
#include "pxr/base/vt/array.h"

USTC_CG_NAMESPACE_OPEN_SCOPE

class GEOMETRY_API InstancerComponent final : public GeometryComponent {
   public:
    explicit InstancerComponent(Geometry* attached_operand);

    ~InstancerComponent() override;
    GeometryComponentHandle copy(Geometry* operand) const override;
    std::string to_string() const override;
    void apply_transform(const glm::mat4& transform) override;

    void add_instance(const glm::mat4& instance);
    const std::vector<glm::mat4>& get_instances() const;
    pxr::VtArray<int> get_proto_indices();

   private:
    std::vector<glm::mat4> instances_;
};

USTC_CG_NAMESPACE_CLOSE_SCOPE
