#pragma once
#include <pxr/usd/usdSkel/skeleton.h>

#include <string>

#include "GCore/Components.h"
#include "GCore/GOP.h"

USTC_CG_NAMESPACE_OPEN_SCOPE
struct GEOMETRY_API SkelComponent : public GeometryComponent {
    explicit SkelComponent(Geometry* attached_operand)
        : GeometryComponent(attached_operand)
    {
    }

    void apply_transform(const glm::mat4& transform) override
    {
        throw std::runtime_error("Not implemented");
    }

    std::string to_string() const override;

    pxr::VtTokenArray jointOrder;
    pxr::UsdSkelTopology topology;
    pxr::VtArray<pxr::GfMatrix4f> localTransforms;
    pxr::VtArray<pxr::GfMatrix4d> bindTransforms;
    pxr::VtArray<float> jointWeight;
    pxr::VtArray<int> jointIndices;

    GeometryComponentHandle copy(Geometry* operand) const override;
};

USTC_CG_NAMESPACE_CLOSE_SCOPE
