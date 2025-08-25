#pragma once

#include <string>

#include "GCore/Components.h"
#include "GCore/GOP.h"

USTC_CG_NAMESPACE_OPEN_SCOPE
struct GEOMETRY_API PointsComponent : public GeometryComponent {
    explicit PointsComponent(Geometry* attached_operand);

    std::string to_string() const override;

    void apply_transform(const glm::mat4& transform) override;

    GeometryComponentHandle copy(Geometry* operand) const override;

    [[nodiscard]] std::vector<glm::vec3> get_vertices() const
    {
        return vertices;
    }

    [[nodiscard]] std::vector<glm::vec3> get_display_color() const
    {
        return displayColor;
    }

    [[nodiscard]] std::vector<float> get_width() const
    {
        return width;
    }

    void set_vertices(const std::vector<glm::vec3>& vertices)
    {
        this->vertices = vertices;
    }

    void set_display_color(const std::vector<glm::vec3>& display_color)
    {
        this->displayColor = display_color;
    }

    void set_width(const std::vector<float>& width)
    {
        this->width = width;
    }

    void set_normals(const std::vector<glm::vec3>& normals)
    {
        this->normals = normals;
    }

    [[nodiscard]] std::vector<glm::vec3> get_normals() const
    {
        return normals;
    }

   private:
    std::vector<glm::vec3> vertices;
    std::vector<glm::vec3> normals;
    std::vector<glm::vec3> displayColor;
    std::vector<float> width;
};

USTC_CG_NAMESPACE_CLOSE_SCOPE
