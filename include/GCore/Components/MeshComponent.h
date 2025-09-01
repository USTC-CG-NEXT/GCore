#pragma once

#include <Eigen/Eigen>
#include <string>

#include "GCore/Components.h"
#include "GCore/GOP.h"
#include "MeshViews.h"

USTC_CG_NAMESPACE_OPEN_SCOPE

struct GEOMETRY_API MeshComponent : public GeometryComponent {
    explicit MeshComponent(Geometry* attached_operand);

    ~MeshComponent() override;

    void apply_transform(const glm::mat4& transform) override;

    std::string to_string() const override;
    GeometryComponentHandle copy(Geometry* operand) const override;

    MeshIGLView get_igl_view();
    ConstMeshIGLView get_igl_view() const;
#ifdef GEOM_USD_EXTENSION
    MeshUSDView get_usd_view();
    ConstMeshUSDView get_usd_view() const;
#endif

    [[nodiscard]] std::vector<glm::vec3> get_vertices() const
    {
        return vertices;
    }

    [[nodiscard]] std::vector<int> get_face_vertex_counts() const
    {
        return faceVertexCounts;
    }

    [[nodiscard]] std::vector<int> get_face_vertex_indices() const
    {
        return faceVertexIndices;
    }

    [[nodiscard]] std::vector<glm::vec3> get_normals() const
    {
        return normals;
    }

    [[nodiscard]] std::vector<glm::vec3> get_display_color() const
    {
        return displayColor;
    }

    [[nodiscard]] std::vector<glm::vec2> get_texcoords_array() const
    {
        return texcoordsArray;
    }

    void set_vertices(const std::vector<glm::vec3>& vertices)
    {
        this->vertices = vertices;
    }

    void set_face_vertex_counts(const std::vector<int>& face_vertex_counts)
    {
        this->faceVertexCounts = face_vertex_counts;
    }

    void set_face_vertex_indices(const std::vector<int>& face_vertex_indices)
    {
        this->faceVertexIndices = face_vertex_indices;
    }

    void set_normals(const std::vector<glm::vec3>& normals)
    {
#if USE_USD_SCRATCH_BUFFER
        mesh.CreateNormalsAttr().Set(normals);
#else
        this->normals = normals;
#endif
    }

    void set_texcoords_array(const std::vector<glm::vec2>& texcoords_array)
    {
#if USE_USD_SCRATCH_BUFFER
        auto PrimVarAPI = pxr::UsdGeomPrimvarsAPI(mesh);
        auto primvar = PrimVarAPI.CreatePrimvar(
            pxr::TfToken("UVMap"), pxr::SdfValueTypeNames->TexCoord2fArray);
        primvar.Set(texcoords_array);

        if (get_texcoords_array().size() == get_vertices().size()) {
            primvar.SetInterpolation(pxr::UsdGeomTokens->vertex);
        }
        else {
            primvar.SetInterpolation(pxr::UsdGeomTokens->faceVarying);
        }
#else
        this->texcoordsArray = texcoords_array;
#endif
    }

    void set_display_color(const std::vector<glm::vec3>& display_color)
    {
#if USE_USD_SCRATCH_BUFFER
        auto PrimVarAPI = pxr::UsdGeomPrimvarsAPI(mesh);
        pxr::UsdGeomPrimvar colorPrimvar = PrimVarAPI.CreatePrimvar(
            pxr::TfToken("displayColor"), pxr::SdfValueTypeNames->Color3fArray);
        colorPrimvar.SetInterpolation(pxr::UsdGeomTokens->vertex);
        colorPrimvar.Set(display_color);
#else
        this->displayColor = display_color;
#endif
    }

    void append_mesh(const std::shared_ptr<MeshComponent>& mesh);

   private:
    // Local cache for mesh attributes when USD cache is not enabled
    std::vector<glm::vec3> vertices;
    std::vector<int> faceVertexCounts;
    std::vector<int> faceVertexIndices;
    std::vector<glm::vec3> normals;
    std::vector<glm::vec3> displayColor;
    std::vector<glm::vec2> texcoordsArray;
};

USTC_CG_NAMESPACE_CLOSE_SCOPE
