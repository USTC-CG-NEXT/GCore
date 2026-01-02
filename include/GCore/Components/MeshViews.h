#pragma once
#include <GCore/api.h>

#include <Eigen/Eigen>
#include <string>

#include "glm/vec2.hpp"
#include "glm/vec3.hpp"

#ifdef GEOM_USD_EXTENSION
#include <pxr/base/gf/matrix3d.h>
#include <pxr/base/gf/matrix3f.h>
#include <pxr/base/gf/matrix4d.h>
#include <pxr/base/gf/matrix4f.h>
#include <pxr/base/gf/vec2f.h>
#include <pxr/base/gf/vec3f.h>
#include <pxr/base/gf/vec4f.h>
#include <pxr/base/vt/array.h>
#endif

RUZINO_NAMESPACE_OPEN_SCOPE
struct MeshComponent;

struct GEOMETRY_API ConstMeshIGLView {
    explicit ConstMeshIGLView(const MeshComponent& mesh);

    // Vertex positions as Nx3 matrix (using float for better compatibility)
    Eigen::MatrixXf get_vertices() const;

    // Face indices as Fx3 matrix (assuming triangular faces)
    Eigen::MatrixXi get_faces() const;

    // Face vertex counts
    Eigen::VectorXi get_face_vertex_counts() const;

    // Face vertex indices (for general polygons)
    Eigen::VectorXi get_face_vertex_indices() const;

    // Normals as Nx3 matrix
    Eigen::MatrixXf get_normals() const;

    // UV coordinates as Nx2 matrix
    Eigen::MatrixXf get_uv_coordinates() const;

    // Display colors as Nx3 matrix
    Eigen::MatrixXf get_display_colors() const;

    // Get scalar quantities
    Eigen::VectorXf get_vertex_scalar_quantity(const std::string& name) const;
    Eigen::VectorXf get_face_scalar_quantity(const std::string& name) const;

    // Get vector quantities
    Eigen::MatrixXf get_vertex_vector_quantity(const std::string& name) const;
    Eigen::MatrixXf get_face_vector_quantity(const std::string& name) const;

    // Get color quantities
    Eigen::MatrixXf get_vertex_color_quantity(const std::string& name) const;
    Eigen::MatrixXf get_face_color_quantity(const std::string& name) const;

    // Get parameterization quantities
    Eigen::MatrixXf get_vertex_parameterization_quantity(
        const std::string& name) const;
    Eigen::MatrixXf get_face_corner_parameterization_quantity(
        const std::string& name) const;

   protected:
    const MeshComponent& mesh_;

    // Utility functions for conversion using Eigen Map for zero-copy when
    // possible
    static Eigen::MatrixXf vec3f_array_to_matrix(
        const std::vector<glm::vec3>& array);
    static Eigen::MatrixXf vec2f_array_to_matrix(
        const std::vector<glm::vec2>& array);
    static Eigen::VectorXf float_array_to_vector(
        const std::vector<float>& array);
    static Eigen::VectorXi int_array_to_vector(const std::vector<int>& array);
};

struct GEOMETRY_API MeshIGLView : public ConstMeshIGLView {
    MeshIGLView(MeshComponent& mesh);

    // Set vertex positions from Nx3 matrix
    void set_vertices(const Eigen::MatrixXf& vertices);

    // Set face indices from Fx3 matrix (for triangular meshes)
    void set_faces(const Eigen::MatrixXi& faces);

    // Set face vertex counts and indices (for general polygons)
    void set_face_topology(
        const Eigen::VectorXi& face_vertex_counts,
        const Eigen::VectorXi& face_vertex_indices);

    // Set normals from Nx3 matrix
    void set_normals(const Eigen::MatrixXf& normals);

    // Set UV coordinates from Nx2 matrix
    void set_uv_coordinates(const Eigen::MatrixXf& uv_coords);

    // Set display colors from Nx3 matrix
    void set_display_colors(const Eigen::MatrixXf& colors);

    // Set scalar quantities
    void set_vertex_scalar_quantity(
        const std::string& name,
        const Eigen::VectorXf& values);
    void set_face_scalar_quantity(
        const std::string& name,
        const Eigen::VectorXf& values);

    // Set vector quantities
    void set_vertex_vector_quantity(
        const std::string& name,
        const Eigen::MatrixXf& vectors);
    void set_face_vector_quantity(
        const std::string& name,
        const Eigen::MatrixXf& vectors);

    // Set color quantities
    void set_vertex_color_quantity(
        const std::string& name,
        const Eigen::MatrixXf& colors);
    void set_face_color_quantity(
        const std::string& name,
        const Eigen::MatrixXf& colors);

    // Set parameterization quantities
    void set_vertex_parameterization_quantity(
        const std::string& name,
        const Eigen::MatrixXf& params);
    void set_face_corner_parameterization_quantity(
        const std::string& name,
        const Eigen::MatrixXf& params);

   private:
    MeshComponent& mutable_mesh_;

    // Utility functions for conversion
    static std::vector<glm::vec3> matrix_to_vec3f_array(
        const Eigen::MatrixXf& matrix);
    static std::vector<glm::vec2> matrix_to_vec2f_array(
        const Eigen::MatrixXf& matrix);
    static std::vector<float> vector_to_float_array(
        const Eigen::VectorXf& vector);
    static std::vector<int> vector_to_int_array(const Eigen::VectorXi& vector);
};
#ifdef GEOM_USD_EXTENSION

// USD View for working with USD types
struct GEOMETRY_API ConstMeshUSDView {
    explicit ConstMeshUSDView(const MeshComponent& mesh);

    // Get vertex positions as VtArray
    pxr::VtArray<pxr::GfVec3f> get_vertices() const;

    // Get face vertex counts and indices
    pxr::VtArray<int> get_face_vertex_counts() const;
    pxr::VtArray<int> get_face_vertex_indices() const;

    // Get normals as VtArray
    pxr::VtArray<pxr::GfVec3f> get_normals() const;

    // Get UV coordinates as VtArray
    pxr::VtArray<pxr::GfVec2f> get_uv_coordinates() const;

    // Get display colors as VtArray
    pxr::VtArray<pxr::GfVec3f> get_display_colors() const;

    // Get scalar quantities
    pxr::VtArray<float> get_vertex_scalar_quantity(
        const std::string& name) const;
    pxr::VtArray<float> get_face_scalar_quantity(const std::string& name) const;

    // Get vector quantities
    pxr::VtArray<pxr::GfVec3f> get_vertex_vector_quantity(
        const std::string& name) const;
    pxr::VtArray<pxr::GfVec3f> get_face_vector_quantity(
        const std::string& name) const;

    // Get color quantities
    pxr::VtArray<pxr::GfVec3f> get_vertex_color_quantity(
        const std::string& name) const;
    pxr::VtArray<pxr::GfVec3f> get_face_color_quantity(
        const std::string& name) const;

    // Get parameterization quantities
    pxr::VtArray<pxr::GfVec2f> get_vertex_parameterization_quantity(
        const std::string& name) const;
    pxr::VtArray<pxr::GfVec2f> get_face_corner_parameterization_quantity(
        const std::string& name) const;

   protected:
    const MeshComponent& mesh_;

    // Utility functions for conversion
    static pxr::VtArray<pxr::GfVec3f> vec3f_array_to_vt_array(
        const std::vector<glm::vec3>& array);
    static pxr::VtArray<pxr::GfVec2f> vec2f_array_to_vt_array(
        const std::vector<glm::vec2>& array);
    static pxr::VtArray<float> float_array_to_vt_array(
        const std::vector<float>& array);
    static pxr::VtArray<int> int_array_to_vt_array(
        const std::vector<int>& array);
};

struct GEOMETRY_API MeshUSDView : public ConstMeshUSDView {
    MeshUSDView(MeshComponent& mesh);

    // Set vertex positions from VtArray
    void set_vertices(const pxr::VtArray<pxr::GfVec3f>& vertices);

    // Set face topology from VtArrays
    void set_face_topology(
        const pxr::VtArray<int>& face_vertex_counts,
        const pxr::VtArray<int>& face_vertex_indices);

    // Set normals from VtArray
    void set_normals(const pxr::VtArray<pxr::GfVec3f>& normals);

    // Set UV coordinates from VtArray
    void set_uv_coordinates(const pxr::VtArray<pxr::GfVec2f>& uv_coords);

    // Set display colors from VtArray
    void set_display_colors(const pxr::VtArray<pxr::GfVec3f>& colors);

    // Set scalar quantities
    void set_vertex_scalar_quantity(
        const std::string& name,
        const pxr::VtArray<float>& values);
    void set_face_scalar_quantity(
        const std::string& name,
        const pxr::VtArray<float>& values);

    // Set vector quantities
    void set_vertex_vector_quantity(
        const std::string& name,
        const pxr::VtArray<pxr::GfVec3f>& vectors);
    void set_face_vector_quantity(
        const std::string& name,
        const pxr::VtArray<pxr::GfVec3f>& vectors);

    // Set color quantities
    void set_vertex_color_quantity(
        const std::string& name,
        const pxr::VtArray<pxr::GfVec3f>& colors);
    void set_face_color_quantity(
        const std::string& name,
        const pxr::VtArray<pxr::GfVec3f>& colors);

    // Set parameterization quantities
    void set_vertex_parameterization_quantity(
        const std::string& name,
        const pxr::VtArray<pxr::GfVec2f>& params);
    void set_face_corner_parameterization_quantity(
        const std::string& name,
        const pxr::VtArray<pxr::GfVec2f>& params);

   private:
    MeshComponent& mutable_mesh_;

    // Utility functions for conversion
    static std::vector<glm::vec3> vt_array_to_vec3f_array(
        const pxr::VtArray<pxr::GfVec3f>& array);
    static std::vector<glm::vec2> vt_array_to_vec2f_array(
        const pxr::VtArray<pxr::GfVec2f>& array);
    static std::vector<float> vt_array_to_float_array(
        const pxr::VtArray<float>& array);
    static std::vector<int> vt_array_to_int_array(
        const pxr::VtArray<int>& array);
};
#endif

RUZINO_NAMESPACE_CLOSE_SCOPE
