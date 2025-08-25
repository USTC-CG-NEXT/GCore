#include <GCore/Components/MeshViews.h>

#include "GCore/Components/MeshComponent.h"

USTC_CG_NAMESPACE_OPEN_SCOPE

// ConstMeshIGLView Implementation
ConstMeshIGLView::ConstMeshIGLView(const MeshComponent& mesh) : mesh_(mesh)
{
}

Eigen::MatrixXf ConstMeshIGLView::get_vertices() const
{
    return vec3f_array_to_matrix(mesh_.get_vertices());
}

Eigen::MatrixXi ConstMeshIGLView::get_faces() const
{
    auto face_vertex_counts = mesh_.get_face_vertex_counts();
    auto face_vertex_indices = mesh_.get_face_vertex_indices();

    // Count triangular faces only
    int triangle_count = 0;
    for (int count : face_vertex_counts) {
        if (count == 3)
            triangle_count++;
    }

    Eigen::MatrixXi faces(triangle_count, 3);
    int face_idx = 0;
    int vertex_idx = 0;

    for (int count : face_vertex_counts) {
        if (count == 3) {
            faces(face_idx, 0) = face_vertex_indices[vertex_idx];
            faces(face_idx, 1) = face_vertex_indices[vertex_idx + 1];
            faces(face_idx, 2) = face_vertex_indices[vertex_idx + 2];
            face_idx++;
        }
        vertex_idx += count;
    }

    return faces;
}

Eigen::VectorXi ConstMeshIGLView::get_face_vertex_counts() const
{
    return int_array_to_vector(mesh_.get_face_vertex_counts());
}

Eigen::VectorXi ConstMeshIGLView::get_face_vertex_indices() const
{
    return int_array_to_vector(mesh_.get_face_vertex_indices());
}

Eigen::MatrixXf ConstMeshIGLView::get_normals() const
{
    return vec3f_array_to_matrix(mesh_.get_normals());
}

Eigen::MatrixXf ConstMeshIGLView::get_uv_coordinates() const
{
    return vec2f_array_to_matrix(mesh_.get_texcoords_array());
}

Eigen::MatrixXf ConstMeshIGLView::get_display_colors() const
{
    return vec3f_array_to_matrix(mesh_.get_display_color());
}

Eigen::VectorXf ConstMeshIGLView::get_vertex_scalar_quantity(
    const std::string& name) const
{
    return float_array_to_vector(mesh_.get_vertex_scalar_quantity(name));
}

Eigen::VectorXf ConstMeshIGLView::get_face_scalar_quantity(
    const std::string& name) const
{
    return float_array_to_vector(mesh_.get_face_scalar_quantity(name));
}

Eigen::MatrixXf ConstMeshIGLView::get_vertex_vector_quantity(
    const std::string& name) const
{
    return vec3f_array_to_matrix(mesh_.get_vertex_vector_quantity(name));
}

Eigen::MatrixXf ConstMeshIGLView::get_face_vector_quantity(
    const std::string& name) const
{
    return vec3f_array_to_matrix(mesh_.get_face_vector_quantity(name));
}

Eigen::MatrixXf ConstMeshIGLView::get_vertex_color_quantity(
    const std::string& name) const
{
    return vec3f_array_to_matrix(mesh_.get_vertex_color_quantity(name));
}

Eigen::MatrixXf ConstMeshIGLView::get_face_color_quantity(
    const std::string& name) const
{
    return vec3f_array_to_matrix(mesh_.get_face_color_quantity(name));
}

Eigen::MatrixXf ConstMeshIGLView::get_vertex_parameterization_quantity(
    const std::string& name) const
{
    return vec2f_array_to_matrix(
        mesh_.get_vertex_parameterization_quantity(name));
}

Eigen::MatrixXf ConstMeshIGLView::get_face_corner_parameterization_quantity(
    const std::string& name) const
{
    return vec2f_array_to_matrix(
        mesh_.get_face_corner_parameterization_quantity(name));
}

// Utility functions for ConstMeshIGLView
Eigen::MatrixXf ConstMeshIGLView::vec3f_array_to_matrix(
    const std::vector<glm::vec3>& array)
{
    if (array.empty()) {
        return Eigen::MatrixXf(0, 3);
    }

    // Use Eigen Map for zero-copy when possible
    // Note: glm::vec3 uses floats, so we can map directly
    Eigen::MatrixXf matrix(array.size(), 3);
    for (size_t i = 0; i < array.size(); ++i) {
        matrix.row(i) = Eigen::Map<const Eigen::Vector3f>(&array[i][0]);
    }
    return matrix;
}

Eigen::MatrixXf ConstMeshIGLView::vec2f_array_to_matrix(
    const std::vector<glm::vec2>& array)
{
    if (array.empty()) {
        return Eigen::MatrixXf(0, 2);
    }

    Eigen::MatrixXf matrix(array.size(), 2);
    for (size_t i = 0; i < array.size(); ++i) {
        matrix.row(i) = Eigen::Map<const Eigen::Vector2f>(&array[i][0]);
    }
    return matrix;
}

Eigen::VectorXf ConstMeshIGLView::float_array_to_vector(
    const std::vector<float>& array)
{
    if (array.empty()) {
        return Eigen::VectorXf(0);
    }

    // Use Eigen Map for zero-copy
    return Eigen::Map<const Eigen::VectorXf>(array.data(), array.size());
}

Eigen::VectorXi ConstMeshIGLView::int_array_to_vector(
    const std::vector<int>& array)
{
    Eigen::VectorXi vector(array.size());
    std::memcpy(vector.data(), array.data(), array.size() * sizeof(int));
    return vector;
}

// MeshIGLView Implementation
MeshIGLView::MeshIGLView(MeshComponent& mesh)
    : ConstMeshIGLView(mesh),
      mutable_mesh_(mesh)
{
}

void MeshIGLView::set_vertices(const Eigen::MatrixXf& vertices)
{
    mutable_mesh_.set_vertices(matrix_to_vec3f_array(vertices));
}

void MeshIGLView::set_faces(const Eigen::MatrixXi& faces)
{
    // Convert triangular faces to face vertex counts and indices
    std::vector<int> face_vertex_counts(faces.rows(), 3);
    std::vector<int> face_vertex_indices(faces.rows() * 3);

    for (int i = 0; i < faces.rows(); ++i) {
        face_vertex_indices[i * 3] = faces(i, 0);
        face_vertex_indices[i * 3 + 1] = faces(i, 1);
        face_vertex_indices[i * 3 + 2] = faces(i, 2);
    }

    mutable_mesh_.set_face_vertex_counts(face_vertex_counts);
    mutable_mesh_.set_face_vertex_indices(face_vertex_indices);
}

void MeshIGLView::set_face_topology(
    const Eigen::VectorXi& face_vertex_counts,
    const Eigen::VectorXi& face_vertex_indices)
{
    mutable_mesh_.set_face_vertex_counts(
        vector_to_int_array(face_vertex_counts));
    mutable_mesh_.set_face_vertex_indices(
        vector_to_int_array(face_vertex_indices));
}

void MeshIGLView::set_normals(const Eigen::MatrixXf& normals)
{
    mutable_mesh_.set_normals(matrix_to_vec3f_array(normals));
}

void MeshIGLView::set_uv_coordinates(const Eigen::MatrixXf& uv_coords)
{
    mutable_mesh_.set_texcoords_array(matrix_to_vec2f_array(uv_coords));
}

void MeshIGLView::set_display_colors(const Eigen::MatrixXf& colors)
{
    mutable_mesh_.set_display_color(matrix_to_vec3f_array(colors));
}

void MeshIGLView::set_vertex_scalar_quantity(
    const std::string& name,
    const Eigen::VectorXf& values)
{
    mutable_mesh_.add_vertex_scalar_quantity(
        name, vector_to_float_array(values));
}

void MeshIGLView::set_face_scalar_quantity(
    const std::string& name,
    const Eigen::VectorXf& values)
{
    mutable_mesh_.add_face_scalar_quantity(name, vector_to_float_array(values));
}

void MeshIGLView::set_vertex_vector_quantity(
    const std::string& name,
    const Eigen::MatrixXf& vectors)
{
    mutable_mesh_.add_vertex_vector_quantity(
        name, matrix_to_vec3f_array(vectors));
}

void MeshIGLView::set_face_vector_quantity(
    const std::string& name,
    const Eigen::MatrixXf& vectors)
{
    mutable_mesh_.add_face_vector_quantity(
        name, matrix_to_vec3f_array(vectors));
}

void MeshIGLView::set_vertex_color_quantity(
    const std::string& name,
    const Eigen::MatrixXf& colors)
{
    mutable_mesh_.add_vertex_color_quantity(
        name, matrix_to_vec3f_array(colors));
}

void MeshIGLView::set_face_color_quantity(
    const std::string& name,
    const Eigen::MatrixXf& colors)
{
    mutable_mesh_.add_face_color_quantity(name, matrix_to_vec3f_array(colors));
}

void MeshIGLView::set_vertex_parameterization_quantity(
    const std::string& name,
    const Eigen::MatrixXf& params)
{
    mutable_mesh_.add_vertex_parameterization_quantity(
        name, matrix_to_vec2f_array(params));
}

void MeshIGLView::set_face_corner_parameterization_quantity(
    const std::string& name,
    const Eigen::MatrixXf& params)
{
    mutable_mesh_.add_face_corner_parameterization_quantity(
        name, matrix_to_vec2f_array(params));
}

// Utility functions for MeshIGLView
std::vector<glm::vec3> MeshIGLView::matrix_to_vec3f_array(
    const Eigen::MatrixXf& matrix)
{
    std::vector<glm::vec3> array(matrix.rows());
    for (int i = 0; i < matrix.rows(); ++i) {
        array[i] = glm::vec3(matrix(i, 0), matrix(i, 1), matrix(i, 2));
    }
    return array;
}

std::vector<glm::vec2> MeshIGLView::matrix_to_vec2f_array(
    const Eigen::MatrixXf& matrix)
{
    std::vector<glm::vec2> array(matrix.rows());
    for (int i = 0; i < matrix.rows(); ++i) {
        array[i] = glm::vec2(matrix(i, 0), matrix(i, 1));
    }
    return array;
}

std::vector<float> MeshIGLView::vector_to_float_array(
    const Eigen::VectorXf& vector)
{
    std::vector<float> array(vector.size());
    // Use zero-copy with Map when possible
    std::memcpy(array.data(), vector.data(), vector.size() * sizeof(float));
    return array;
}

std::vector<int> MeshIGLView::vector_to_int_array(const Eigen::VectorXi& vector)
{
    std::vector<int> array(vector.size());
    std::memcpy(array.data(), vector.data(), vector.size() * sizeof(int));
    return array;
}

MeshIGLView MeshComponent::get_igl_view()
{
    return MeshIGLView(*this);
}

ConstMeshIGLView MeshComponent::get_igl_view() const
{
    return ConstMeshIGLView(*this);
}
#ifdef GEOM_USD_EXTENSION

// ConstMeshUSDView Implementation
ConstMeshUSDView::ConstMeshUSDView(const MeshComponent& mesh) : mesh_(mesh)
{
}

pxr::VtArray<pxr::GfVec3f> ConstMeshUSDView::get_vertices() const
{
    return vec3f_array_to_vt_array(mesh_.get_vertices());
}

pxr::VtArray<int> ConstMeshUSDView::get_face_vertex_counts() const
{
    return int_array_to_vt_array(mesh_.get_face_vertex_counts());
}

pxr::VtArray<int> ConstMeshUSDView::get_face_vertex_indices() const
{
    return int_array_to_vt_array(mesh_.get_face_vertex_indices());
}

pxr::VtArray<pxr::GfVec3f> ConstMeshUSDView::get_normals() const
{
    return vec3f_array_to_vt_array(mesh_.get_normals());
}

pxr::VtArray<pxr::GfVec2f> ConstMeshUSDView::get_uv_coordinates() const
{
    return vec2f_array_to_vt_array(mesh_.get_texcoords_array());
}

pxr::VtArray<pxr::GfVec3f> ConstMeshUSDView::get_display_colors() const
{
    return vec3f_array_to_vt_array(mesh_.get_display_color());
}

pxr::VtArray<float> ConstMeshUSDView::get_vertex_scalar_quantity(
    const std::string& name) const
{
    return float_array_to_vt_array(mesh_.get_vertex_scalar_quantity(name));
}

pxr::VtArray<float> ConstMeshUSDView::get_face_scalar_quantity(
    const std::string& name) const
{
    return float_array_to_vt_array(mesh_.get_face_scalar_quantity(name));
}

pxr::VtArray<pxr::GfVec3f> ConstMeshUSDView::get_vertex_vector_quantity(
    const std::string& name) const
{
    return vec3f_array_to_vt_array(mesh_.get_vertex_vector_quantity(name));
}

pxr::VtArray<pxr::GfVec3f> ConstMeshUSDView::get_face_vector_quantity(
    const std::string& name) const
{
    return vec3f_array_to_vt_array(mesh_.get_face_vector_quantity(name));
}

pxr::VtArray<pxr::GfVec3f> ConstMeshUSDView::get_vertex_color_quantity(
    const std::string& name) const
{
    return vec3f_array_to_vt_array(mesh_.get_vertex_color_quantity(name));
}

pxr::VtArray<pxr::GfVec3f> ConstMeshUSDView::get_face_color_quantity(
    const std::string& name) const
{
    return vec3f_array_to_vt_array(mesh_.get_face_color_quantity(name));
}

pxr::VtArray<pxr::GfVec2f>
ConstMeshUSDView::get_vertex_parameterization_quantity(
    const std::string& name) const
{
    return vec2f_array_to_vt_array(
        mesh_.get_vertex_parameterization_quantity(name));
}

pxr::VtArray<pxr::GfVec2f>
ConstMeshUSDView::get_face_corner_parameterization_quantity(
    const std::string& name) const
{
    return vec2f_array_to_vt_array(
        mesh_.get_face_corner_parameterization_quantity(name));
}

// Utility functions for ConstMeshUSDView
pxr::VtArray<pxr::GfVec3f> ConstMeshUSDView::vec3f_array_to_vt_array(
    const std::vector<glm::vec3>& array)
{
    pxr::VtArray<pxr::GfVec3f> vt_array(array.size());
    for (size_t i = 0; i < array.size(); ++i) {
        vt_array[i] = pxr::GfVec3f(array[i].x, array[i].y, array[i].z);
    }
    return vt_array;
}

pxr::VtArray<pxr::GfVec2f> ConstMeshUSDView::vec2f_array_to_vt_array(
    const std::vector<glm::vec2>& array)
{
    pxr::VtArray<pxr::GfVec2f> vt_array(array.size());
    for (size_t i = 0; i < array.size(); ++i) {
        vt_array[i] = pxr::GfVec2f(array[i].x, array[i].y);
    }
    return vt_array;
}

pxr::VtArray<float> ConstMeshUSDView::float_array_to_vt_array(
    const std::vector<float>& array)
{
    pxr::VtArray<float> vt_array(array.size());
    std::memcpy(vt_array.data(), array.data(), array.size() * sizeof(float));
    return vt_array;
}

pxr::VtArray<int> ConstMeshUSDView::int_array_to_vt_array(
    const std::vector<int>& array)
{
    pxr::VtArray<int> vt_array(array.size());
    std::memcpy(vt_array.data(), array.data(), array.size() * sizeof(int));
    return vt_array;
}

// MeshUSDView Implementation
MeshUSDView::MeshUSDView(MeshComponent& mesh)
    : ConstMeshUSDView(mesh),
      mutable_mesh_(mesh)
{
}

void MeshUSDView::set_vertices(const pxr::VtArray<pxr::GfVec3f>& vertices)
{
    mutable_mesh_.set_vertices(vt_array_to_vec3f_array(vertices));
}

void MeshUSDView::set_face_topology(
    const pxr::VtArray<int>& face_vertex_counts,
    const pxr::VtArray<int>& face_vertex_indices)
{
    mutable_mesh_.set_face_vertex_counts(
        vt_array_to_int_array(face_vertex_counts));
    mutable_mesh_.set_face_vertex_indices(
        vt_array_to_int_array(face_vertex_indices));
}

void MeshUSDView::set_normals(const pxr::VtArray<pxr::GfVec3f>& normals)
{
    mutable_mesh_.set_normals(vt_array_to_vec3f_array(normals));
}

void MeshUSDView::set_uv_coordinates(
    const pxr::VtArray<pxr::GfVec2f>& uv_coords)
{
    mutable_mesh_.set_texcoords_array(vt_array_to_vec2f_array(uv_coords));
}

void MeshUSDView::set_display_colors(const pxr::VtArray<pxr::GfVec3f>& colors)
{
    mutable_mesh_.set_display_color(vt_array_to_vec3f_array(colors));
}

void MeshUSDView::set_vertex_scalar_quantity(
    const std::string& name,
    const pxr::VtArray<float>& values)
{
    mutable_mesh_.add_vertex_scalar_quantity(
        name, vt_array_to_float_array(values));
}

void MeshUSDView::set_face_scalar_quantity(
    const std::string& name,
    const pxr::VtArray<float>& values)
{
    mutable_mesh_.add_face_scalar_quantity(
        name, vt_array_to_float_array(values));
}

void MeshUSDView::set_vertex_vector_quantity(
    const std::string& name,
    const pxr::VtArray<pxr::GfVec3f>& vectors)
{
    mutable_mesh_.add_vertex_vector_quantity(
        name, vt_array_to_vec3f_array(vectors));
}

void MeshUSDView::set_face_vector_quantity(
    const std::string& name,
    const pxr::VtArray<pxr::GfVec3f>& vectors)
{
    mutable_mesh_.add_face_vector_quantity(
        name, vt_array_to_vec3f_array(vectors));
}

void MeshUSDView::set_vertex_color_quantity(
    const std::string& name,
    const pxr::VtArray<pxr::GfVec3f>& colors)
{
    mutable_mesh_.add_vertex_color_quantity(
        name, vt_array_to_vec3f_array(colors));
}

void MeshUSDView::set_face_color_quantity(
    const std::string& name,
    const pxr::VtArray<pxr::GfVec3f>& colors)
{
    mutable_mesh_.add_face_color_quantity(
        name, vt_array_to_vec3f_array(colors));
}

void MeshUSDView::set_vertex_parameterization_quantity(
    const std::string& name,
    const pxr::VtArray<pxr::GfVec2f>& params)
{
    mutable_mesh_.add_vertex_parameterization_quantity(
        name, vt_array_to_vec2f_array(params));
}

void MeshUSDView::set_face_corner_parameterization_quantity(
    const std::string& name,
    const pxr::VtArray<pxr::GfVec2f>& params)
{
    mutable_mesh_.add_face_corner_parameterization_quantity(
        name, vt_array_to_vec2f_array(params));
}

// Utility functions for MeshUSDView
std::vector<glm::vec3> MeshUSDView::vt_array_to_vec3f_array(
    const pxr::VtArray<pxr::GfVec3f>& array)
{
    std::vector<glm::vec3> result(array.size());
    for (size_t i = 0; i < array.size(); ++i) {
        result[i] = glm::vec3(array[i][0], array[i][1], array[i][2]);
    }
    return result;
}

std::vector<glm::vec2> MeshUSDView::vt_array_to_vec2f_array(
    const pxr::VtArray<pxr::GfVec2f>& array)
{
    std::vector<glm::vec2> result(array.size());
    for (size_t i = 0; i < array.size(); ++i) {
        result[i] = glm::vec2(array[i][0], array[i][1]);
    }
    return result;
}

std::vector<float> MeshUSDView::vt_array_to_float_array(
    const pxr::VtArray<float>& array)
{
    std::vector<float> result(array.size());
    std::memcpy(result.data(), array.data(), array.size() * sizeof(float));
    return result;
}

std::vector<int> MeshUSDView::vt_array_to_int_array(
    const pxr::VtArray<int>& array)
{
    std::vector<int> result(array.size());
    std::memcpy(result.data(), array.data(), array.size() * sizeof(int));
    return result;
}

MeshUSDView MeshComponent::get_usd_view()
{
    return MeshUSDView(*this);
}

ConstMeshUSDView MeshComponent::get_usd_view() const
{
    return ConstMeshUSDView(*this);
}

#endif

USTC_CG_NAMESPACE_CLOSE_SCOPE
