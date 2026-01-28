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

#if RUZINO_WITH_CUDA
#include <RHI/internal/cuda_extension.hpp>
#endif
#include <nvrhi/nvrhi.h>

RUZINO_NAMESPACE_OPEN_SCOPE
struct MeshComponent;

struct GEOMETRY_API ConstMeshIGLView {
    explicit ConstMeshIGLView(const MeshComponent& mesh);
    ~ConstMeshIGLView();

    // Delete copy and move to prevent bypassing view lock
    ConstMeshIGLView(const ConstMeshIGLView&) = delete;
    ConstMeshIGLView& operator=(const ConstMeshIGLView&) = delete;
    ConstMeshIGLView(ConstMeshIGLView&&) = delete;
    ConstMeshIGLView& operator=(ConstMeshIGLView&&) = delete;

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
    ~MeshIGLView();

    // Delete copy and move to prevent bypassing view lock
    MeshIGLView(const MeshIGLView&) = delete;
    MeshIGLView& operator=(const MeshIGLView&) = delete;
    MeshIGLView(MeshIGLView&&) = delete;
    MeshIGLView& operator=(MeshIGLView&&) = delete;

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
    ~ConstMeshUSDView();

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
    ~MeshUSDView();

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

#if RUZINO_WITH_CUDA
// CUDA View for working with GPU buffers via CudaLinearBuffer
struct GEOMETRY_API ConstMeshCUDAView {
    explicit ConstMeshCUDAView(const MeshComponent& mesh);
    ~ConstMeshCUDAView();

    // Get vertex positions as CUDA buffer (lazy creation)
    cuda::CUDALinearBufferHandle get_vertices() const;

    // Get face vertex counts and indices
    cuda::CUDALinearBufferHandle get_face_vertex_counts() const;
    cuda::CUDALinearBufferHandle get_face_vertex_indices() const;

    // Get normals as CUDA buffer
    cuda::CUDALinearBufferHandle get_normals() const;

    // Get UV coordinates as CUDA buffer
    cuda::CUDALinearBufferHandle get_uv_coordinates() const;

    // Get display colors as CUDA buffer
    cuda::CUDALinearBufferHandle get_display_colors() const;

    // Get scalar quantities
    cuda::CUDALinearBufferHandle get_vertex_scalar_quantity(
        const std::string& name) const;
    cuda::CUDALinearBufferHandle get_face_scalar_quantity(
        const std::string& name) const;

    // Get vector quantities
    cuda::CUDALinearBufferHandle get_vertex_vector_quantity(
        const std::string& name) const;
    cuda::CUDALinearBufferHandle get_face_vector_quantity(
        const std::string& name) const;

    // Get parameterization quantities
    cuda::CUDALinearBufferHandle get_vertex_parameterization_quantity(
        const std::string& name) const;
    cuda::CUDALinearBufferHandle get_face_corner_parameterization_quantity(
        const std::string& name) const;

   protected:
    const MeshComponent& mesh_;

    // Lazy-created CUDA buffers
    mutable std::map<std::string, cuda::CUDALinearBufferHandle> cuda_buffers_;

    cuda::CUDALinearBufferHandle get_or_create_buffer(
        const std::string& key,
        const void* data,
        size_t element_count,
        size_t element_size) const;
};

struct GEOMETRY_API MeshCUDAView : public ConstMeshCUDAView {
    MeshCUDAView(MeshComponent& mesh);
    ~MeshCUDAView();  // RAII: sync data back to CPU on destruction

    // Set vertex positions from CUDA buffer
    void set_vertices(cuda::CUDALinearBufferHandle vertices);

    // Set face topology from CUDA buffers
    void set_face_topology(
        cuda::CUDALinearBufferHandle face_vertex_counts,
        cuda::CUDALinearBufferHandle face_vertex_indices);

    // Set normals from CUDA buffer
    void set_normals(cuda::CUDALinearBufferHandle normals);

    // Set UV coordinates from CUDA buffer
    void set_uv_coordinates(cuda::CUDALinearBufferHandle uv_coords);

    // Set display colors from CUDA buffer
    void set_display_colors(cuda::CUDALinearBufferHandle colors);

    // Set scalar quantities
    void set_vertex_scalar_quantity(
        const std::string& name,
        cuda::CUDALinearBufferHandle values);
    void set_face_scalar_quantity(
        const std::string& name,
        cuda::CUDALinearBufferHandle values);

    // Set vector quantities
    void set_vertex_vector_quantity(
        const std::string& name,
        cuda::CUDALinearBufferHandle vectors);
    void set_face_vector_quantity(
        const std::string& name,
        cuda::CUDALinearBufferHandle vectors);

    // Set parameterization quantities
    void set_vertex_parameterization_quantity(
        const std::string& name,
        cuda::CUDALinearBufferHandle params);
    void set_face_corner_parameterization_quantity(
        const std::string& name,
        cuda::CUDALinearBufferHandle params);

   private:
    MeshComponent& mutable_mesh_;

    // Track modified buffers for sync back to CPU
    mutable std::set<std::string> modified_buffers_;

    void sync_buffer_to_cpu(const std::string& buffer_name);
    void sync_all_to_cpu();
};
#endif

// NVRHI View for working with nvrhi::IBuffer
struct GEOMETRY_API ConstMeshNVRHIView {
    explicit ConstMeshNVRHIView(
        const MeshComponent& mesh,
        nvrhi::IDevice* device);
    ~ConstMeshNVRHIView();

    // Get vertex positions as NVRHI buffer (lazy creation)
    // commandList is optional - if provided, buffer will be initialized
    // immediately
    nvrhi::IBuffer* get_vertices(
        nvrhi::ICommandList* commandList = nullptr) const;

    // Get face vertex counts and indices
    nvrhi::IBuffer* get_face_vertex_counts(
        nvrhi::ICommandList* commandList = nullptr) const;
    nvrhi::IBuffer* get_face_vertex_indices(
        nvrhi::ICommandList* commandList = nullptr) const;

    // Get normals as NVRHI buffer
    nvrhi::IBuffer* get_normals(
        nvrhi::ICommandList* commandList = nullptr) const;

    // Get UV coordinates as NVRHI buffer
    nvrhi::IBuffer* get_uv_coordinates(
        nvrhi::ICommandList* commandList = nullptr) const;

    // Get display colors as NVRHI buffer
    nvrhi::IBuffer* get_display_colors(
        nvrhi::ICommandList* commandList = nullptr) const;

    // Get scalar quantities
    nvrhi::IBuffer* get_vertex_scalar_quantity(
        const std::string& name,
        nvrhi::ICommandList* commandList = nullptr) const;
    nvrhi::IBuffer* get_face_scalar_quantity(
        const std::string& name,
        nvrhi::ICommandList* commandList = nullptr) const;

    // Get vector quantities
    nvrhi::IBuffer* get_vertex_vector_quantity(
        const std::string& name,
        nvrhi::ICommandList* commandList = nullptr) const;
    nvrhi::IBuffer* get_face_vector_quantity(
        const std::string& name,
        nvrhi::ICommandList* commandList = nullptr) const;

    // Get parameterization quantities
    nvrhi::IBuffer* get_vertex_parameterization_quantity(
        const std::string& name,
        nvrhi::ICommandList* commandList = nullptr) const;
    nvrhi::IBuffer* get_face_corner_parameterization_quantity(
        const std::string& name,
        nvrhi::ICommandList* commandList = nullptr) const;

   protected:
    const MeshComponent& mesh_;
    nvrhi::IDevice* device_;

    // Lazy-created NVRHI buffers
    mutable std::map<std::string, nvrhi::BufferHandle> nvrhi_buffers_;

    nvrhi::BufferHandle get_or_create_buffer(
        const std::string& key,
        const void* data,
        size_t element_count,
        size_t element_size,
        nvrhi::ICommandList* commandList) const;
};

struct GEOMETRY_API MeshNVRHIView : public ConstMeshNVRHIView {
    MeshNVRHIView(MeshComponent& mesh, nvrhi::IDevice* device);
    ~MeshNVRHIView();  // RAII: clears modified buffer tracking

    // Set vertex positions from NVRHI buffer
    void set_vertices(nvrhi::IBuffer* vertices);

    // Set face topology from NVRHI buffers
    void set_face_topology(
        nvrhi::IBuffer* face_vertex_counts,
        nvrhi::IBuffer* face_vertex_indices);

    // Set normals from NVRHI buffer
    void set_normals(nvrhi::IBuffer* normals);

    // Set UV coordinates from NVRHI buffer
    void set_uv_coordinates(nvrhi::IBuffer* uv_coords);

    // Set display colors from NVRHI buffer
    void set_display_colors(nvrhi::IBuffer* colors);

    // Set scalar quantities
    void set_vertex_scalar_quantity(
        const std::string& name,
        nvrhi::IBuffer* values);
    void set_face_scalar_quantity(
        const std::string& name,
        nvrhi::IBuffer* values);

    // Set vector quantities
    void set_vertex_vector_quantity(
        const std::string& name,
        nvrhi::IBuffer* vectors);
    void set_face_vector_quantity(
        const std::string& name,
        nvrhi::IBuffer* vectors);

    // Set parameterization quantities
    void set_vertex_parameterization_quantity(
        const std::string& name,
        nvrhi::IBuffer* params);
    void set_face_corner_parameterization_quantity(
        const std::string& name,
        nvrhi::IBuffer* params);

   private:
    MeshComponent& mutable_mesh_;

    // Track modified buffers - user needs to manually sync if needed
    mutable std::set<std::string> modified_buffers_;
};

RUZINO_NAMESPACE_CLOSE_SCOPE
