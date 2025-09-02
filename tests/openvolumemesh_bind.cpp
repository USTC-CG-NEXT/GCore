#include <gtest/gtest.h>

#include <glm/glm.hpp>

#include "GCore/Components/MeshComponent.h"
#include "GCore/util_openmesh_bind.h"


using namespace USTC_CG;

class OpenVolumeMeshBindTest : public ::testing::Test {
   protected:
    void SetUp() override
    {
        // Create a simple tetrahedron for testing
        geometry = std::make_shared<Geometry>();
        mesh = std::make_shared<MeshComponent>(geometry.get());
        geometry->attach_component(mesh);
    }

    void CreateSingleTetrahedron()
    {
        // Define vertices of a regular tetrahedron
        std::vector<glm::vec3> vertices = {
            { 0.0f, 0.0f, 0.0f },     // v0
            { 1.0f, 0.0f, 0.0f },     // v1
            { 0.5f, 0.866f, 0.0f },   // v2
            { 0.5f, 0.289f, 0.816f }  // v3
        };

        // Define the 4 triangular faces of the tetrahedron
        // Each face is oriented outward
        std::vector<int> faceVertexIndices = {
            // Face 0: v0, v1, v2 (bottom face)
            0,
            1,
            2,
            // Face 1: v0, v3, v1 (side face)
            0,
            3,
            1,
            // Face 2: v1, v3, v2 (side face)
            1,
            3,
            2,
            // Face 3: v2, v3, v0 (side face)
            2,
            3,
            0
        };

        std::vector<int> faceVertexCounts = { 3, 3, 3, 3 };

        mesh->set_vertices(vertices);
        mesh->set_face_vertex_indices(faceVertexIndices);
        mesh->set_face_vertex_counts(faceVertexCounts);
    }

    void CreateTwoAdjacentTetrahedra()
    {
        // Define vertices for two adjacent tetrahedra sharing a triangular face
        std::vector<glm::vec3> vertices = {
            { 0.0f, 0.0f, 0.0f },      // v0
            { 1.0f, 0.0f, 0.0f },      // v1
            { 0.5f, 0.866f, 0.0f },    // v2
            { 0.5f, 0.289f, 0.816f },  // v3 (apex of first tet)
            { 0.5f, 0.289f, -0.816f }  // v4 (apex of second tet)
        };

        // Define triangular faces for both tetrahedra
        // First tetrahedron: v0, v1, v2, v3
        // Second tetrahedron: v0, v1, v2, v4 (shares face v0,v1,v2 with first)
        std::vector<int> faceVertexIndices = {
            // First tetrahedron faces
            0,
            1,
            2,  // shared face
            0,
            3,
            1,  // side face
            1,
            3,
            2,  // side face
            2,
            3,
            0,  // side face
            // Second tetrahedron faces (excluding shared face)
            0,
            4,
            1,  // side face
            1,
            4,
            2,  // side face
            2,
            4,
            0  // side face
        };

        std::vector<int> faceVertexCounts = { 3, 3, 3, 3, 3, 3, 3 };

        mesh->set_vertices(vertices);
        mesh->set_face_vertex_indices(faceVertexIndices);
        mesh->set_face_vertex_counts(faceVertexCounts);
    }

    std::shared_ptr<Geometry> geometry;
    std::shared_ptr<MeshComponent> mesh;
};

TEST_F(OpenVolumeMeshBindTest, SingleTetrahedronReconstruction)
{
    CreateSingleTetrahedron();

    // Convert to OpenVolumeMesh
    auto volume_mesh = operand_to_openvolulemesh_from_faces(geometry.get());

    // Verify the conversion
    EXPECT_EQ(volume_mesh->n_vertices(), 4);  // Should have 4 vertices
    EXPECT_EQ(volume_mesh->n_cells(), 1);     // Should have 1 tetrahedron
    EXPECT_EQ(
        volume_mesh->n_faces(), 4);  // Should have 4 faces (auto-generated)

    // Convert back to operand
    auto reconstructed = openvolulemesh_to_operand(volume_mesh.get());
    auto reconstructed_mesh = reconstructed->get_component<MeshComponent>();

    // Verify round-trip conversion
    EXPECT_EQ(reconstructed_mesh->get_vertices().size(), 4);
    // The face representation might differ after round-trip due to internal
    // structure
}

TEST_F(OpenVolumeMeshBindTest, EmptyMeshHandling)
{
    // Test with empty mesh
    auto volume_mesh = operand_to_openvolulemesh_from_faces(geometry.get());

    EXPECT_EQ(volume_mesh->n_vertices(), 0);
    EXPECT_EQ(volume_mesh->n_cells(), 0);
}

TEST_F(OpenVolumeMeshBindTest, InvalidTrianglesHandling)
{
    // Create mesh with non-triangular faces
    std::vector<glm::vec3> vertices = {
        { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 }
    };
    std::vector<int> faceVertexIndices = {
        0, 1, 2, 3
    };  // One quad instead of triangles
    std::vector<int> faceVertexCounts = { 4 };

    mesh->set_vertices(vertices);
    mesh->set_face_vertex_indices(faceVertexIndices);
    mesh->set_face_vertex_counts(faceVertexCounts);

    auto volume_mesh = operand_to_openvolulemesh_from_faces(geometry.get());

    // Should handle gracefully - no tetrahedra should be created from quad
    EXPECT_EQ(volume_mesh->n_cells(), 0);
}

TEST_F(OpenVolumeMeshBindTest, TwoAdjacentTetrahedraReconstruction)
{
    CreateTwoAdjacentTetrahedra();

    // Convert to OpenVolumeMesh
    auto volume_mesh = operand_to_openvolulemesh_from_faces(geometry.get());

    // Verify the conversion
    EXPECT_EQ(volume_mesh->n_vertices(), 5);  // Should have 5 vertices
    EXPECT_EQ(volume_mesh->n_cells(), 2);     // Should have 2 tetrahedra
    // Note: n_faces() will be auto-computed by OpenVolumeMesh and may differ

    // Convert back to operand
    auto reconstructed = openvolulemesh_to_operand(volume_mesh.get());
    auto reconstructed_mesh = reconstructed->get_component<MeshComponent>();

    // Verify round-trip conversion
    EXPECT_EQ(reconstructed_mesh->get_vertices().size(), 5);
    // The face representation will change - OpenVolumeMesh will create
    // internal data structure with proper adjacency information
}
