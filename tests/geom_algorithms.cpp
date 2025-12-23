#include <gtest/gtest.h>

#include <memory>

#include "GCore/Components/MeshComponent.h"
#include "GCore/algorithms/intersection.h"

#ifdef GPU_GEOM_ALGORITHM

#include "RHI/rhi.hpp"
#include "glm/ext/matrix_transform.hpp"

using namespace USTC_CG;

class IntersectionTests : public ::testing::Test {
   protected:
    static void SetUpTestSuite()
    {
        RHI::init();
        init_gpu_geometry_algorithms();
    }

    static void TearDownTestSuite()
    {
        deinit_gpu_geometry_algorithms();
        RHI::shutdown();
    }

    // Helper to create a simple triangle mesh
    Geometry CreateTriangleMesh()
    {
        Geometry mesh = Geometry::CreateMesh();
        auto meshComp = mesh.get_component<MeshComponent>();

        // Create a simple triangle
        std::vector<glm::vec3> vertices;
        vertices.push_back(glm::vec3(0.0f, 0.0f, 0.0f));
        vertices.push_back(glm::vec3(1.0f, 0.0f, 0.0f));
        vertices.push_back(glm::vec3(0.0f, 1.0f, 0.0f));

        std::vector<int> faceVertexCounts;
        faceVertexCounts.push_back(3);

        std::vector<int> faceVertexIndices;
        faceVertexIndices.push_back(0);
        faceVertexIndices.push_back(1);
        faceVertexIndices.push_back(2);

        std::vector<glm::vec3> normals;
        normals.push_back(glm::vec3(0.0f, 0.0f, 1.0f));
        normals.push_back(glm::vec3(0.0f, 0.0f, 1.0f));
        normals.push_back(glm::vec3(0.0f, 0.0f, 1.0f));

        meshComp->set_vertices(vertices);
        meshComp->set_face_vertex_counts(faceVertexCounts);
        meshComp->set_face_vertex_indices(faceVertexIndices);
        meshComp->set_normals(normals);

        return mesh;
    }
};

TEST_F(IntersectionTests, BasicRayIntersection)
{
    // Create a simple mesh
    Geometry mesh = CreateTriangleMesh();

    // Create rays for intersection
    std::vector<glm::ray> rays;

    // Ray that should hit the triangle (inside the triangle bounds)
    rays.push_back(
        glm::ray(glm::vec3(0.25f, 0.25f, 1.0f), glm::vec3(0.0f, 0.0f, -1.0f)));

    // Ray that should miss the triangle
    rays.push_back(
        glm::ray(glm::vec3(-0.5f, 0.5f, 1.0f), glm::vec3(0.0f, 0.0f, -1.0f)));

    // Test intersection
    std::vector<PointSample> samples = Intersect(rays, mesh);

    // Check results
    ASSERT_EQ(samples.size(), 2);

    // First ray should hit
    EXPECT_TRUE(samples[0].valid);
    EXPECT_NEAR(samples[0].position[0], 0.25f, 1e-4);
    EXPECT_NEAR(samples[0].position[1], 0.25f, 1e-4);
    EXPECT_NEAR(samples[0].position[2], 0.0f, 1e-4);

    // Normal should be correct
    EXPECT_NEAR(samples[0].normal[0], 0.0f, 1e-4);
    EXPECT_NEAR(samples[0].normal[1], 0.0f, 1e-4);
    EXPECT_NEAR(samples[0].normal[2], 1.0f, 1e-4);

    // Second ray should miss
    EXPECT_FALSE(samples[1].valid);
}

TEST_F(IntersectionTests, RayParallelToTriangle)
{
    Geometry mesh = CreateTriangleMesh();

    std::vector<glm::ray> rays;
    // Ray parallel to triangle
    rays.push_back(
        glm::ray(glm::vec3(0.5f, 0.5f, 1.0f), glm::vec3(1.0f, 0.0f, 0.0f)));

    std::vector<PointSample> samples = Intersect(rays, mesh);

    ASSERT_EQ(samples.size(), 1);
    EXPECT_FALSE(samples[0].valid);
}

TEST_F(IntersectionTests, MultipleTriangles)
{
    Geometry mesh = Geometry::CreateMesh();
    auto meshComp = mesh.get_component<MeshComponent>();

    // Create two triangles
    std::vector<glm::vec3> vertices;
    vertices.push_back(glm::vec3(0.0f, 0.0f, 0.0f));
    vertices.push_back(glm::vec3(1.0f, 0.0f, 0.0f));
    vertices.push_back(glm::vec3(0.0f, 1.0f, 0.0f));
    vertices.push_back(glm::vec3(1.0f, 1.0f, 0.0f));

    std::vector<int> faceVertexCounts;
    faceVertexCounts.push_back(3);
    faceVertexCounts.push_back(3);

    std::vector<int> faceVertexIndices;
    faceVertexIndices.push_back(0);
    faceVertexIndices.push_back(1);
    faceVertexIndices.push_back(2);
    faceVertexIndices.push_back(1);
    faceVertexIndices.push_back(3);
    faceVertexIndices.push_back(2);

    std::vector<glm::vec3> normals;
    for (int i = 0; i < 6; i++) {
        normals.push_back(glm::vec3(0.0f, 0.0f, 1.0f));
    }

    meshComp->set_vertices(vertices);
    meshComp->set_face_vertex_counts(faceVertexCounts);
    meshComp->set_face_vertex_indices(faceVertexIndices);
    meshComp->set_normals(normals);

    // Create rays for intersection
    std::vector<glm::ray> rays;
    // Hit first triangle
    rays.push_back(
        glm::ray(glm::vec3(0.25f, 0.25f, 1.0f), glm::vec3(0.0f, 0.0f, -1.0f)));
    // Hit second triangle
    rays.push_back(
        glm::ray(glm::vec3(0.75f, 0.75f, 1.0f), glm::vec3(0.0f, 0.0f, -1.0f)));

    std::vector<PointSample> samples = Intersect(rays, mesh);

    ASSERT_EQ(samples.size(), 2);
    EXPECT_TRUE(samples[0].valid);
    EXPECT_TRUE(samples[1].valid);

    EXPECT_NEAR(samples[0].position[0], 0.25f, 1e-4);
    EXPECT_NEAR(samples[0].position[1], 0.25f, 1e-4);
    EXPECT_NEAR(samples[0].position[2], 0.0f, 1e-4);

    EXPECT_NEAR(samples[1].position[0], 0.75f, 1e-4);
    EXPECT_NEAR(samples[1].position[1], 0.75f, 1e-4);
    EXPECT_NEAR(samples[1].position[2], 0.0f, 1e-4);
}

TEST_F(IntersectionTests, EmptyMesh)
{
    Geometry mesh = Geometry::CreateMesh();

    std::vector<glm::ray> rays;
    rays.push_back(
        glm::ray(glm::vec3(0.5f, 0.5f, 1.0f), glm::vec3(0.0f, 0.0f, -1.0f)));

    std::vector<PointSample> samples = Intersect(rays, mesh);

    ASSERT_EQ(samples.size(), 1);
    EXPECT_FALSE(samples[0].valid);
}

TEST_F(IntersectionTests, TransformedMesh)
{
    Geometry mesh = CreateTriangleMesh();
    auto meshComp = mesh.get_component<MeshComponent>();

    // Apply a transform to the mesh (move it to z=-1)
    glm::mat4 transform;
    transform =
        glm::translate(glm::identity<glm::mat4>(), glm::vec3(0.0, 0.0, -1.0));
    meshComp->apply_transform(transform);

    std::vector<glm::ray> rays;
    rays.push_back(
        glm::ray(glm::vec3(0.5f, 0.5f, 1.0f), glm::vec3(0.0f, 0.0f, -1.0f)));

    std::vector<PointSample> samples = Intersect(rays, mesh);

    ASSERT_EQ(samples.size(), 1);
    EXPECT_TRUE(samples[0].valid);

    EXPECT_NEAR(samples[0].position[0], 0.5f, 1e-4);
    EXPECT_NEAR(samples[0].position[1], 0.5f, 1e-4);
    EXPECT_NEAR(samples[0].position[2], -1.0f, 1e-4);
}

TEST_F(IntersectionTests, FindNeighborsClosePoints)
{
    // Create a point cloud with some close points
    Geometry point_cloud = Geometry::CreateMesh();
    auto meshComp = point_cloud.get_component<MeshComponent>();

    std::vector<glm::vec3> vertices;
    // Point 0 and 1 are close (distance = 0.01)
    vertices.push_back(glm::vec3(0.0f, 0.0f, 0.0f));
    vertices.push_back(glm::vec3(0.01f, 0.0f, 0.0f));

    // Point 2 is far from others
    vertices.push_back(glm::vec3(1.0f, 1.0f, 0.0f));

    // Point 3 and 4 are close to each other (distance = 0.015)
    vertices.push_back(glm::vec3(2.0f, 0.0f, 0.0f));
    vertices.push_back(glm::vec3(2.015f, 0.0f, 0.0f));

    meshComp->set_vertices(vertices);

    // Find neighbors
    unsigned pair_count = 0;
    float search_radius = 0.02f;  // Use 0.02 as the search radius
    std::vector<PointPairs> pairs =
        FindNeighbors(point_cloud, search_radius, pair_count);

    // Should find at least 2 pairs: (0,1) and (3,4)
    EXPECT_GT(pair_count, 0);
    EXPECT_EQ(pairs.size(), pair_count);

    std::cout << "Found " << pair_count << " neighbor pairs:" << std::endl;
}

TEST_F(IntersectionTests, FindNeighborsNoClosePoints)
{
    // Create a point cloud with points far apart
    Geometry point_cloud = Geometry::CreateMesh();
    auto meshComp = point_cloud.get_component<MeshComponent>();

    std::vector<glm::vec3> vertices;
    // All points are far from each other (> 0.02)
    vertices.push_back(glm::vec3(0.0f, 0.0f, 0.0f));
    vertices.push_back(glm::vec3(1.0f, 0.0f, 0.0f));
    vertices.push_back(glm::vec3(0.0f, 1.0f, 0.0f));
    vertices.push_back(glm::vec3(1.0f, 1.0f, 0.0f));

    meshComp->set_vertices(vertices);

    unsigned pair_count = 0;
    float search_radius = 0.02f;
    std::vector<PointPairs> pairs =
        FindNeighbors(point_cloud, search_radius, pair_count);

    // Should find no pairs
    EXPECT_EQ(pair_count, 0);
    EXPECT_EQ(pairs.size(), 0);

    std::cout << "No close neighbors found (as expected)" << std::endl;
}

TEST_F(IntersectionTests, FindNeighborsGrid)
{
    // Create a regular grid of points
    Geometry point_cloud = Geometry::CreateMesh();
    auto meshComp = point_cloud.get_component<MeshComponent>();

    std::vector<glm::vec3> vertices;
    float spacing = 0.015f;  // Close enough to be neighbors

    // Create a 3x3 grid
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            vertices.push_back(glm::vec3(i * spacing, j * spacing, 0.0f));
        }
    }

    meshComp->set_vertices(vertices);

    unsigned pair_count = 0;
    float search_radius = 0.022f;  // Use 0.022 to include diagonal neighbors
                                   // (sqrt(2)*0.015 ≈ 0.0212)
    std::vector<PointPairs> pairs =
        FindNeighbors(point_cloud, search_radius, pair_count);

    std::cout << "Grid test: Found " << pair_count
              << " neighbor pairs in 3x3 grid" << std::endl;

    // In a 3x3 grid with spacing 0.015, each internal point should have
    // neighbors
    EXPECT_GT(pair_count, 0);

    // Verify all pairs are valid
    for (unsigned i = 0; i < pair_count; ++i) {
        EXPECT_LT(pairs[i].p1, vertices.size());
        EXPECT_LT(pairs[i].p2, vertices.size());

        // Calculate actual distance
        glm::vec3 pos1 = vertices[pairs[i].p1];
        glm::vec3 pos2 = vertices[pairs[i].p2];
        float dist = glm::length(pos1 - pos2);

        std::cout << "  Pair " << i << ": (" << pairs[i].p1 << ", "
                  << pairs[i].p2 << ") distance = " << dist << std::endl;

        // Distance should be within the detection range
        EXPECT_LE(dist, search_radius);
    }
}

TEST_F(IntersectionTests, FindNeighborsEmptyPointCloud)
{
    // Create an empty point cloud
    Geometry point_cloud = Geometry::CreateMesh();

    unsigned pair_count = 0;
    float search_radius = 0.02f;
    std::vector<PointPairs> pairs =
        FindNeighbors(point_cloud, search_radius, pair_count);

    // Should return no pairs
    EXPECT_EQ(pair_count, 0);
    EXPECT_EQ(pairs.size(), 0);
}

TEST_F(IntersectionTests, FindNeighborsBoxDistance_0_2)
{
    // Test box distance with threshold 0.2
    // Box distance means: at least one axis distance < threshold
    Geometry point_cloud = Geometry::CreateMesh();
    auto meshComp = point_cloud.get_component<MeshComponent>();

    std::vector<glm::vec3> vertices;
    // Point 0 at origin
    vertices.push_back(glm::vec3(0.0f, 0.0f, 0.0f));

    // Point 1: x-axis distance = 0.15 (< 0.2), y = 0.5, z = 0.5
    // Euclidean distance = sqrt(0.15^2 + 0.5^2 + 0.5^2) ≈ 0.72 >> 0.2
    // Should be neighbor because x-axis distance < 0.2
    vertices.push_back(glm::vec3(0.15f, 0.5f, 0.5f));

    // Point 2: x = 0.5, y-axis distance = 0.18 (< 0.2), z = 0.5
    // Euclidean distance = sqrt(0.5^2 + 0.18^2 + 0.5^2) ≈ 0.73 >> 0.2
    // Should be neighbor because y-axis distance < 0.2
    vertices.push_back(glm::vec3(0.5f, 0.18f, 0.5f));

    // Point 3: x = 0.5, y = 0.5, z-axis distance = 0.19 (< 0.2)
    // Euclidean distance = sqrt(0.5^2 + 0.5^2 + 0.19^2) ≈ 0.74 >> 0.2
    // Should be neighbor because z-axis distance < 0.2
    vertices.push_back(glm::vec3(0.5f, 0.5f, 0.19f));

    // Point 4: All axis distances >= 0.2
    // x = 0.25, y = 0.25, z = 0.25, all > 0.2
    // Should NOT be neighbor
    vertices.push_back(glm::vec3(0.25f, 0.25f, 0.25f));

    meshComp->set_vertices(vertices);

    unsigned pair_count = 0;
    float search_radius = 0.2f;
    std::vector<PointPairs> pairs =
        FindNeighbors(point_cloud, search_radius, pair_count);

    std::cout << "Box distance test (threshold=0.2): Found " << pair_count
              << " neighbor pairs" << std::endl;

    // Should find 3 pairs: (1,0), (2,0), (3,0)
    EXPECT_EQ(pair_count, 3);

    // Verify each pair and check box distance property
    for (unsigned i = 0; i < pair_count; ++i) {
        glm::vec3 pos1 = vertices[pairs[i].p1];
        glm::vec3 pos2 = vertices[pairs[i].p2];

        float dx = std::abs(pos1.x - pos2.x);
        float dy = std::abs(pos1.y - pos2.y);
        float dz = std::abs(pos1.z - pos2.z);
        float euclidean_dist = glm::length(pos1 - pos2);

        std::cout << "  Pair " << i << ": (" << pairs[i].p1 << ", "
                  << pairs[i].p2 << ")" << std::endl;
        std::cout << "    Axis distances: dx=" << dx << ", dy=" << dy
                  << ", dz=" << dz << std::endl;
        std::cout << "    Euclidean distance: " << euclidean_dist << std::endl;

        // At least one axis distance should be < threshold
        bool has_close_axis = (dx < search_radius) || (dy < search_radius) ||
                              (dz < search_radius);
        EXPECT_TRUE(has_close_axis);

        // Verify it's actually a box distance, not Euclidean
        // (some pairs will have large Euclidean distance)
        std::cout << "    Box distance satisfied: "
                  << (has_close_axis ? "YES" : "NO") << std::endl;
    }
}

TEST_F(IntersectionTests, FindNeighborsBoxDistance_0_1)
{
    // Test box distance with threshold 0.1
    Geometry point_cloud = Geometry::CreateMesh();
    auto meshComp = point_cloud.get_component<MeshComponent>();

    std::vector<glm::vec3> vertices;
    // Point 0 at origin
    vertices.push_back(glm::vec3(0.0f, 0.0f, 0.0f));

    // Point 1: x-axis distance = 0.08 (< 0.1), y = 0.3, z = 0.3
    // Euclidean distance ≈ 0.43 >> 0.1
    // Should be neighbor because x-axis distance < 0.1
    vertices.push_back(glm::vec3(0.08f, 0.3f, 0.3f));

    // Point 2: x = 0.3, y-axis distance = 0.09 (< 0.1), z = 0.3
    // Euclidean distance ≈ 0.43 >> 0.1
    // Should be neighbor because y-axis distance < 0.1
    vertices.push_back(glm::vec3(0.3f, 0.09f, 0.3f));

    // Point 3: x = 0.3, y = 0.3, z-axis distance = 0.095 (< 0.1)
    // Euclidean distance ≈ 0.44 >> 0.1
    // Should be neighbor because z-axis distance < 0.1
    vertices.push_back(glm::vec3(0.3f, 0.3f, 0.095f));

    // Point 4: All axis distances >= 0.1
    // x = 0.15, y = 0.15, z = 0.15
    // Should NOT be neighbor
    vertices.push_back(glm::vec3(0.15f, 0.15f, 0.15f));

    // Point 5: Close on all axes
    // x = 0.05, y = 0.06, z = 0.07, all < 0.1
    // Should be neighbor
    vertices.push_back(glm::vec3(0.05f, 0.06f, 0.07f));

    meshComp->set_vertices(vertices);

    unsigned pair_count = 0;
    float search_radius = 0.1f;
    std::vector<PointPairs> pairs =
        FindNeighbors(point_cloud, search_radius, pair_count);

    std::cout << "Box distance test (threshold=0.1): Found " << pair_count
              << " neighbor pairs" << std::endl;

    // Should find 4 pairs: (1,0), (2,0), (3,0), (5,0)
    EXPECT_EQ(pair_count, 4);

    // Verify each pair satisfies box distance property
    for (unsigned i = 0; i < pair_count; ++i) {
        glm::vec3 pos1 = vertices[pairs[i].p1];
        glm::vec3 pos2 = vertices[pairs[i].p2];

        float dx = std::abs(pos1.x - pos2.x);
        float dy = std::abs(pos1.y - pos2.y);
        float dz = std::abs(pos1.z - pos2.z);
        float euclidean_dist = glm::length(pos1 - pos2);

        std::cout << "  Pair " << i << ": (" << pairs[i].p1 << ", "
                  << pairs[i].p2 << ")" << std::endl;
        std::cout << "    Axis distances: dx=" << dx << ", dy=" << dy
                  << ", dz=" << dz << std::endl;
        std::cout << "    Euclidean distance: " << euclidean_dist << std::endl;

        // At least one axis distance should be < threshold
        bool has_close_axis = (dx < search_radius) || (dy < search_radius) ||
                              (dz < search_radius);
        EXPECT_TRUE(has_close_axis);

        std::cout << "    Box distance satisfied: "
                  << (has_close_axis ? "YES" : "NO") << std::endl;
    }
}

#endif
