#include <gtest/gtest.h>

#include <ctime>
#include <memory>

#include "GCore/Components/MeshComponent.h"
#include "GCore/Components/PointsComponent.h"
#include "GCore/algorithms/intersection.h"

#ifdef GPU_GEOM_ALGORITHM

#include "RHI/rhi.hpp"
#include "glm/ext/matrix_transform.hpp"

using namespace Ruzino;

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
    Geometry point_cloud = Geometry::CreatePoints();
    auto pointsComp = point_cloud.get_component<PointsComponent>();

    std::vector<glm::vec3> vertices;
    // Point 0 and 1 are close (distance = 0.01)
    vertices.push_back(glm::vec3(0.0f, 0.0f, 0.0f));
    vertices.push_back(glm::vec3(0.01f, 0.0f, 0.0f));

    // Point 2 is far from others
    vertices.push_back(glm::vec3(1.0f, 1.0f, 0.0f));

    // Point 3 and 4 are close to each other (distance = 0.015)
    vertices.push_back(glm::vec3(2.0f, 0.0f, 0.0f));
    vertices.push_back(glm::vec3(2.015f, 0.0f, 0.0f));

    pointsComp->set_vertices(vertices);

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
    Geometry point_cloud = Geometry::CreatePoints();
    auto pointsComp = point_cloud.get_component<PointsComponent>();

    std::vector<glm::vec3> vertices;
    // All points are far from each other (> 0.02)
    vertices.push_back(glm::vec3(0.0f, 0.0f, 0.0f));
    vertices.push_back(glm::vec3(1.0f, 0.0f, 0.0f));
    vertices.push_back(glm::vec3(0.0f, 1.0f, 0.0f));
    vertices.push_back(glm::vec3(1.0f, 1.0f, 0.0f));

    pointsComp->set_vertices(vertices);

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
    Geometry point_cloud = Geometry::CreatePoints();
    auto pointsComp = point_cloud.get_component<PointsComponent>();

    std::vector<glm::vec3> vertices;
    float spacing = 0.015f;  // Close enough to be neighbors

    // Create a 3x3 grid
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            vertices.push_back(glm::vec3(i * spacing, j * spacing, 0.0f));
        }
    }

    pointsComp->set_vertices(vertices);

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
    Geometry point_cloud = Geometry::CreatePoints();

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
    // Test Euclidean distance with threshold 0.2
    Geometry point_cloud = Geometry::CreatePoints();
    auto pointsComp = point_cloud.get_component<PointsComponent>();

    // Generate random points
    std::vector<glm::vec3> vertices;
    const int num_points = 1000;
    const float search_radius = 0.2f;

    std::srand(static_cast<unsigned>(std::time(nullptr)));  // Use time as seed
    for (int i = 0; i < num_points; ++i) {
        float x = (std::rand() % 1000) / 1000.0f;  // [0, 1)
        float y = (std::rand() % 1000) / 1000.0f;
        float z = (std::rand() % 1000) / 1000.0f;
        vertices.push_back(glm::vec3(x, y, z));
    }

    pointsComp->set_vertices(vertices);

    // Compute expected pairs on CPU using Euclidean distance
    std::set<std::pair<unsigned, unsigned>> expected_pairs;
    for (unsigned i = 0; i < vertices.size(); ++i) {
        for (unsigned j = i + 1; j < vertices.size(); ++j) {
            glm::vec3 pos1 = vertices[i];
            glm::vec3 pos2 = vertices[j];

            float dist = glm::length(pos1 - pos2);

            // Use Euclidean distance
            if (dist < search_radius) {
                expected_pairs.insert(std::make_pair(i, j));
            }
        }
    }

    // Run GPU implementation
    unsigned pair_count = 0;
    std::vector<PointPairs> pairs =
        FindNeighbors(point_cloud, search_radius, pair_count);

    std::cout << "Euclidean distance test (threshold=" << search_radius
              << "): " << std::endl;
    std::cout << "  CPU found " << expected_pairs.size() << " pairs"
              << std::endl;
    std::cout << "  GPU found " << pair_count << " pairs" << std::endl;

    // Convert GPU results to set for comparison
    std::set<std::pair<unsigned, unsigned>> gpu_pairs;
    for (unsigned i = 0; i < pair_count; ++i) {
        unsigned p1 = std::min(pairs[i].p1, pairs[i].p2);
        unsigned p2 = std::max(pairs[i].p1, pairs[i].p2);
        gpu_pairs.insert(std::make_pair(p1, p2));
    }

    // Verify GPU and CPU results match
    EXPECT_EQ(gpu_pairs.size(), expected_pairs.size());

    // Check all expected pairs are found by GPU
    for (const auto& pair : expected_pairs) {
        if (gpu_pairs.count(pair) == 0) {
            glm::vec3 pos1 = vertices[pair.first];
            glm::vec3 pos2 = vertices[pair.second];
            float dist = glm::length(pos1 - pos2);
            std::cout << "  Missing pair (" << pair.first << ", " << pair.second
                      << "):" << std::endl;
            std::cout << "    pos1=(" << pos1.x << ", " << pos1.y << ", "
                      << pos1.z << ")" << std::endl;
            std::cout << "    pos2=(" << pos2.x << ", " << pos2.y << ", "
                      << pos2.z << ")" << std::endl;
            std::cout << "    distance=" << dist << std::endl;
        }
        EXPECT_TRUE(gpu_pairs.count(pair) > 0)
            << "CPU found pair (" << pair.first << ", " << pair.second
            << ") but GPU didn't";
    }

    // Check GPU didn't find extra pairs
    for (const auto& pair : gpu_pairs) {
        if (expected_pairs.count(pair) == 0) {
            glm::vec3 pos1 = vertices[pair.first];
            glm::vec3 pos2 = vertices[pair.second];
            float dist = glm::length(pos1 - pos2);
            std::cout << "  Extra pair (" << pair.first << ", " << pair.second
                      << "):" << std::endl;
            std::cout << "    pos1=(" << pos1.x << ", " << pos1.y << ", "
                      << pos1.z << ")" << std::endl;
            std::cout << "    pos2=(" << pos2.x << ", " << pos2.y << ", "
                      << pos2.z << ")" << std::endl;
            std::cout << "    distance=" << dist << std::endl;
        }
        EXPECT_TRUE(expected_pairs.count(pair) > 0)
            << "GPU found pair (" << pair.first << ", " << pair.second
            << ") but CPU didn't";
    }
}

TEST_F(IntersectionTests, FindNeighborsBoxDistance_0_1)
{
    // Test Euclidean distance with threshold 0.1
    Geometry point_cloud = Geometry::CreatePoints();
    auto pointsComp = point_cloud.get_component<PointsComponent>();

    // Generate random points
    std::vector<glm::vec3> vertices;
    const int num_points = 1000;
    const float search_radius = 0.1f;

    std::srand(static_cast<unsigned>(std::time(nullptr)));  // Use time as seed
    for (int i = 0; i < num_points; ++i) {
        float x = (std::rand() % 1000) / 1000.0f;  // [0, 1)
        float y = (std::rand() % 1000) / 1000.0f;
        float z = (std::rand() % 1000) / 1000.0f;
        vertices.push_back(glm::vec3(x, y, z));
    }

    pointsComp->set_vertices(vertices);

    // Compute expected pairs on CPU using Euclidean distance
    std::set<std::pair<unsigned, unsigned>> expected_pairs;
    for (unsigned i = 0; i < vertices.size(); ++i) {
        for (unsigned j = i + 1; j < vertices.size(); ++j) {
            glm::vec3 pos1 = vertices[i];
            glm::vec3 pos2 = vertices[j];

            float dist = glm::length(pos1 - pos2);

            // Use Euclidean distance
            if (dist < search_radius) {
                expected_pairs.insert(std::make_pair(i, j));
            }
        }
    }

    // Run GPU implementation
    unsigned pair_count = 0;
    std::vector<PointPairs> pairs =
        FindNeighbors(point_cloud, search_radius, pair_count);

    std::cout << "Euclidean distance test (threshold=" << search_radius
              << "): " << std::endl;
    std::cout << "  CPU found " << expected_pairs.size() << " pairs"
              << std::endl;
    std::cout << "  GPU found " << pair_count << " pairs" << std::endl;

    // Convert GPU results to set for comparison
    std::set<std::pair<unsigned, unsigned>> gpu_pairs;
    for (unsigned i = 0; i < pair_count; ++i) {
        unsigned p1 = std::min(pairs[i].p1, pairs[i].p2);
        unsigned p2 = std::max(pairs[i].p1, pairs[i].p2);
        gpu_pairs.insert(std::make_pair(p1, p2));
    }

    // Verify GPU and CPU results match
    EXPECT_EQ(gpu_pairs.size(), expected_pairs.size());

    // Check all expected pairs are found by GPU
    for (const auto& pair : expected_pairs) {
        if (gpu_pairs.count(pair) == 0) {
            glm::vec3 pos1 = vertices[pair.first];
            glm::vec3 pos2 = vertices[pair.second];
            float dist = glm::length(pos1 - pos2);
            std::cout << "  Missing pair (" << pair.first << ", " << pair.second
                      << "):" << std::endl;
            std::cout << "    pos1=(" << pos1.x << ", " << pos1.y << ", "
                      << pos1.z << ")" << std::endl;
            std::cout << "    pos2=(" << pos2.x << ", " << pos2.y << ", "
                      << pos2.z << ")" << std::endl;
            std::cout << "    distance=" << dist << std::endl;
        }
        EXPECT_TRUE(gpu_pairs.count(pair) > 0)
            << "CPU found pair (" << pair.first << ", " << pair.second
            << ") but GPU didn't";
    }

    // Check GPU didn't find extra pairs
    for (const auto& pair : gpu_pairs) {
        if (expected_pairs.count(pair) == 0) {
            glm::vec3 pos1 = vertices[pair.first];
            glm::vec3 pos2 = vertices[pair.second];
            float dist = glm::length(pos1 - pos2);
            std::cout << "  Extra pair (" << pair.first << ", " << pair.second
                      << "):" << std::endl;
            std::cout << "    pos1=(" << pos1.x << ", " << pos1.y << ", "
                      << pos1.z << ")" << std::endl;
            std::cout << "    pos2=(" << pos2.x << ", " << pos2.y << ", "
                      << pos2.z << ")" << std::endl;
            std::cout << "    distance=" << dist << std::endl;
        }
        EXPECT_TRUE(expected_pairs.count(pair) > 0)
            << "GPU found pair (" << pair.first << ", " << pair.second
            << ") but CPU didn't";
    }
}

#endif
