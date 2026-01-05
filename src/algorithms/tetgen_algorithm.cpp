#include "GCore/algorithms/tetgen_algorithm.h"
#include "GCore/Components/MeshComponent.h"
#include "tetgen.h"

RUZINO_NAMESPACE_OPEN_SCOPE

namespace geom_algorithm {

Geometry tetrahedralize(const Geometry& geometry, const TetgenParams& params)
{
    Geometry input_copy = geometry;
    input_copy.apply_transform();
    
    auto mesh_component = input_copy.get_const_component<MeshComponent>();
    if (!mesh_component) {
        throw std::runtime_error("No mesh component found in input geometry");
    }

    // Get mesh data
    const auto& vertices = mesh_component->get_vertices();
    const auto& indices = mesh_component->get_face_vertex_indices();
    const auto& face_counts = mesh_component->get_face_vertex_counts();

    if (vertices.empty() || indices.empty()) {
        throw std::runtime_error("Input mesh is empty");
    }

    // Prepare TetGen input/output - use scope to ensure cleanup
    tetgenio input, output;
    tetgenbehavior behavior;

    try {
        // Set up points
        input.numberofpoints = static_cast<int>(vertices.size());
        input.pointlist = new REAL[input.numberofpoints * 3];

        for (size_t i = 0; i < vertices.size(); ++i) {
            input.pointlist[i * 3] = vertices[i][0];
            input.pointlist[i * 3 + 1] = vertices[i][1];
            input.pointlist[i * 3 + 2] = vertices[i][2];
        }

        // Set up facets (only triangles supported)
        input.numberoffacets = 0;
        for (size_t i = 0; i < face_counts.size(); ++i) {
            if (face_counts[i] == 3) {
                input.numberoffacets++;
            } else {
                throw std::runtime_error(
                    "TetGen only supports triangular faces. Please triangulate the mesh first.");
            }
        }

        input.facetlist = new tetgenio::facet[input.numberoffacets];
        input.facetmarkerlist = new int[input.numberoffacets];

        size_t face_idx = 0;
        size_t vertex_offset = 0;

        for (size_t i = 0; i < face_counts.size(); ++i) {
            if (face_counts[i] == 3) {
                tetgenio::facet& f = input.facetlist[face_idx];
                f.numberofpolygons = 1;
                f.polygonlist = new tetgenio::polygon[1];
                f.numberofholes = 0;
                f.holelist = nullptr;

                tetgenio::polygon& p = f.polygonlist[0];
                p.numberofvertices = 3;
                p.vertexlist = new int[3];

                p.vertexlist[0] = indices[vertex_offset];
                p.vertexlist[1] = indices[vertex_offset + 1];
                p.vertexlist[2] = indices[vertex_offset + 2];

                input.facetmarkerlist[face_idx] = 1;
                face_idx++;
            }
            vertex_offset += face_counts[i];
        }

        // Set TetGen behavior
        behavior.plc = 1;  // Piecewise Linear Complex
        behavior.quality = params.refine ? 1 : 0;
        behavior.minratio = params.quality_ratio;
        behavior.fixedvolume = 1;
        behavior.maxvolume = params.max_volume;
        behavior.quiet = params.quiet ? 1 : 0;
        behavior.facesout = 1;  // Output boundary faces

        if (params.conforming_delaunay) {
            behavior.cdt = 1;  // Conforming Delaunay
        }

        // Run TetGen
        tetrahedralize(&behavior, &input, &output);

        if (output.numberoftetrahedra == 0) {
            throw std::runtime_error("TetGen failed to generate tetrahedra");
        }

        // Create output geometry
        Geometry output_geometry;
        std::shared_ptr<MeshComponent> output_mesh =
            std::make_shared<MeshComponent>(&output_geometry);
        output_geometry.attach_component(output_mesh);

        // Convert output points
        std::vector<glm::vec3> output_points;
        output_points.reserve(output.numberofpoints);

        for (int i = 0; i < output.numberofpoints; ++i) {
            output_points.push_back(glm::vec3(
                output.pointlist[i * 3],
                output.pointlist[i * 3 + 1],
                output.pointlist[i * 3 + 2]));
        }

        // Extract all tetrahedra faces as triangle soup
        std::vector<int> output_indices;
        std::vector<int> output_face_counts;
        
        output_indices.reserve(output.numberoftetrahedra * 12);
        output_face_counts.reserve(output.numberoftetrahedra * 4);

        for (int i = 0; i < output.numberoftetrahedra; ++i) {
            int* tet = &output.tetrahedronlist[i * 4];

            // All 4 faces of the tetrahedron
            // Face 0: (v1, v2, v3) - opposite to v0
            output_face_counts.push_back(3);
            output_indices.push_back(tet[1]);
            output_indices.push_back(tet[2]);
            output_indices.push_back(tet[3]);

            // Face 1: (v0, v3, v2) - opposite to v1
            output_face_counts.push_back(3);
            output_indices.push_back(tet[0]);
            output_indices.push_back(tet[3]);
            output_indices.push_back(tet[2]);

            // Face 2: (v0, v1, v3) - opposite to v2
            output_face_counts.push_back(3);
            output_indices.push_back(tet[0]);
            output_indices.push_back(tet[1]);
            output_indices.push_back(tet[3]);

            // Face 3: (v0, v2, v1) - opposite to v3
            output_face_counts.push_back(3);
            output_indices.push_back(tet[0]);
            output_indices.push_back(tet[2]);
            output_indices.push_back(tet[1]);
        }

        // Set mesh data
        output_mesh->set_vertices(output_points);
        output_mesh->set_face_vertex_indices(output_indices);
        output_mesh->set_face_vertex_counts(output_face_counts);

        // tetgenio destructors will clean up automatically
        return output_geometry;
    }
    catch (const std::exception& e) {
        // tetgenio destructors will still be called for cleanup
        throw std::runtime_error(std::string("TetGen error: ") + e.what());
    }
}

}  // namespace geom_algorithm

RUZINO_NAMESPACE_CLOSE_SCOPE
