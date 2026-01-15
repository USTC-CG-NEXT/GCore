#include "GCore/algorithms/tetgen_algorithm.h"

#include <iostream>

#include "GCore/Components/MeshComponent.h"
#include "spdlog/spdlog.h"
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
            }
            else {
                throw std::runtime_error(
                    "TetGen only supports triangular faces. Please triangulate "
                    "the mesh first.");
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

        // Set minimum dihedral angle constraint for better quality tetrahedra
        // mindihedral is always applied when quality mesh is enabled
        behavior.mindihedral = params.min_dihedral_angle;

        if (params.conforming_delaunay) {
            behavior.cdt = 1;  // Conforming Delaunay
        }

        // Additional quality improvements
        behavior.varvolume = 0;  // Don't use variable volume constraints
        behavior.diagnose = 0;   // Don't just diagnose, actually fix issues

        // Run TetGen
        tetrahedralize(&behavior, &input, &output);

        if (output.numberoftetrahedra == 0) {
            throw std::runtime_error("TetGen failed to generate tetrahedra");
        }

        spdlog::info(
            "[TetGen] Generated {} tetrahedra from {} input faces",
            output.numberoftetrahedra,
            input.numberoffacets);

        // Create output geometry
        Geometry output_geometry;
        std::shared_ptr<MeshComponent> output_mesh =
            std::make_shared<MeshComponent>(&output_geometry);
        output_geometry.attach_component(output_mesh);

        // Convert output points
        std::vector<glm::vec3> output_points;
        output_points.reserve(output.numberofpoints);

        for (int i = 0; i < output.numberofpoints; ++i) {
            output_points.push_back(
                glm::vec3(
                    output.pointlist[i * 3],
                    output.pointlist[i * 3 + 1],
                    output.pointlist[i * 3 + 2]));
        }

        // Use TetGen's triangle face output (includes boundary and interior
        // faces)
        std::vector<int> output_indices;
        std::vector<int> output_face_counts;
        std::vector<float>
            surface_markers;  // Non-zero = boundary face, zero = interior face

        if (output.numberoftrifaces > 0) {
            output_indices.reserve(output.numberoftrifaces * 3);
            output_face_counts.reserve(output.numberoftrifaces);
            surface_markers.reserve(output.numberoftrifaces);

            for (int i = 0; i < output.numberoftrifaces; ++i) {
                // Get triangle vertices
                int* tri = &output.trifacelist[i * 3];

                output_face_counts.push_back(3);
                output_indices.push_back(tri[0]);
                output_indices.push_back(tri[1]);
                output_indices.push_back(tri[2]);

                // Get surface marker (non-zero = boundary/surface face)
                float marker =
                    (output.trifacemarkerlist != nullptr)
                        ? static_cast<float>(output.trifacemarkerlist[i])
                        : 0.0f;
                surface_markers.push_back(marker);
            }
        }
        else {
            // Fallback: extract all tetrahedra faces as triangle soup
            output_indices.reserve(output.numberoftetrahedra * 12);
            output_face_counts.reserve(output.numberoftetrahedra * 4);
            surface_markers.reserve(output.numberoftetrahedra * 4);

            for (int i = 0; i < output.numberoftetrahedra; ++i) {
                int* tet = &output.tetrahedronlist[i * 4];

                // All 4 faces of the tetrahedron
                // Face 0: (v1, v2, v3) - opposite to v0
                output_face_counts.push_back(3);
                output_indices.push_back(tet[1]);
                output_indices.push_back(tet[2]);
                output_indices.push_back(tet[3]);
                surface_markers.push_back(0.0f);  // Unknown if surface

                // Face 1: (v0, v3, v2) - opposite to v1
                output_face_counts.push_back(3);
                output_indices.push_back(tet[0]);
                output_indices.push_back(tet[3]);
                output_indices.push_back(tet[2]);
                surface_markers.push_back(0.0f);

                // Face 2: (v0, v1, v3) - opposite to v2
                output_face_counts.push_back(3);
                output_indices.push_back(tet[0]);
                output_indices.push_back(tet[1]);
                output_indices.push_back(tet[3]);
                surface_markers.push_back(0.0f);

                // Face 3: (v0, v2, v1) - opposite to v3
                output_face_counts.push_back(3);
                output_indices.push_back(tet[0]);
                output_indices.push_back(tet[2]);
                output_indices.push_back(tet[1]);
                surface_markers.push_back(0.0f);
            }
        }

        // Validate tetrahedra quality
        if (!params.quiet && output.numberoftetrahedra > 0) {
            double min_volume = std::numeric_limits<double>::max();
            double max_volume = 0.0;
            double avg_volume = 0.0;
            int degenerate_count = 0;
            const double VOLUME_THRESHOLD = 1e-10;

            for (int i = 0; i < output.numberoftetrahedra; ++i) {
                int* tet = &output.tetrahedronlist[i * 4];

                // Get vertices
                glm::dvec3 v0(
                    output.pointlist[tet[0] * 3 + 0],
                    output.pointlist[tet[0] * 3 + 1],
                    output.pointlist[tet[0] * 3 + 2]);
                glm::dvec3 v1(
                    output.pointlist[tet[1] * 3 + 0],
                    output.pointlist[tet[1] * 3 + 1],
                    output.pointlist[tet[1] * 3 + 2]);
                glm::dvec3 v2(
                    output.pointlist[tet[2] * 3 + 0],
                    output.pointlist[tet[2] * 3 + 1],
                    output.pointlist[tet[2] * 3 + 2]);
                glm::dvec3 v3(
                    output.pointlist[tet[3] * 3 + 0],
                    output.pointlist[tet[3] * 3 + 1],
                    output.pointlist[tet[3] * 3 + 2]);

                // Calculate volume
                glm::dvec3 e1 = v1 - v0;
                glm::dvec3 e2 = v2 - v0;
                glm::dvec3 e3 = v3 - v0;
                double volume =
                    std::abs(glm::dot(e1, glm::cross(e2, e3))) / 6.0;

                min_volume = std::min(min_volume, volume);
                max_volume = std::max(max_volume, volume);
                avg_volume += volume;

                if (volume < VOLUME_THRESHOLD) {
                    degenerate_count++;
                }
            }

            avg_volume /= output.numberoftetrahedra;

            std::cout << "=== TetGen Quality Report ===" << std::endl;
            std::cout << "Total tetrahedra: " << output.numberoftetrahedra
                      << std::endl;
            std::cout << "Volume statistics:" << std::endl;
            std::cout << "  Min: " << min_volume << std::endl;
            std::cout << "  Max: " << max_volume << std::endl;
            std::cout << "  Avg: " << avg_volume << std::endl;
            std::cout << "  Degenerate count (volume < " << VOLUME_THRESHOLD
                      << "): " << degenerate_count << std::endl;

            if (degenerate_count > 0) {
                std::cout << "WARNING: Found " << degenerate_count
                          << " degenerate tetrahedra! Consider adjusting "
                             "TetGen parameters."
                          << std::endl;
            }
            std::cout << "============================" << std::endl;
        }

        // Set mesh data
        output_mesh->set_vertices(output_points);
        output_mesh->set_face_vertex_indices(output_indices);
        output_mesh->set_face_vertex_counts(output_face_counts);

        // Add surface marker as face scalar quantity
        output_mesh->add_face_scalar_quantity("surf_of_vol", surface_markers);

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
