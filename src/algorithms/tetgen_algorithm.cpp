#include "GCore/algorithms/tetgen_algorithm.h"

#include <cmath>
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

    // Check if we need to preserve face attributes
    std::vector<float> input_face_attributes;
    bool has_preserve_attribute = false;

    if (!params.preserve_face_attribute.empty()) {
        auto attr = mesh_component->get_face_scalar_quantity(
            params.preserve_face_attribute);
        if (!attr.empty()) {
            if (attr.size() != face_counts.size()) {
                throw std::runtime_error(
                    "Face attribute size mismatch: " +
                    params.preserve_face_attribute);
            }
            input_face_attributes = attr;
            has_preserve_attribute = true;
            spdlog::info(
                "[TetGen] Preserving face attribute: {}",
                params.preserve_face_attribute);
        }
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

                // Set facet marker: use preserved attribute if available
                if (has_preserve_attribute) {
                    // Convert float attribute to integer marker
                    // Scale by 1000 to preserve up to 3 decimal places
                    input.facetmarkerlist[face_idx] = static_cast<int>(
                        std::round(input_face_attributes[i] * 1000.0f));
                }
                else {
                    input.facetmarkerlist[face_idx] = 1;
                }

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
            boundary_markers;  // 1.0 = boundary face, 0.0 = interior face
        std::vector<float> preserved_attributes;  // Inherited from input

        if (output.numberoftrifaces > 0) {
            output_indices.reserve(output.numberoftrifaces * 3);
            output_face_counts.reserve(output.numberoftrifaces);
            boundary_markers.reserve(output.numberoftrifaces);

            if (has_preserve_attribute) {
                preserved_attributes.reserve(output.numberoftrifaces);
            }

            for (int i = 0; i < output.numberoftrifaces; ++i) {
                // Get triangle vertices
                int* tri = &output.trifacelist[i * 3];

                output_face_counts.push_back(3);
                output_indices.push_back(tri[0]);
                output_indices.push_back(tri[1]);
                output_indices.push_back(tri[2]);

                // Get marker from TetGen
                int marker_int = (output.trifacemarkerlist != nullptr)
                                     ? output.trifacemarkerlist[i]
                                     : 0;

                if (has_preserve_attribute) {
                    // Preserved attribute: decode from marker
                    // We need to distinguish between boundary markers and
                    // preserved attributes TetGen sets marker to 0 for interior
                    // faces For boundary faces, it uses the input
                    // facetmarkerlist value
                    if (marker_int == 0) {
                        // Interior face - no preserved attribute
                        preserved_attributes.push_back(0.0f);
                        boundary_markers.push_back(0.0f);
                    }
                    else {
                        // Boundary face - decode preserved attribute
                        preserved_attributes.push_back(marker_int / 1000.0f);
                        boundary_markers.push_back(1.0f);
                    }
                }
                else {
                    // Just use marker as boundary indicator
                    boundary_markers.push_back(marker_int != 0 ? 1.0f : 0.0f);
                }
            }
        }

        // Set mesh data
        output_mesh->set_vertices(output_points);
        output_mesh->set_face_vertex_indices(output_indices);
        output_mesh->set_face_vertex_counts(output_face_counts);

        // Add boundary marker as face scalar quantity (1.0 = boundary, 0.0 =
        // interior)
        output_mesh->add_face_scalar_quantity("boundary", boundary_markers);

        // Add preserved attribute if it was requested
        if (has_preserve_attribute) {
            output_mesh->add_face_scalar_quantity(
                params.preserve_face_attribute, preserved_attributes);
            spdlog::info(
                "[TetGen] Preserved attribute '{}' inherited to {} faces",
                params.preserve_face_attribute,
                preserved_attributes.size());
        }

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
