#include "GCore/util_openmesh_bind.h"

#include "GCore/Components/MeshComponent.h"
USTC_CG_NAMESPACE_OPEN_SCOPE
std::shared_ptr<PolyMesh> operand_to_openmesh(Geometry* mesh_oeprand)
{
    auto openmesh = std::make_shared<PolyMesh>();
    auto topology = mesh_oeprand->get_component<MeshComponent>();

    // Get mesh data
    const auto& vertices = topology->get_vertices();
    for (const auto& vv : vertices) {
        OpenMesh::Vec3f v;
        v[0] = vv[0];
        v[1] = vv[1];
        v[2] = vv[2];
        openmesh->add_vertex(v);
    }

    auto faceVertexIndices = topology->get_face_vertex_indices();
    auto faceVertexCounts = topology->get_face_vertex_counts();
    auto normals = topology->get_normals();
    auto texcoords = topology->get_texcoords_array();
    auto colors = topology->get_display_color();

    bool hasNormals = !normals.empty();
    bool perVertexNormals = hasNormals && (normals.size() == vertices.size());
    bool hasTexcoords = !texcoords.empty();
    bool perVertexTexcoords =
        hasTexcoords && (texcoords.size() == vertices.size());
    bool hasColors = !colors.empty();
    bool perVertexColors = hasColors && (colors.size() == vertices.size());

    // Request vertex normals, texcoords, and colors only when they exist
    if (hasNormals)
        openmesh->request_vertex_normals();
    if (hasTexcoords)
        openmesh->request_vertex_texcoords2D();
    if (hasColors)
        openmesh->request_vertex_colors();

    int vertexIndex = 0;
    for (int i = 0; i < faceVertexCounts.size(); i++) {
        // Create a vector of vertex handles for the face
        std::vector<PolyMesh::VertexHandle> face_vhandles;
        for (int j = 0; j < faceVertexCounts[i]; j++) {
            int index = faceVertexIndices[vertexIndex];
            // Get the vertex handle from the index
            PolyMesh::VertexHandle vh = openmesh->vertex_handle(index);
            // Add it to the vector
            face_vhandles.push_back(vh);

            // Set normal if available
            if (hasNormals) {
                if (perVertexNormals) {
                    // Use per-vertex normals
                    OpenMesh::Vec3f n(
                        normals[index][0],
                        normals[index][1],
                        normals[index][2]);
                    openmesh->set_normal(vh, n);
                }
                else {
                    // Use per-face-vertex normals
                    OpenMesh::Vec3f n(
                        normals[vertexIndex][0],
                        normals[vertexIndex][1],
                        normals[vertexIndex][2]);
                    openmesh->set_normal(vh, n);
                }
            }

            // Set texcoords if available
            if (hasTexcoords) {
                if (perVertexTexcoords) {
                    // Use per-vertex texcoords
                    OpenMesh::Vec2f t(texcoords[index][0], texcoords[index][1]);
                    openmesh->set_texcoord2D(vh, t);
                }
                else {
                    // Use per-face-vertex texcoords
                    OpenMesh::Vec2f t(
                        texcoords[vertexIndex][0], texcoords[vertexIndex][1]);
                    openmesh->set_texcoord2D(vh, t);
                }
            }

            // Set colors if available
            if (hasColors) {
                if (perVertexColors) {
                    // Use per-vertex colors
                    OpenMesh::Vec3f c(
                        colors[index][0], colors[index][1], colors[index][2]);
                    openmesh->set_color(
                        vh,
                        PolyMesh::Color(
                            static_cast<unsigned char>(c[0] * 255),
                            static_cast<unsigned char>(c[1] * 255),
                            static_cast<unsigned char>(c[2] * 255)));
                }
                else {
                    // Use per-face-vertex colors
                    OpenMesh::Vec3f c(
                        colors[vertexIndex][0],
                        colors[vertexIndex][1],
                        colors[vertexIndex][2]);
                    openmesh->set_color(
                        vh,
                        PolyMesh::Color(
                            static_cast<unsigned char>(c[0] * 255),
                            static_cast<unsigned char>(c[1] * 255),
                            static_cast<unsigned char>(c[2] * 255)));
                }
            }

            vertexIndex++;
        }
        // Add the face to the mesh
        openmesh->add_face(face_vhandles);
    }

    if (hasNormals) {
        // Update or compute vertex normals if necessary
        openmesh->update_normals();
    }

    return openmesh;
}
std::shared_ptr<Geometry> openmesh_to_operand(PolyMesh* openmesh)
{
    auto geometry = std::make_shared<Geometry>();
    std::shared_ptr<MeshComponent> mesh =
        std::make_shared<MeshComponent>(geometry.get());
    geometry->attach_component(mesh);

    std::vector<glm::vec3> points;
    std::vector<int> faceVertexIndices;
    std::vector<int> faceVertexCounts;
    std::vector<glm::vec3> normals;
    std::vector<glm::vec2> texcoords;
    std::vector<glm::vec3> colors;

    bool hasNormals = openmesh->has_vertex_normals();
    bool hasTexcoords = openmesh->has_vertex_texcoords2D();
    bool hasColors = openmesh->has_vertex_colors();

    if (hasNormals && !openmesh->has_vertex_normals())
        openmesh->request_vertex_normals();
    if (hasTexcoords && !openmesh->has_vertex_texcoords2D())
        openmesh->request_vertex_texcoords2D();
    if (hasColors && !openmesh->has_vertex_colors())
        openmesh->request_vertex_colors();

    // Ensure normals are updated
    if (hasNormals)
        openmesh->update_normals();

    // Set the points
    for (const auto& v : openmesh->vertices()) {
        const auto& p = openmesh->point(v);
        points.push_back(glm::vec3(p[0], p[1], p[2]));
    }

    // Determine if we should use per-vertex attributes
    bool usePerVertexNormals = hasNormals;
    bool usePerVertexTexcoords = hasTexcoords;
    bool usePerVertexColors = hasColors;

    // If using per-vertex attributes, populate them now
    if (usePerVertexNormals) {
        for (const auto& v : openmesh->vertices()) {
            const auto& n = openmesh->normal(v);
            normals.push_back(glm::vec3(n[0], n[1], n[2]));
        }
    }

    if (usePerVertexTexcoords) {
        for (const auto& v : openmesh->vertices()) {
            const auto& t = openmesh->texcoord2D(v);
            texcoords.push_back(glm::vec2(t[0], t[1]));
        }
    }

    if (usePerVertexColors) {
        for (const auto& v : openmesh->vertices()) {
            const auto& c = openmesh->color(v);
            colors.push_back(
                glm::vec3(c[0] / 255.0f, c[1] / 255.0f, c[2] / 255.0f));
        }
    }

    // If we need per-face-vertex attributes, prepare for them
    std::vector<glm::vec3> faceVertexNormals;
    std::vector<glm::vec2> faceVertexTexcoords;
    std::vector<glm::vec3> faceVertexColors;
    bool hasFaceVertexNormals = hasNormals && !usePerVertexNormals;
    bool hasFaceVertexTexcoords = hasTexcoords && !usePerVertexTexcoords;
    bool hasFaceVertexColors = hasColors && !usePerVertexColors;

    // Set the topology
    for (const auto& f : openmesh->faces()) {
        size_t count = 0;
        for (const auto& vf : f.vertices()) {
            faceVertexIndices.push_back(vf.idx());
            count += 1;

            // If we need per-face-vertex attributes, collect them here
            if (hasFaceVertexNormals) {
                const auto& n = openmesh->normal(vf);
                faceVertexNormals.push_back(glm::vec3(n[0], n[1], n[2]));
            }

            if (hasFaceVertexTexcoords) {
                const auto& t = openmesh->texcoord2D(vf);
                faceVertexTexcoords.push_back(glm::vec2(t[0], t[1]));
            }

            if (hasFaceVertexColors) {
                const auto& c = openmesh->color(vf);
                faceVertexColors.push_back(
                    glm::vec3(c[0] / 255.0f, c[1] / 255.0f, c[2] / 255.0f));
            }
        }
        faceVertexCounts.push_back(count);
    }

    mesh->set_vertices(points);
    mesh->set_face_vertex_indices(faceVertexIndices);
    mesh->set_face_vertex_counts(faceVertexCounts);

    if (hasNormals) {
        mesh->set_normals(usePerVertexNormals ? normals : faceVertexNormals);
    }

    if (hasTexcoords) {
        mesh->set_texcoords_array(
            usePerVertexTexcoords ? texcoords : faceVertexTexcoords);
    }

    if (hasColors) {
        mesh->set_display_color(usePerVertexColors ? colors : faceVertexColors);
    }

    return geometry;
}

std::shared_ptr<VolumeMesh> operand_to_openvolulemesh(Geometry* mesh_operand)
{
    auto volumemesh = std::make_shared<VolumeMesh>();
    auto topology = mesh_operand->get_component<MeshComponent>();

    // Get mesh data
    const auto& vertices = topology->get_vertices();

    // Add vertices to volume mesh
    for (const auto& vv : vertices) {
        VolumeMesh::PointT point(vv[0], vv[1], vv[2]);
        volumemesh->add_vertex(point);
    }

    auto faceVertexIndices = topology->get_face_vertex_indices();
    auto faceVertexCounts = topology->get_face_vertex_counts();

    // For tetrahedral mesh construction from triangle faces, we need to assume
    // a specific data format or use tetrahedralization algorithms

    // Check what type of data we have
    bool allTriangles = true;
    bool hasTetrahedra = false;

    for (size_t i = 0; i < faceVertexCounts.size(); i++) {
        if (faceVertexCounts[i] != 3)
            allTriangles = false;
        if (faceVertexCounts[i] == 4)
            hasTetrahedra = true;
    }

    if (hasTetrahedra && !allTriangles) {
        // Special case: Data actually represents tetrahedra as 4-vertex "faces"
        size_t vertexIndex = 0;
        for (size_t i = 0; i < faceVertexCounts.size(); i++) {
            if (faceVertexCounts[i] == 4) {
                std::vector<OpenVolumeMesh::VertexHandle> tet_vertices;
                for (int j = 0; j < 4; j++) {
                    int index = faceVertexIndices[vertexIndex + j];
                    tet_vertices.push_back(OpenVolumeMesh::VertexHandle(index));
                }
                volumemesh->add_cell(tet_vertices);
            }
            vertexIndex += faceVertexCounts[i];
        }
    }
    else if (allTriangles) {
        // Standard case: Triangle faces representing surface mesh
        // We need a different approach here:
        // Option 1: Assume this is boundary of tetrahedral mesh and we need
        //           additional tetrahedral connectivity data
        // Option 2: Apply Delaunay tetrahedralization
        // Option 3: Assume specific data encoding

        // For now, we'll assume the MeshComponent actually stores tetrahedral
        // connectivity in a non-standard way, or we create a simple example

        // Simple fallback: if we have exactly 4 vertices, create one
        // tetrahedron
        if (vertices.size() == 4 && faceVertexCounts.size() == 4) {
            // Assume the 4 triangular faces define a single tetrahedron
            std::vector<OpenVolumeMesh::VertexHandle> tet_vertices;
            for (size_t i = 0; i < 4; i++) {
                tet_vertices.push_back(OpenVolumeMesh::VertexHandle(i));
            }
            volumemesh->add_cell(tet_vertices);
        }
        else {
            // TODO: Implement proper tetrahedralization
            // For now, throw an error or warning
            // This requires additional tetrahedral connectivity information
            // that's not available in the standard triangle face representation
        }
    }

    return volumemesh;
}

std::shared_ptr<VolumeMesh> operand_to_openvolulemesh_with_tets(
    Geometry* mesh_operand,
    const std::vector<std::array<int, 4>>& tetrahedral_indices)
{
    auto volumemesh = std::make_shared<VolumeMesh>();
    auto topology = mesh_operand->get_component<MeshComponent>();

    // Get mesh data
    const auto& vertices = topology->get_vertices();

    // Add vertices to volume mesh
    for (const auto& vv : vertices) {
        VolumeMesh::PointT point(vv[0], vv[1], vv[2]);
        volumemesh->add_vertex(point);
    }

    // Add tetrahedra using explicit connectivity
    for (const auto& tet : tetrahedral_indices) {
        std::vector<OpenVolumeMesh::VertexHandle> tet_vertices;
        for (int i = 0; i < 4; i++) {
            tet_vertices.push_back(OpenVolumeMesh::VertexHandle(tet[i]));
        }
        volumemesh->add_cell(tet_vertices);
    }

    return volumemesh;
}

std::shared_ptr<Geometry> openvolulemesh_to_operand(VolumeMesh* volumemesh)
{
    auto geometry = std::make_shared<Geometry>();
    std::shared_ptr<MeshComponent> mesh =
        std::make_shared<MeshComponent>(geometry.get());
    geometry->attach_component(mesh);

    std::vector<glm::vec3> points;
    std::vector<int> cellVertexIndices;
    std::vector<int> cellVertexCounts;

    // Extract vertices
    for (auto v_it = volumemesh->vertices_begin();
         v_it != volumemesh->vertices_end();
         ++v_it) {
        const auto& point = volumemesh->vertex(*v_it);
        points.push_back(glm::vec3(point[0], point[1], point[2]));
    }

    // Extract cells (tetrahedra) - using face interface for compatibility
    for (auto c_it = volumemesh->cells_begin(); c_it != volumemesh->cells_end();
         ++c_it) {
        auto cell_vertices = volumemesh->get_cell_vertices(*c_it);
        cellVertexCounts.push_back(cell_vertices.size());

        for (const auto& cv : cell_vertices) {
            cellVertexIndices.push_back(cv.idx());
        }
    }

    mesh->set_vertices(points);
    // Using face_vertex methods as MeshComponent interface - these represent
    // cells
    mesh->set_face_vertex_indices(cellVertexIndices);
    mesh->set_face_vertex_counts(cellVertexCounts);

    return geometry;
}

std::shared_ptr<PolyhedralVolumeMesh> operand_to_polyhedralmesh(
    Geometry* mesh_operand)
{
    auto polymesh = std::make_shared<PolyhedralVolumeMesh>();
    auto topology = mesh_operand->get_component<MeshComponent>();

    // Get mesh data
    const auto& vertices = topology->get_vertices();

    // Add vertices to polyhedral mesh
    for (const auto& vv : vertices) {
        PolyhedralVolumeMesh::PointT point(vv[0], vv[1], vv[2]);
        polymesh->add_vertex(point);
    }

    // Get face data - these could represent faces of polyhedral cells
    auto faceVertexIndices = topology->get_face_vertex_indices();
    auto faceVertexCounts = topology->get_face_vertex_counts();

    // For PolyhedralMesh, we need to build faces first, then cells
    std::vector<OpenVolumeMesh::FaceHandle> faces;

    size_t vertexIndex = 0;
    for (size_t i = 0; i < faceVertexCounts.size(); i++) {
        std::vector<OpenVolumeMesh::VertexHandle> face_vertices;
        for (int j = 0; j < faceVertexCounts[i]; j++) {
            int index = faceVertexIndices[vertexIndex + j];
            face_vertices.push_back(OpenVolumeMesh::VertexHandle(index));
        }

        // Add face to mesh
        auto face_handle = polymesh->add_face(face_vertices);
        if (face_handle.is_valid()) {
            faces.push_back(face_handle);
        }

        vertexIndex += faceVertexCounts[i];
    }

    // TODO: Add logic to group faces into polyhedral cells
    // For now, we assume each face is part of a single cell
    // This needs more sophisticated logic based on your data structure

    return polymesh;
}

std::shared_ptr<Geometry> polyhedralmesh_to_operand(
    PolyhedralVolumeMesh* polyhedralmesh)
{
    auto geometry = std::make_shared<Geometry>();
    std::shared_ptr<MeshComponent> mesh =
        std::make_shared<MeshComponent>(geometry.get());
    geometry->attach_component(mesh);

    std::vector<glm::vec3> points;
    std::vector<int> faceVertexIndices;
    std::vector<int> faceVertexCounts;

    // Extract vertices
    for (auto v_it = polyhedralmesh->vertices_begin();
         v_it != polyhedralmesh->vertices_end();
         ++v_it) {
        const auto& point = polyhedralmesh->vertex(*v_it);
        points.push_back(glm::vec3(point[0], point[1], point[2]));
    }

    // Extract faces from polyhedral mesh
    for (auto f_it = polyhedralmesh->faces_begin();
         f_it != polyhedralmesh->faces_end();
         ++f_it) {
        std::vector<OpenVolumeMesh::VertexHandle> face_vertices;

        // Get face vertices through face-vertex iteration
        for (auto fv_it = polyhedralmesh->fv_iter(*f_it); fv_it.valid();
             ++fv_it) {
            face_vertices.push_back(*fv_it);
        }

        faceVertexCounts.push_back(face_vertices.size());

        for (const auto& fv : face_vertices) {
            faceVertexIndices.push_back(fv.idx());
        }
    }

    mesh->set_vertices(points);
    mesh->set_face_vertex_indices(faceVertexIndices);
    mesh->set_face_vertex_counts(faceVertexCounts);

    return geometry;
}

std::shared_ptr<VolumeMesh> operand_to_openvolulemesh_from_faces(
    Geometry* mesh_operand)
{
    auto volumemesh = std::make_shared<VolumeMesh>();
    auto topology = mesh_operand->get_component<MeshComponent>();

    // Get mesh data
    const auto& vertices = topology->get_vertices();
    auto faceVertexIndices = topology->get_face_vertex_indices();
    auto faceVertexCounts = topology->get_face_vertex_counts();

    // Add vertices to volume mesh
    for (const auto& vv : vertices) {
        VolumeMesh::PointT point(vv[0], vv[1], vv[2]);
        volumemesh->add_vertex(point);
    }

    // Build triangle faces and reconstruct tetrahedra
    std::vector<std::array<int, 3>> triangles;
    size_t vertexIndex = 0;

    // Extract all triangular faces
    for (size_t i = 0; i < faceVertexCounts.size(); i++) {
        if (faceVertexCounts[i] == 3) {
            std::array<int, 3> triangle;
            for (int j = 0; j < 3; j++) {
                triangle[j] = faceVertexIndices[vertexIndex + j];
            }
            triangles.push_back(triangle);
        }
        vertexIndex += faceVertexCounts[i];
    }

    // Build adjacency information: edge -> adjacent triangles
    std::map<std::pair<int, int>, std::vector<size_t>> edge_to_triangles;

    for (size_t t = 0; t < triangles.size(); t++) {
        const auto& tri = triangles[t];
        // Add each edge of the triangle
        for (int i = 0; i < 3; i++) {
            int v1 = tri[i];
            int v2 = tri[(i + 1) % 3];
            // Ensure consistent edge ordering
            if (v1 > v2)
                std::swap(v1, v2);
            edge_to_triangles[{ v1, v2 }].push_back(t);
        }
    }

    // Find tetrahedra by looking for sets of 4 triangles that share the same 4
    // vertices
    std::map<std::set<int>, std::vector<size_t>> vertex_set_to_triangles;

    for (size_t t = 0; t < triangles.size(); t++) {
        const auto& tri = triangles[t];
        std::set<int> tri_vertices(tri.begin(), tri.end());
        vertex_set_to_triangles[tri_vertices].push_back(t);
    }

    // For each unique set of 4 vertices, find all triangles that use only these
    // vertices
    std::set<std::set<int>> processed_tets;

    for (const auto& tri : triangles) {
        std::set<int> tri_vertices(tri.begin(), tri.end());

        // Look for other triangles that could form a tetrahedron with this one
        for (const auto& other_tri : triangles) {
            std::set<int> combined_vertices(tri_vertices);
            combined_vertices.insert(other_tri.begin(), other_tri.end());

            // If we have exactly 4 vertices, this could be part of a
            // tetrahedron
            if (combined_vertices.size() == 4 &&
                processed_tets.find(combined_vertices) ==
                    processed_tets.end()) {
                // Check if we have exactly 4 triangular faces using these 4
                // vertices
                std::vector<size_t> candidate_triangles;

                for (size_t t = 0; t < triangles.size(); t++) {
                    std::set<int> t_vertices(
                        triangles[t].begin(), triangles[t].end());
                    // Check if this triangle uses only vertices from our set of
                    // 4
                    if (std::includes(
                            combined_vertices.begin(),
                            combined_vertices.end(),
                            t_vertices.begin(),
                            t_vertices.end())) {
                        candidate_triangles.push_back(t);
                    }
                }

                // If we found exactly 4 triangles, we have a tetrahedron
                if (candidate_triangles.size() == 4) {
                    processed_tets.insert(combined_vertices);

                    // Create tetrahedron from the 4 vertices
                    std::vector<OpenVolumeMesh::VertexHandle> tet_vertices;
                    for (int v : combined_vertices) {
                        tet_vertices.push_back(OpenVolumeMesh::VertexHandle(v));
                    }
                    volumemesh->add_cell(tet_vertices);
                }
            }
        }
    }

    return volumemesh;
}

USTC_CG_NAMESPACE_CLOSE_SCOPE
