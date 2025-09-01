#include "GCore/Components/MeshComponent.h"

#include "GCore/Components/MeshViews.h"
#include "GCore/GOP.h"

USTC_CG_NAMESPACE_OPEN_SCOPE
MeshComponent::MeshComponent(Geometry* attached_operand)
    : GeometryComponent(attached_operand)
{
}

MeshComponent::~MeshComponent()
{
}

void MeshComponent::apply_transform(const glm::mat4& transform)
{
    auto vertices = get_vertices();
    for (auto& vertex : vertices) {
        glm::vec4 homogeneous = glm::vec4(vertex, 1.0f);
        vertex = glm::vec3(transform * homogeneous);
    }

    auto normals = get_normals();
    if (!normals.empty()) {
        for (auto& normal : normals) {
            // Transform normals with the inverse transpose to preserve
            // orthogonality
            glm::mat3 normalTransform =
                glm::mat3(glm::transpose(glm::inverse(transform)));
            normal = glm::normalize(normalTransform * normal);
        }
        set_normals(normals);
    }

    set_vertices(vertices);
}

std::string MeshComponent::to_string() const
{
    std::ostringstream out;
    // Loop over the faces and vertices and print the data
    out << "Topology component. "
        << "Vertices count " << get_vertices().size()
        << ". Face vertices count " << get_face_vertex_counts().size()
        << ". Face vertex indices " << get_face_vertex_indices().size() << ".";
    return out.str();
}

GeometryComponentHandle MeshComponent::copy(Geometry* operand) const
{
    auto ret = std::make_shared<MeshComponent>(operand);
#if USE_USD_SCRATCH_BUFFER
    copy_prim(mesh.GetPrim(), ret->mesh.GetPrim());
    pxr::UsdGeomImageable(mesh).MakeInvisible();
#else
    ret->set_vertices(this->vertices);
    ret->set_face_vertex_counts(this->faceVertexCounts);
    ret->set_face_vertex_indices(this->faceVertexIndices);
    ret->set_normals(this->normals);
    ret->set_display_color(this->displayColor);
    ret->set_texcoords_array(this->texcoordsArray);
#endif
    return ret;
}

void MeshComponent::append_mesh(const std::shared_ptr<MeshComponent>& mesh)

{
    auto this_vertices = get_vertices();
    auto this_face_vertex_indices = get_face_vertex_indices();

    auto that_vertices = mesh->get_vertices();

    auto that_face_vertex_indices = mesh->get_face_vertex_indices();

    int this_index_offset = this_vertices.size();

    this_vertices.resize(this_vertices.size() + that_vertices.size());
    memcpy(
        this_vertices.data() + this_index_offset,
        that_vertices.data(),
        that_vertices.size() * sizeof(glm::vec3));

    // Append face vertex indices
    for (auto& index : that_face_vertex_indices) {
        this_face_vertex_indices.push_back(index + this_index_offset);
    }

    set_vertices(this_vertices);
    set_face_vertex_indices(this_face_vertex_indices);

    auto this_vertex_counts = get_face_vertex_counts();
    auto this_vertex_counts_size = this_vertex_counts.size();
    auto that_vertex_counts = mesh->get_face_vertex_counts();

    this_vertex_counts.resize(
        this_vertex_counts.size() + that_vertex_counts.size());

    memcpy(
        this_vertex_counts.data() + this_vertex_counts_size,
        that_vertex_counts.data(),
        that_vertex_counts.size() * sizeof(int));

    set_face_vertex_counts(this_vertex_counts);
}

USTC_CG_NAMESPACE_CLOSE_SCOPE
