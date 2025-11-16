#include <nanobind/nanobind.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/map.h>
#include <nanobind/ndarray.h>

#include "GCore/GOP.h"
#include "GCore/Components.h"
#include "GCore/Components/MeshComponent.h"
#include "GCore/Components/PointsComponent.h"
#include "GCore/Components/CurveComponent.h"
#include "GCore/Components/XformComponent.h"

// Include glm bindings
#include <glm/glm.hpp>

// Include for meta_any interop
#include <entt/meta/meta.hpp>

namespace nb = nanobind;
using namespace USTC_CG;

// Helper function to convert Python list to glm::vec3
glm::vec3 to_vec3(const nb::tuple& t) {
    if (t.size() != 3) {
        throw std::runtime_error("Expected tuple of size 3 for vec3");
    }
    return glm::vec3(
        nb::cast<float>(t[0]),
        nb::cast<float>(t[1]),
        nb::cast<float>(t[2])
    );
}

// Helper function to convert Python list to glm::vec2
glm::vec2 to_vec2(const nb::tuple& t) {
    if (t.size() != 2) {
        throw std::runtime_error("Expected tuple of size 2 for vec2");
    }
    return glm::vec2(
        nb::cast<float>(t[0]),
        nb::cast<float>(t[1])
    );
}

NB_MODULE(geometry_py, m)
{
    m.doc() = "Ruzino Geometry Python bindings";

    // Bind glm types
    nb::class_<glm::vec3>(m, "vec3")
        .def(nb::init<float, float, float>())
        .def_rw("x", &glm::vec3::x)
        .def_rw("y", &glm::vec3::y)
        .def_rw("z", &glm::vec3::z)
        .def("__repr__", [](const glm::vec3& v) {
            return "vec3(" + std::to_string(v.x) + ", " + 
                   std::to_string(v.y) + ", " + 
                   std::to_string(v.z) + ")";
        })
        .def("__getitem__", [](const glm::vec3& v, size_t i) {
            if (i >= 3) throw std::out_of_range("vec3 index out of range");
            return v[i];
        })
        .def("__setitem__", [](glm::vec3& v, size_t i, float value) {
            if (i >= 3) throw std::out_of_range("vec3 index out of range");
            v[i] = value;
        });

    nb::class_<glm::vec2>(m, "vec2")
        .def(nb::init<float, float>())
        .def_rw("x", &glm::vec2::x)
        .def_rw("y", &glm::vec2::y)
        .def("__repr__", [](const glm::vec2& v) {
            return "vec2(" + std::to_string(v.x) + ", " + 
                   std::to_string(v.y) + ")";
        })
        .def("__getitem__", [](const glm::vec2& v, size_t i) {
            if (i >= 2) throw std::out_of_range("vec2 index out of range");
            return v[i];
        })
        .def("__setitem__", [](glm::vec2& v, size_t i, float value) {
            if (i >= 2) throw std::out_of_range("vec2 index out of range");
            v[i] = value;
        });

    // GeometryComponent base class
    nb::class_<GeometryComponent>(m, "GeometryComponent")
        .def("to_string", &GeometryComponent::to_string)
        .def("hash", &GeometryComponent::hash)
        .def("get_attached_operand", &GeometryComponent::get_attached_operand,
            nb::rv_policy::reference);

    // MeshComponent
    nb::class_<MeshComponent, GeometryComponent>(m, "MeshComponent")
        .def("to_string", &MeshComponent::to_string)
        .def("get_vertices", &MeshComponent::get_vertices)
        .def("get_face_vertex_counts", &MeshComponent::get_face_vertex_counts)
        .def("get_face_vertex_indices", &MeshComponent::get_face_vertex_indices)
        .def("get_normals", &MeshComponent::get_normals)
        .def("get_display_color", &MeshComponent::get_display_color)
        .def("get_texcoords_array", &MeshComponent::get_texcoords_array)
        .def("set_vertices", &MeshComponent::set_vertices)
        .def("set_face_vertex_counts", &MeshComponent::set_face_vertex_counts)
        .def("set_face_vertex_indices", &MeshComponent::set_face_vertex_indices)
        .def("set_normals", &MeshComponent::set_normals)
        .def("set_display_color", &MeshComponent::set_display_color)
        .def("set_texcoords_array", &MeshComponent::set_texcoords_array)
        // Vertex scalar quantities
        .def("get_vertex_scalar_quantity", &MeshComponent::get_vertex_scalar_quantity)
        .def("get_vertex_scalar_quantity_names", &MeshComponent::get_vertex_scalar_quantity_names)
        .def("add_vertex_scalar_quantity", &MeshComponent::add_vertex_scalar_quantity)
        .def("set_vertex_scalar_quantities", &MeshComponent::set_vertex_scalar_quantities)
        // Face scalar quantities
        .def("get_face_scalar_quantity", &MeshComponent::get_face_scalar_quantity)
        .def("get_face_scalar_quantity_names", &MeshComponent::get_face_scalar_quantity_names)
        .def("add_face_scalar_quantity", &MeshComponent::add_face_scalar_quantity)
        .def("set_face_scalar_quantities", &MeshComponent::set_face_scalar_quantities)
        // Vertex color quantities
        .def("get_vertex_color_quantity", &MeshComponent::get_vertex_color_quantity)
        .def("get_vertex_color_quantity_names", &MeshComponent::get_vertex_color_quantity_names)
        .def("add_vertex_color_quantity", &MeshComponent::add_vertex_color_quantity)
        .def("set_vertex_color_quantities", &MeshComponent::set_vertex_color_quantities)
        // Face color quantities
        .def("get_face_color_quantity", &MeshComponent::get_face_color_quantity)
        .def("get_face_color_quantity_names", &MeshComponent::get_face_color_quantity_names)
        .def("add_face_color_quantity", &MeshComponent::add_face_color_quantity)
        .def("set_face_color_quantities", &MeshComponent::set_face_color_quantities)
        // Vertex vector quantities
        .def("get_vertex_vector_quantity", &MeshComponent::get_vertex_vector_quantity)
        .def("get_vertex_vector_quantity_names", &MeshComponent::get_vertex_vector_quantity_names)
        .def("add_vertex_vector_quantity", &MeshComponent::add_vertex_vector_quantity)
        .def("set_vertex_vector_quantities", &MeshComponent::set_vertex_vector_quantities)
        // Face vector quantities
        .def("get_face_vector_quantity", &MeshComponent::get_face_vector_quantity)
        .def("get_face_vector_quantity_names", &MeshComponent::get_face_vector_quantity_names)
        .def("add_face_vector_quantity", &MeshComponent::add_face_vector_quantity)
        .def("set_face_vector_quantities", &MeshComponent::set_face_vector_quantities)
        // Parameterization quantities
        .def("get_face_corner_parameterization_quantity", 
            &MeshComponent::get_face_corner_parameterization_quantity)
        .def("get_face_corner_parameterization_quantity_names", 
            &MeshComponent::get_face_corner_parameterization_quantity_names)
        .def("add_face_corner_parameterization_quantity", 
            &MeshComponent::add_face_corner_parameterization_quantity)
        .def("set_face_corner_parameterization_quantities", 
            &MeshComponent::set_face_corner_parameterization_quantities)
        .def("get_vertex_parameterization_quantity", 
            &MeshComponent::get_vertex_parameterization_quantity)
        .def("get_vertex_parameterization_quantity_names", 
            &MeshComponent::get_vertex_parameterization_quantity_names)
        .def("add_vertex_parameterization_quantity", 
            &MeshComponent::add_vertex_parameterization_quantity)
        .def("set_vertex_parameterization_quantities", 
            &MeshComponent::set_vertex_parameterization_quantities)
        .def("append_mesh", &MeshComponent::append_mesh);

    // PointsComponent
    nb::class_<PointsComponent, GeometryComponent>(m, "PointsComponent")
        .def("to_string", &PointsComponent::to_string)
        .def("get_vertices", &PointsComponent::get_vertices)
        .def("get_display_color", &PointsComponent::get_display_color)
        .def("get_width", &PointsComponent::get_width)
        .def("get_normals", &PointsComponent::get_normals)
        .def("set_vertices", &PointsComponent::set_vertices)
        .def("set_display_color", &PointsComponent::set_display_color)
        .def("set_width", &PointsComponent::set_width)
        .def("set_normals", &PointsComponent::set_normals);

    // CurveComponent
    nb::class_<CurveComponent, GeometryComponent>(m, "CurveComponent")
        .def("to_string", &CurveComponent::to_string)
        .def("get_vertices", &CurveComponent::get_vertices)
        .def("get_width", &CurveComponent::get_width)
        .def("get_widths", &CurveComponent::get_widths)
        .def("get_vert_count", &CurveComponent::get_vert_count)
        .def("get_curve_counts", &CurveComponent::get_curve_counts)
        .def("get_display_color", &CurveComponent::get_display_color)
        .def("get_periodic", &CurveComponent::get_periodic)
        .def("get_curve_normals", &CurveComponent::get_curve_normals)
        .def("get_type", &CurveComponent::get_type)
        .def("set_vertices", &CurveComponent::set_vertices)
        .def("set_width", &CurveComponent::set_width)
        .def("set_widths", &CurveComponent::set_widths)
        .def("set_vert_count", &CurveComponent::set_vert_count)
        .def("set_curve_counts", &CurveComponent::set_curve_counts)
        .def("set_display_color", &CurveComponent::set_display_color)
        .def("set_periodic", &CurveComponent::set_periodic)
        .def("set_curve_normals", &CurveComponent::set_curve_normals)
        .def("set_type", &CurveComponent::set_type);

    // CurveType enum
    nb::enum_<CurveComponent::CurveType>(m, "CurveType")
        .value("Linear", CurveComponent::CurveType::Linear)
        .value("Cubic", CurveComponent::CurveType::Cubic);

    // XformComponent
    nb::class_<XformComponent, GeometryComponent>(m, "XformComponent")
        .def("to_string", &XformComponent::to_string);

    // Geometry class
    nb::class_<Geometry>(m, "Geometry")
        .def(nb::init<>())
        .def(nb::init<const Geometry&>())
        .def("apply_transform", &Geometry::apply_transform)
        .def("to_string", &Geometry::to_string)
        .def("hash", &Geometry::hash)
        .def("attach_component", &Geometry::attach_component)
        .def("detach_component", &Geometry::detach_component)
        .def("get_components", &Geometry::get_components)
        // get_component with type - need Python-friendly interface
        .def("get_mesh_component", [](const Geometry& g, size_t idx) {
            return g.get_component<MeshComponent>(idx);
        }, nb::arg("idx") = 0, "Get MeshComponent at index")
        .def("get_points_component", [](const Geometry& g, size_t idx) {
            return g.get_component<PointsComponent>(idx);
        }, nb::arg("idx") = 0, "Get PointsComponent at index")
        .def("get_curve_component", [](const Geometry& g, size_t idx) {
            return g.get_component<CurveComponent>(idx);
        }, nb::arg("idx") = 0, "Get CurveComponent at index")
        .def("get_xform_component", [](const Geometry& g, size_t idx) {
            return g.get_component<XformComponent>(idx);
        }, nb::arg("idx") = 0, "Get XformComponent at index")
        .def("__eq__", [](const Geometry& a, const Geometry& b) { return a == b; })
        .def("__ne__", [](const Geometry& a, const Geometry& b) { return a != b; })
        .def("__repr__", &Geometry::to_string);

    // Static factory methods
    m.def("CreateMesh", &Geometry::CreateMesh, "Create a mesh geometry");
    m.def("CreatePoints", &Geometry::CreatePoints, "Create a points geometry");
    m.def("CreateCurve", &Geometry::CreateCurve, "Create a curve geometry");
#ifdef GEOM_USD_EXTENSION
    m.def("CreateVolume", &Geometry::CreateVolume, "Create a volume geometry");
#endif

    // Helper functions for creating geometries from Python data
    m.def("create_mesh_from_arrays", 
        [](const std::vector<glm::vec3>& vertices,
           const std::vector<int>& face_vertex_counts,
           const std::vector<int>& face_vertex_indices) {
            auto geom = Geometry::CreateMesh();
            auto mesh = geom.get_component<MeshComponent>();
            if (mesh) {
                mesh->set_vertices(vertices);
                mesh->set_face_vertex_counts(face_vertex_counts);
                mesh->set_face_vertex_indices(face_vertex_indices);
            }
            return geom;
        },
        nb::arg("vertices"),
        nb::arg("face_vertex_counts"),
        nb::arg("face_vertex_indices"),
        "Create a mesh from vertex and face arrays");

    m.def("create_points_from_array",
        [](const std::vector<glm::vec3>& vertices) {
            auto geom = Geometry::CreatePoints();
            auto points = geom.get_component<PointsComponent>();
            if (points) {
                points->set_vertices(vertices);
            }
            return geom;
        },
        nb::arg("vertices"),
        "Create points from vertex array");
    
    // NumPy array support functions for zero-copy data transfer
    // Get vertices as numpy array (zero-copy with ownership management)
    m.def("get_vertices_as_array",
        [](std::shared_ptr<MeshComponent> mesh) -> nb::ndarray<nb::numpy, float, nb::shape<-1, 3>> {
            auto verts = mesh->get_vertices();
            size_t size = verts.size();
            float* data = new float[size * 3];
            for (size_t i = 0; i < size; ++i) {
                data[i*3] = verts[i].x;
                data[i*3+1] = verts[i].y;
                data[i*3+2] = verts[i].z;
            }
            size_t shape[2] = {size, 3};
            nb::capsule owner(data, [](void* p) noexcept { delete[] (float*)p; });
            return nb::ndarray<nb::numpy, float, nb::shape<-1, 3>>(
                data, 2, shape, owner);
        },
        nb::arg("mesh"),
        "Get vertices as numpy array (zero-copy)");
    
    // Set vertices from numpy array
    m.def("set_vertices_from_array",
        [](std::shared_ptr<MeshComponent> mesh, nb::ndarray<float, nb::shape<-1, 3>> arr) {
            size_t n = arr.shape(0);
            std::vector<glm::vec3> vertices(n);
            const float* data = arr.data();
            for (size_t i = 0; i < n; ++i) {
                vertices[i] = glm::vec3(data[i*3], data[i*3+1], data[i*3+2]);
            }
            mesh->set_vertices(vertices);
        },
        nb::arg("mesh"),
        nb::arg("array"),
        "Set vertices from numpy array");
    
    // Get face indices as numpy array (zero-copy)
    m.def("get_face_indices_as_array",
        [](std::shared_ptr<MeshComponent> mesh) -> nb::ndarray<nb::numpy, int> {
            auto indices = mesh->get_face_vertex_indices();
            size_t size = indices.size();
            int* data = new int[size];
            std::copy(indices.begin(), indices.end(), data);
            size_t shape[1] = {size};
            nb::capsule owner(data, [](void* p) noexcept { delete[] (int*)p; });
            return nb::ndarray<nb::numpy, int>(data, 1, shape, owner);
        },
        nb::arg("mesh"),
        "Get face indices as numpy array (zero-copy)");
    
    // Helper function to extract Geometry from entt::meta_any
    // This is useful when getting outputs from node graphs
    m.def("extract_geometry_from_meta_any",
        [](const entt::meta_any& any) -> std::shared_ptr<Geometry> {
            if (!any) {
                throw std::runtime_error("meta_any is empty");
            }
            
            // Try shared_ptr first
            if (auto ptr = any.try_cast<std::shared_ptr<Geometry>>()) {
                return *ptr;
            }
            
            // Try raw pointer
            if (auto ptr = any.try_cast<Geometry*>()) {
                // Wrap in shared_ptr (non-owning)
                return std::shared_ptr<Geometry>(*ptr, [](Geometry*){});
            }
            
            // Try reference
            if (auto ptr = any.try_cast<Geometry>()) {
                // Wrap in shared_ptr (non-owning)
                return std::shared_ptr<Geometry>(const_cast<Geometry*>(&(*ptr)), [](Geometry*){});
            }
            
            throw std::runtime_error("meta_any does not contain a Geometry (type: " + std::string(any.type().info().name()) + ")");
        },
        nb::arg("meta_any"),
        "Extract Geometry from entt::meta_any (from node graph outputs)");
}
