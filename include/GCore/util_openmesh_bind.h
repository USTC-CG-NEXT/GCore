#pragma once

#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralMesh.hh>
#include <memory>

#include "GCore/GOP.h"

USTC_CG_NAMESPACE_OPEN_SCOPE
using PolyMesh = OpenMesh::PolyMesh_ArrayKernelT<>;
using VolumeMesh = OpenVolumeMesh::GeometricTetrahedralMeshV3d;
using PolyhedralVolumeMesh = OpenVolumeMesh::GeometricPolyhedralMeshV3d;

// Surface mesh binding
GEOMETRY_API std::shared_ptr<PolyMesh> operand_to_openmesh(
    Geometry* mesh_oeprand);

GEOMETRY_API std::shared_ptr<Geometry> openmesh_to_operand(PolyMesh* openmesh);

// Volume mesh binding
GEOMETRY_API std::shared_ptr<VolumeMesh> operand_to_openvolulemesh(
    Geometry* mesh_operand);

GEOMETRY_API std::shared_ptr<Geometry> openvolulemesh_to_operand(
    VolumeMesh* volumemesh);

// Volume mesh binding from triangle faces (reconstructs tetrahedra from all
// faces)
GEOMETRY_API std::shared_ptr<VolumeMesh> operand_to_openvolulemesh_from_faces(
    Geometry* mesh_operand);

// Volume mesh binding with explicit tetrahedral connectivity
GEOMETRY_API std::shared_ptr<VolumeMesh> operand_to_openvolulemesh_with_tets(
    Geometry* mesh_operand,
    const std::vector<std::array<int, 4>>& tetrahedral_indices);

// Polyhedral volume mesh binding
GEOMETRY_API std::shared_ptr<PolyhedralVolumeMesh> operand_to_polyhedralmesh(
    Geometry* mesh_operand);

GEOMETRY_API std::shared_ptr<Geometry> polyhedralmesh_to_operand(
    PolyhedralVolumeMesh* polyhedralmesh);

USTC_CG_NAMESPACE_CLOSE_SCOPE
