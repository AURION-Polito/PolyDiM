// _LICENSE_HEADER_
//
// Copyright (C) 2019 - 2025.
// Terms register on the GPL-3.0 license.
//
// This file can be redistributed and/or modified under the license terms.
//
// See top level LICENSE file for more details.
//
// This file can be used citing references in CITATION.cff file.

#include "PDE_Mesh_Utilities.hpp"

#include "Gedim_Macro.hpp"

namespace Polydim
{
namespace PDETools
{
namespace Mesh
{

void PDE_Mesh_Utilities::create_mesh_2D(const Gedim::GeometryUtilities &geometry_utilities,
                                        const Gedim::MeshUtilities &mesh_utilities,
                                        const Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D &mesh_type,
                                        const Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D &pde_domain,
                                        const double &max_relative_area,
                                        Gedim::MeshMatricesDAO &mesh)
{
    switch (mesh_type)
    {
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::Triangular: {
#if ENABLE_TRIANGLE == 0
        throw std::runtime_error("Triangle library not active");
#endif
        const double max_cell_area = pde_domain.area * max_relative_area;
        mesh_utilities.CreateTriangularMesh(pde_domain.vertices, max_cell_area, mesh);
    }
    break;
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::Minimal: {
        mesh_utilities.Mesh2DFromPolygon(pde_domain.vertices, {}, {}, mesh);
    }
    break;
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::Polygonal: {
#if ENABLE_VORO == 0
        throw std::runtime_error("Voro library not active");
#endif
        const unsigned num_cells = static_cast<unsigned int>(std::max(1.0, 1.0 / max_relative_area));

        mesh_utilities.CreatePolygonalMesh(geometry_utilities, pde_domain.vertices, num_cells, 10, mesh, 10);
    }
    break;
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::Squared: {
        switch (pde_domain.shape_type)
        {
        case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D::Domain_Shape_Types::Parallelogram:
            break;
        default:
            throw std::runtime_error("Squared mesh cannot be created");
        }

        const double max_cell_edge = sqrt(pde_domain.area * max_relative_area);

        const Eigen::Vector3d domain_origin = pde_domain.vertices.col(0);
        const Eigen::Vector3d domain_base_tangent = pde_domain.vertices.col(1) - domain_origin;
        const Eigen::Vector3d domain_height_tangent = pde_domain.vertices.rightCols(1) - domain_origin;
        const unsigned int num_cells_base = ceil(domain_base_tangent.norm() / max_cell_edge);
        const unsigned int num_cells_height = ceil(domain_height_tangent.norm() / max_cell_edge);

        mesh_utilities.CreateRectangleMesh(domain_origin,
                                           domain_base_tangent,
                                           domain_height_tangent,
                                           geometry_utilities.EquispaceCoordinates(num_cells_base + 1, 0.0, 1.0, true),
                                           geometry_utilities.EquispaceCoordinates(num_cells_height + 1, 0.0, 1.0, true),
                                           mesh);
    }
    break;
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::RandomDistorted: {
        switch (pde_domain.shape_type)
        {
        case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D::Domain_Shape_Types::Parallelogram:
            break;
        default:
            throw std::runtime_error("Squared mesh cannot be created");
        }

        const double max_cell_edge = sqrt(pde_domain.area * max_relative_area);

        const Eigen::Vector3d domain_origin = pde_domain.vertices.col(0);
        const Eigen::Vector3d domain_base_tangent = pde_domain.vertices.col(1) - domain_origin;
        const Eigen::Vector3d domain_height_tangent = pde_domain.vertices.rightCols(1) - domain_origin;
        const unsigned int num_cells_base = ceil(domain_base_tangent.norm() / max_cell_edge);
        const unsigned int num_cells_height = ceil(domain_height_tangent.norm() / max_cell_edge);

        mesh_utilities.CreateRandomlyDeformedQuadrilaterals(geometry_utilities,
                                                            domain_origin,
                                                            domain_base_tangent,
                                                            domain_height_tangent,
                                                            num_cells_base,
                                                            num_cells_height,
                                                            0.2,
                                                            0.2,
                                                            mesh);
    }
    break;
    default:
        throw std::runtime_error("MeshGenerator_Types_2D " + std::to_string((unsigned int)mesh_type) + " not supported");
    }
}

void PDE_Mesh_Utilities::create_mesh_3D(const Gedim::GeometryUtilities &geometry_utilities,
                                        const Gedim::MeshUtilities &mesh_utilities,
                                        const Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_3D &mesh_type,
                                        const Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_3D &pde_domain,
                                        const double &max_relative_volume,
                                        Gedim::MeshMatricesDAO &mesh)
{
    switch (mesh_type)
    {
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_3D::Tetrahedral: {
#if ENABLE_TETGEN == 0
        throw std::runtime_error("Tetgen library not active");
#endif
        const double max_cell_volume = pde_domain.volume * max_relative_volume;
        mesh_utilities.CreateTetrahedralMesh(pde_domain.vertices, pde_domain.edges, pde_domain.faces, max_cell_volume, mesh);
    }
    break;
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_3D::Minimal: {
        mesh_utilities.Mesh3DFromPolyhedron(pde_domain.vertices, pde_domain.edges, pde_domain.faces, {}, {}, {}, mesh);
    }
    break;
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_3D::Polyhedral: {
#if ENABLE_VORO == 0
        throw std::runtime_error("Voro library not active");
#endif
        const unsigned num_cells = static_cast<unsigned int>(std::max(1.0, 1.0 / max_relative_volume));

        mesh_utilities.CreatePolyhedralMesh(geometry_utilities, pde_domain.vertices, pde_domain.edges, pde_domain.faces, num_cells, 10, mesh, 10);
    }
    break;
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_3D::Cubic: {
        switch (pde_domain.shape_type)
        {
        case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_3D::Domain_Shape_Types::Parallelepiped:
            break;
        default:
            throw std::runtime_error("Cubic mesh cannot be created");
        }

        const double max_cell_edge = std::cbrt(pde_domain.volume * max_relative_volume);

        const Eigen::Vector3d domain_origin = pde_domain.vertices.col(0);
        const Eigen::Vector3d domain_base_tangent = pde_domain.vertices.col(1) - domain_origin;
        const Eigen::Vector3d domain_width_tangent = pde_domain.vertices.col(4) - domain_origin;
        const Eigen::Vector3d domain_heigth_tangent = pde_domain.vertices.col(3) - domain_origin;
        const unsigned int num_cells_base = round(domain_base_tangent.norm() / max_cell_edge);
        const unsigned int num_cells_width = round(domain_width_tangent.norm() / max_cell_edge);
        const unsigned int num_cells_height = round(domain_heigth_tangent.norm() / max_cell_edge);

        mesh_utilities.CreateParallelepipedMesh(domain_origin,
                                                domain_base_tangent,
                                                domain_heigth_tangent,
                                                domain_width_tangent,
                                                geometry_utilities.EquispaceCoordinates(num_cells_base + 1, 0.0, 1.0, true),
                                                geometry_utilities.EquispaceCoordinates(num_cells_height + 1, 0.0, 1.0, true),
                                                geometry_utilities.EquispaceCoordinates(num_cells_width + 1, 0.0, 1.0, true),
                                                mesh);
    }
    break;
    default:
        throw std::runtime_error("MeshGenerator_Types_3D " + std::to_string((unsigned int)mesh_type) + " not supported");
    }
}

} // namespace Mesh
} // namespace PDETools
} // namespace Polydim
