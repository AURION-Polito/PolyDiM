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

#ifndef __PDETOOLS_MESH_PDE_Mesh_Utilities_HPP
#define __PDETOOLS_MESH_PDE_Mesh_Utilities_HPP

#include "MeshDAOImporterFromCsv.hpp"
#include "MeshFromCsvUtilities.hpp"
#include "MeshMatricesDAO.hpp"
#include "MeshUtilities.hpp"

namespace Polydim
{
namespace PDETools
{
namespace Mesh
{

/// @brief Domain definitions and mesh generation/import helpers for PDE problems.
///
/// This namespace bridges the abstract computational domains used by the examples
/// (@ref PDE_Domain_1D, @ref PDE_Domain_2D, @ref PDE_Domain_3D) to the GeDiM mesh
/// generators and importers. It offers a uniform, dimension-templated way to create a
/// mesh from a domain, import one from file, and precompute the per-cell geometric
/// data required by the local spaces. Optional generators throw at run time when the
/// corresponding third-party library (Triangle, Voro++, TetGen) is not enabled.
namespace PDE_Mesh_Utilities
{

/// @brief One-dimensional computational domain (a segment).
class PDE_Domain_1D final
{
  public:
    Eigen::MatrixXd vertices; ///< Segment endpoints (one per column).
    double length;            ///< Length of the segment.
};

/// @brief Two-dimensional computational domain.
class PDE_Domain_2D final
{
  public:
    /// @brief Shape family of a 2D domain.
    enum class Domain_Shape_Types
    {
        Parallelogram = 0, ///< Parallelogram (enables structured meshes).
        Polygon = 1,       ///< Generic polygon.
        Ellipse = 2,       ///< Ellipse (uses the radius/center/rotation fields).
        Unknown = 3        ///< Unspecified shape.
    };

    Eigen::MatrixXd vertices; ///< Domain vertices (one per column); boundary polygon for polygonal shapes.
    double area;              ///< Area of the domain.

    // Ellipse type
    double radius_1;                ///< First semi-axis (ellipse shape only).
    double radius_2;                ///< Second semi-axis (ellipse shape only).
    Eigen::Vector3d center;         ///< Center of the ellipse (ellipse shape only).
    Eigen::Vector3d rotation_angle; ///< Rotation of the ellipse (ellipse shape only).

    Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D::Domain_Shape_Types shape_type; ///< Domain shape family.
};

/// @brief Space–time domain for a 2D time-dependent problem.
struct PDE_Time_Domain_2D final
{
    std::array<double, 2> time_domain; ///< Time interval \f$[t_0, t_1]\f$.
    Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D spatial_domain;
};

/// @brief Three-dimensional computational domain.
class PDE_Domain_3D final
{
  public:
    enum class Domain_Shape_Types
    {
        Parallelepiped = 0,
        Polygon = 1
    };

    Eigen::MatrixXd vertices;           ///< Domain vertices (one per column).
    Eigen::MatrixXi edges;              ///< Domain edges (vertex-index pairs).
    std::vector<Eigen::MatrixXi> faces; ///< Domain faces (per-face vertex/edge indices).
    double volume;                      ///< Volume of the domain.
    Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_3D::Domain_Shape_Types shape_type; ///< Domain shape family.
};

enum class MeshGenerator_Types_1D
{
    Equispaced = 0,  ///< equispaced mesh
    Minimal = 2,     ///< minimal mesh
    CsvImporter = 1, ///< imported csv mesh
};

enum class MeshGenerator_Types_2D
{
    Triangular = 0,               ///< generated triangular mesh
    Minimal = 1,                  ///< generated minimal mesh
    Polygonal = 2,                ///< generated voronoi polygonal mesh
    OFFImporter = 3,              ///< imported off mesh
    CsvImporter = 4,              ///< imported csv mesh
    Squared = 5,                  ///< squared mesh
    RandomDistorted = 6,          ///< random distorted
    TriangularSimpleImporter = 7, ///< import 2D triangular mesh
    StructuredTriangular = 8,     ///<
    QuadFromTriangular = 9        ///< generate quadrilateral mesh starting fromt triangular
};

enum class MeshGenerator_Types_3D
{
    Tetrahedral = 0, ///< generated tetrahedral mesh
    Minimal = 1,     ///< generated minimal mesh
    Polyhedral = 2,  ///< generated voronoi polyhedral mesh
    OVMImporter = 3, ///< imported ovm mesh
    VtkImporter = 4, ///< imported vtk mesh
    CsvImporter = 5, ///< imported csv mesh
    Cubic = 6        ///< cubic mesh
};

/// @brief Generate a 1D mesh on a segment domain.
///
/// Fills @p mesh with a segment discretization according to @p mesh_type. The
/// @c Minimal generator produces a single-element mesh, while @c Equispaced subdivides
/// the segment into cells whose relative length does not exceed @p max_relative_length.
///
/// @param geometry_utilities   GeDiM geometry helper.
/// @param mesh_utilities       GeDiM mesh helper.
/// @param mesh_type            Generator to use.
/// @param pde_domain           Segment domain.
/// @param max_relative_length  Target cell length relative to the segment (Equispaced only).
/// @param[out] mesh            Resulting mesh.
/// @throws std::runtime_error if @p mesh_type is not a supported 1D generator.
inline void create_mesh_1D(const Gedim::GeometryUtilities &geometry_utilities,
                           const Gedim::MeshUtilities &mesh_utilities,
                           const Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_1D &mesh_type,
                           const Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_1D &pde_domain,
                           const double &max_relative_length,
                           Gedim::MeshMatricesDAO &mesh)
{
    switch (mesh_type)
    {
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_1D::Minimal: {
        const Eigen::Vector3d segment_origin = pde_domain.vertices.col(0);
        const Eigen::Vector3d segment_tangent = pde_domain.vertices.col(1) - segment_origin;
        const std::vector<double> coordinates = {0.0, 1.0};
        mesh_utilities.FillMesh1D(geometry_utilities, segment_origin, segment_tangent, coordinates, mesh);
    }
    break;
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_1D::Equispaced: {
        const Eigen::Vector3d segment_origin = pde_domain.vertices.col(0);
        const Eigen::Vector3d segment_tangent = pde_domain.vertices.col(1) - segment_origin;
        const unsigned int num_points = round(1.0 / max_relative_length) + 1;

        mesh_utilities.FillMesh1D(geometry_utilities,
                                  segment_origin,
                                  segment_tangent,
                                  geometry_utilities.EquispaceCoordinates(num_points, 0.0, 1.0, true),
                                  mesh);
    }
    break;
    default:
        throw std::runtime_error("MeshGenerator_Types_1D " + std::to_string((unsigned int)mesh_type) + " not supported");
    }
}

/// @brief Generate a 2D mesh on a polygonal/parallelogram domain.
///
/// Fills @p mesh according to @p mesh_type: unstructured triangular (Triangle),
/// Voronoi polygonal (Voro++), single-polygon minimal, structured squared/triangular
/// grids, randomly distorted quadrilaterals, or quadrilaterals from a triangular mesh.
/// Structured generators require a @c Parallelogram domain. Cell size is controlled by
/// @p max_relative_area (target cell area relative to the domain area).
///
/// @param geometry_utilities GeDiM geometry helper.
/// @param mesh_utilities     GeDiM mesh helper.
/// @param mesh_type          Generator to use.
/// @param pde_domain         2D domain.
/// @param max_relative_area  Target cell area relative to the domain area.
/// @param[out] mesh          Resulting mesh.
/// @throws std::runtime_error if @p mesh_type is unsupported, a required library
///         (Triangle, Voro++) is disabled, or a structured generator is used on a non-parallelogram domain.
inline void create_mesh_2D(const Gedim::GeometryUtilities &geometry_utilities,
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
        const std::vector<unsigned int> markers = {};
        mesh_utilities.Mesh2DFromPolygon(pde_domain.vertices, markers, markers, mesh);
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
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::StructuredTriangular: {
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

        mesh_utilities.CreateStructuredTriangularMesh(domain_origin,
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
                                                            0.4,
                                                            0.4,
                                                            mesh);
    }
    break;
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::QuadFromTriangular: {
        const double max_cell_area = pde_domain.area * max_relative_area;
        mesh_utilities.CreateQuadrilateralMeshFromTriangularMesh(geometry_utilities, pde_domain.vertices, max_cell_area, mesh);
    }
    break;
    default:
        throw std::runtime_error("MeshGenerator_Types_2D " + std::to_string((unsigned int)mesh_type) + " not supported");
    }
}

/// @brief Generate a 3D mesh on a polyhedral/parallelepiped domain.
///
/// Fills @p mesh according to @p mesh_type: unstructured tetrahedral (TetGen), Voronoi
/// polyhedral (Voro++), single-polyhedron minimal, or structured cubic grid. The cubic
/// generator requires a @c Parallelepiped domain. Cell size is controlled by
/// @p max_relative_volume (target cell volume relative to the domain volume).
///
/// @param geometry_utilities  GeDiM geometry helper.
/// @param mesh_utilities      GeDiM mesh helper.
/// @param mesh_type           Generator to use.
/// @param pde_domain          3D domain (vertices, edges and faces).
/// @param max_relative_volume Target cell volume relative to the domain volume.
/// @param[out] mesh           Resulting mesh.
/// @throws std::runtime_error if @p mesh_type is unsupported, a required library
///         (TetGen, Voro++) is disabled, or the cubic generator is used on a non-parallelepiped domain.
inline void create_mesh_3D(const Gedim::GeometryUtilities &geometry_utilities,
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
        const std::vector<unsigned int> markers = {};
        mesh_utilities.Mesh3DFromPolyhedron(pde_domain.vertices, pde_domain.edges, pde_domain.faces, markers, markers, markers, mesh);
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

/// @brief Import a 1D mesh from file.
///
/// @param mesh_type  Importer to use (currently @c CsvImporter).
/// @param file_path  Path to the mesh folder/file.
/// @param[out] mesh  Resulting mesh.
/// @throws std::runtime_error if @p mesh_type is not a supported 1D importer.
inline void import_mesh_1D(const Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_1D &mesh_type,
                           const std::string &file_path,
                           Gedim::MeshMatricesDAO &mesh)
{
    switch (mesh_type)
    {
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_1D::CsvImporter: {
        Gedim::MeshFromCsvUtilities importerUtilities;
        Gedim::MeshFromCsvUtilities::Configuration meshImporterConfiguration;
        meshImporterConfiguration.Folder = file_path;
        meshImporterConfiguration.Separator = ';';
        Gedim::MeshDAOImporterFromCsv importer(importerUtilities);
        importer.Import(meshImporterConfiguration, mesh);
    }
    break;
    default:
        throw std::runtime_error("MeshGenerator_Types_1D " + std::to_string((unsigned int)mesh_type) + " not supported");
    }
}

/// @brief Import a 2D mesh from file.
///
/// Supports CSV import (GeDiM cell files), Object File Format (OFF), and a simple
/// triangular-mesh importer reading @c Cell0Ds / @c Cell2Ds / @c Cell2DsMarker CSV files.
///
/// @param geometry_utilities GeDiM geometry helper.
/// @param mesh_utilities     GeDiM mesh helper.
/// @param mesh_type          Importer to use.
/// @param file_path          Path to the mesh folder/file.
/// @param[out] mesh          Resulting mesh.
/// @throws std::runtime_error if @p mesh_type is not a supported 2D importer.
inline void import_mesh_2D(const Gedim::GeometryUtilities &geometry_utilities,
                           const Gedim::MeshUtilities &mesh_utilities,
                           const Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D &mesh_type,
                           const std::string &file_path,
                           Gedim::MeshMatricesDAO &mesh)
{
    switch (mesh_type)
    {
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::CsvImporter: {
        Gedim::MeshFromCsvUtilities importerUtilities;
        Gedim::MeshFromCsvUtilities::Configuration meshImporterConfiguration;
        meshImporterConfiguration.Folder = file_path;
        meshImporterConfiguration.Separator = ';';
        Gedim::MeshDAOImporterFromCsv importer(importerUtilities);
        importer.Import(meshImporterConfiguration, mesh);
    }
    break;
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::OFFImporter: {
        mesh_utilities.ImportObjectFileFormat(file_path, mesh);
    }
    break;
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::TriangularSimpleImporter: {
        mesh_utilities.ImportTriangularMesh(geometry_utilities,
                                            file_path + "/Cell0Ds.csv",
                                            file_path + "/Cell2Ds.csv",
                                            file_path + "/Cell2DsMarker.csv",
                                            ',',
                                            mesh);
    }
    break;
    default:
        throw std::runtime_error("MeshGenerator_Types_2D " + std::to_string((unsigned int)mesh_type) + " not supported");
    }
}

/// @brief Import a 3D mesh from file.
///
/// Supports CSV import (GeDiM cell files), OpenVolumeMesh (OVM), and VTK import.
///
/// @param mesh_utilities GeDiM mesh helper.
/// @param mesh_type      Importer to use.
/// @param file_path      Path to the mesh folder/file.
/// @param[out] mesh      Resulting mesh.
/// @throws std::runtime_error if @p mesh_type is not a supported 3D importer.
inline void import_mesh_3D(const Gedim::MeshUtilities &mesh_utilities,
                           const Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_3D &mesh_type,
                           const std::string &file_path,
                           Gedim::MeshMatricesDAO &mesh)
{
    switch (mesh_type)
    {
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_3D::CsvImporter: {
        Gedim::MeshFromCsvUtilities importerUtilities;
        Gedim::MeshFromCsvUtilities::Configuration meshImporterConfiguration;
        meshImporterConfiguration.Folder = file_path;
        meshImporterConfiguration.Separator = ';';
        Gedim::MeshDAOImporterFromCsv importer(importerUtilities);
        importer.Import(meshImporterConfiguration, mesh);
    }
    break;
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_3D::OVMImporter: {
        std::vector<std::vector<bool>> meshCell3DsFacesOrientation;
        mesh_utilities.ImportOpenVolumeMesh(file_path, mesh, meshCell3DsFacesOrientation);
    }
    break;
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_3D::VtkImporter: {
        mesh_utilities.ImportVtkMesh3D(file_path, mesh);
    }
    break;
    default:
        throw std::runtime_error("MeshGenerator_Types_3D " + std::to_string((unsigned int)mesh_type) + " not supported");
    }
}

/// @brief Precompute the per-cell geometric data of a 1D mesh.
///
/// @param geometry_utilities GeDiM geometry helper.
/// @param mesh_utilities     GeDiM mesh helper.
/// @param mesh               Input mesh.
/// @return The 1D geometric data required by the local spaces.
inline Gedim::MeshUtilities::MeshGeometricData1D compute_mesh_1D_geometry_data(const Gedim::GeometryUtilities &geometry_utilities,
                                                                               const Gedim::MeshUtilities &mesh_utilities,
                                                                               const Gedim::MeshMatricesDAO &mesh)
{
    return mesh_utilities.FillMesh1DGeometricData(geometry_utilities, mesh);
}

/// @brief Precompute the per-cell geometric data of a 2D mesh.
///
/// All cells are treated as generic (possibly concave) polygons. The set of geometric
/// quantities to compute is controlled by @p mesh_geometric_data_config (all enabled by default).
///
/// @param geometry_utilities         GeDiM geometry helper.
/// @param mesh_utilities             GeDiM mesh helper.
/// @param mesh                       Input mesh.
/// @param mesh_geometric_data_config Which geometric quantities to precompute.
/// @return The 2D geometric data required by the local spaces.
inline Gedim::MeshUtilities::MeshGeometricData2D compute_mesh_2D_geometry_data(
    const Gedim::GeometryUtilities &geometry_utilities,
    const Gedim::MeshUtilities &mesh_utilities,
    const Gedim::MeshMatricesDAO &mesh,
    const Gedim::MeshUtilities::MeshGeometricData2DConfig &mesh_geometric_data_config =
        Gedim::MeshUtilities::MeshGeometricData2DConfig(true, true, true, true, true, true, true, true, true, true, true, true, true))
{
    std::vector<Gedim::GeometryUtilities::PolygonTypes> cell2Ds_types(mesh.Cell2DTotalNumber(),
                                                                      Gedim::GeometryUtilities::PolygonTypes::Generic_Concave);
    return mesh_utilities.FillMesh2DGeometricData(geometry_utilities, mesh, cell2Ds_types, mesh_geometric_data_config);
}

/// @brief Precompute the per-cell geometric data of a 3D mesh.
///
/// Computes the 2D-cell/3D-cell neighbour connectivity before assembling the geometric data.
///
/// @param geometry_utilities GeDiM geometry helper.
/// @param[in,out] mesh       Input mesh (neighbour connectivity is computed in place).
/// @return The 3D geometric data required by the local spaces.
inline Gedim::MeshUtilities::MeshGeometricData3D compute_mesh_3D_geometry_data(const Gedim::GeometryUtilities &geometry_utilities,
                                                                               Gedim::MeshMatricesDAO &mesh)
{
    Gedim::MeshUtilities mesh_utilities;
    mesh_utilities.ComputeCell2DCell3DNeighbours(mesh);
    return mesh_utilities.FillMesh3DGeometricData(geometry_utilities, mesh);
}
} // namespace PDE_Mesh_Utilities
} // namespace Mesh
} // namespace PDETools
} // namespace Polydim

#endif
