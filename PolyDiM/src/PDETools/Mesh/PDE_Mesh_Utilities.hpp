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
namespace PDE_Mesh_Utilities
{
struct PDE_Domain_1D final
{
    Eigen::MatrixXd vertices;
    double length;
};

struct PDE_Domain_2D final
{
    enum struct Domain_Shape_Types
    {
        Parallelogram = 0,
        Polygon = 1
    };

    Eigen::MatrixXd vertices;
    double area;
    Domain_Shape_Types shape_type;
};

struct PDE_Domain_3D final
{
    enum struct Domain_Shape_Types
    {
        Parallelepiped = 0,
        Polygon = 1
    };

    Eigen::MatrixXd vertices;
    Eigen::MatrixXi edges;
    std::vector<Eigen::MatrixXi> faces;
    double volume;
    Domain_Shape_Types shape_type;
};

enum struct MeshGenerator_Types_1D
{
    Equispaced = 0,  ///< equispaced mesh
    Minimal = 2,     ///< minimal mesh
    CsvImporter = 1, ///< imported csv mesh
};

enum struct MeshGenerator_Types_2D
{
    Triangular = 0,  ///< generated triangular mesh
    Minimal = 1,     ///< generated minimal mesh
    Polygonal = 2,   ///< generated voronoi polygonal mesh
    OFFImporter = 3, ///< imported off mesh
    CsvImporter = 4, ///< imported csv mesh
    Squared = 5      ///< squared mesh
};

enum struct MeshGenerator_Types_3D
{
    Tetrahedral = 0, ///< generated tetrahedral mesh
    Minimal = 1,     ///< generated minimal mesh
    Polyhedral = 2,  ///< generated voronoi polyhedral mesh
    OVMImporter = 3, ///< imported ovm mesh
    VtkImporter = 4, ///< imported vtk mesh
    CsvImporter = 5, ///< imported csv mesh
    Cubic = 6        ///< cubic mesh
};

inline void create_mesh_1D(const Gedim::GeometryUtilities &geometry_utilities,
                           const Gedim::MeshUtilities &mesh_utilities,
                           const MeshGenerator_Types_1D &mesh_type,
                           const PDE_Domain_1D &pde_domain,
                           const double &max_relative_length,
                           Gedim::MeshMatricesDAO &mesh)
{
    switch (mesh_type)
    {
    case MeshGenerator_Types_1D::Minimal: {
        const Eigen::Vector3d segment_origin = pde_domain.vertices.col(0);
        const Eigen::Vector3d segment_tangent = pde_domain.vertices.col(1) - segment_origin;

        mesh_utilities.FillMesh1D(geometry_utilities, segment_origin, segment_tangent, {0.0, 1.0}, mesh);
    }
    break;
    case MeshGenerator_Types_1D::Equispaced: {
        const double max_cell_length = pde_domain.length * max_relative_length;
        const Eigen::Vector3d segment_origin = pde_domain.vertices.col(0);
        const Eigen::Vector3d segment_tangent = pde_domain.vertices.col(1) - segment_origin;

        mesh_utilities.FillMesh1D(geometry_utilities,
                                  segment_origin,
                                  segment_tangent,
                                  geometry_utilities.EquispaceCoordinates(max_cell_length, true),
                                  mesh);
    }
    break;
    default:
        throw std::runtime_error("MeshGenerator_Types_1D " + std::to_string((unsigned int)mesh_type) + " not supported");
    }
}

inline void create_mesh_2D(const Gedim::GeometryUtilities &geometry_utilities,
                           const Gedim::MeshUtilities &mesh_utilities,
                           const MeshGenerator_Types_2D &mesh_type,
                           const PDE_Domain_2D &pde_domain,
                           const double &max_relative_area,
                           Gedim::MeshMatricesDAO &mesh)
{
    switch (mesh_type)
    {
    case MeshGenerator_Types_2D::Triangular: {
        const double max_cell_area = pde_domain.area * max_relative_area;
        mesh_utilities.CreateTriangularMesh(pde_domain.vertices, max_cell_area, mesh);
    }
    break;
    case MeshGenerator_Types_2D::Minimal: {
        mesh_utilities.Mesh2DFromPolygon(pde_domain.vertices, {}, {}, mesh);
    }
    break;
    case MeshGenerator_Types_2D::Polygonal: {
        const unsigned num_cells = static_cast<unsigned int>(std::max(1.0, 1.0 / max_relative_area));

        mesh_utilities.CreatePolygonalMesh(geometry_utilities, pde_domain.vertices, num_cells, 10, mesh, 10);
    }
    break;
    case MeshGenerator_Types_2D::Squared: {
        switch (pde_domain.shape_type)
        {
        case PDE_Domain_2D::Domain_Shape_Types::Parallelogram:
            break;
        default:
            throw std::runtime_error("Squared mesh cannot be created");
        }

        const double max_cell_edge = sqrt(pde_domain.area * max_relative_area);

        const Eigen::Vector3d domain_origin = pde_domain.vertices.col(0);
        const Eigen::Vector3d domain_base_tangent = pde_domain.vertices.col(1) - domain_origin;
        const Eigen::Vector3d domain_height_tangent = pde_domain.vertices.rightCols(1) - domain_origin;
        const unsigned num_cells_base = ceil(domain_base_tangent.norm() / max_cell_edge);
        const unsigned num_cells_height = ceil(domain_height_tangent.norm() / max_cell_edge);

        mesh_utilities.CreateRectangleMesh(domain_origin,
                                           domain_base_tangent,
                                           domain_height_tangent,
                                           geometry_utilities.EquispaceCoordinates(num_cells_base + 1, 0.0, 1.0, true),
                                           geometry_utilities.EquispaceCoordinates(num_cells_height + 1, 0.0, 1.0, true),
                                           mesh);
    }
    break;
    default:
        throw std::runtime_error("MeshGenerator_Types_2D " + std::to_string((unsigned int)mesh_type) + " not supported");
    }
}

inline void create_mesh_3D(const Gedim::GeometryUtilities &geometry_utilities,
                           const Gedim::MeshUtilities &mesh_utilities,
                           const MeshGenerator_Types_3D &mesh_type,
                           const PDE_Domain_3D &pde_domain,
                           const double &max_relative_volume,
                           Gedim::MeshMatricesDAO &mesh)
{
    switch (mesh_type)
    {
    case MeshGenerator_Types_3D::Tetrahedral: {
        const double max_cell_volume = pde_domain.volume * max_relative_volume;
        mesh_utilities.CreateTetrahedralMesh(pde_domain.vertices, pde_domain.edges, pde_domain.faces, max_cell_volume, mesh);
    }
    break;
    case MeshGenerator_Types_3D::Minimal: {
        mesh_utilities.Mesh3DFromPolyhedron(pde_domain.vertices, pde_domain.edges, pde_domain.faces, {}, {}, {}, mesh);
    }
    break;
    case MeshGenerator_Types_3D::Polyhedral: {
        const unsigned num_cells = static_cast<unsigned int>(std::max(1.0, 1.0 / max_relative_volume));

        mesh_utilities.CreatePolyhedralMesh(geometry_utilities, pde_domain.vertices, pde_domain.edges, pde_domain.faces, num_cells, 10, mesh, 10);
    }
    break;
    case MeshGenerator_Types_3D::Cubic: {
        switch (pde_domain.shape_type)
        {
        case PDE_Domain_3D::Domain_Shape_Types::Parallelepiped:
            break;
        default:
            throw std::runtime_error("Cubic mesh cannot be created");
        }

        const double max_cell_edge = std::cbrt(pde_domain.volume * max_relative_volume);

        const Eigen::Vector3d domain_origin = pde_domain.vertices.col(0);
        const Eigen::Vector3d domain_base_tangent = pde_domain.vertices.col(1) - domain_origin;
        const Eigen::Vector3d domain_width_tangent = pde_domain.vertices.col(4) - domain_origin;
        const Eigen::Vector3d domain_heigth_tangent = pde_domain.vertices.col(3) - domain_origin;
        const unsigned num_cells_base = ceil(domain_base_tangent.norm() / max_cell_edge);
        const unsigned num_cells_width = ceil(domain_width_tangent.norm() / max_cell_edge);
        const unsigned num_cells_height = ceil(domain_heigth_tangent.norm() / max_cell_edge);

        mesh_utilities.CreateParallelepipedMesh(domain_origin,
                                                domain_base_tangent,
                                                domain_heigth_tangent,
                                                domain_width_tangent,
                                                geometry_utilities.EquispaceCoordinates(num_cells_base + 1, 0.0, 1.0, true),
                                                geometry_utilities.EquispaceCoordinates(num_cells_height + 1, 0.0, 1.0, true),
                                                geometry_utilities.EquispaceCoordinates(num_cells_width + 1, 0.0, 1.0, true),
                                                mesh);

        cout << num_cells_base << " " << num_cells_width << " " << num_cells_height << " " << mesh.Cell3DTotalNumber() << endl;
    }
    break;
    default:
        throw std::runtime_error("MeshGenerator_Types_3D " + std::to_string((unsigned int)mesh_type) + " not supported");
    }
}

inline void import_mesh_1D(const Gedim::GeometryUtilities &geometry_utilities,
                           const Gedim::MeshUtilities &mesh_utilities,
                           const MeshGenerator_Types_1D &mesh_type,
                           const std::string &file_path,
                           Gedim::MeshMatricesDAO &mesh)
{
    switch (mesh_type)
    {
    case MeshGenerator_Types_1D::CsvImporter: {
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

inline void import_mesh_2D(const Gedim::GeometryUtilities &geometry_utilities,
                           const Gedim::MeshUtilities &mesh_utilities,
                           const MeshGenerator_Types_2D &mesh_type,
                           const std::string &file_path,
                           Gedim::MeshMatricesDAO &mesh)
{
    switch (mesh_type)
    {
    case MeshGenerator_Types_2D::CsvImporter: {
        Gedim::MeshFromCsvUtilities importerUtilities;
        Gedim::MeshFromCsvUtilities::Configuration meshImporterConfiguration;
        meshImporterConfiguration.Folder = file_path;
        meshImporterConfiguration.Separator = ';';
        Gedim::MeshDAOImporterFromCsv importer(importerUtilities);
        importer.Import(meshImporterConfiguration, mesh);
    }
    break;
    case MeshGenerator_Types_2D::OFFImporter: {
        mesh_utilities.ImportObjectFileFormat(file_path, mesh);
    }
    break;
    default:
        throw std::runtime_error("MeshGenerator_Types_2D " + std::to_string((unsigned int)mesh_type) + " not supported");
    }
}

inline void import_mesh_3D(const Gedim::GeometryUtilities &geometry_utilities,
                           const Gedim::MeshUtilities &mesh_utilities,
                           const MeshGenerator_Types_3D &mesh_type,
                           const std::string &file_path,
                           Gedim::MeshMatricesDAO &mesh)
{
    switch (mesh_type)
    {
    case MeshGenerator_Types_3D::CsvImporter: {
        Gedim::MeshFromCsvUtilities importerUtilities;
        Gedim::MeshFromCsvUtilities::Configuration meshImporterConfiguration;
        meshImporterConfiguration.Folder = file_path;
        meshImporterConfiguration.Separator = ';';
        Gedim::MeshDAOImporterFromCsv importer(importerUtilities);
        importer.Import(meshImporterConfiguration, mesh);
    }
    break;
    case MeshGenerator_Types_3D::OVMImporter: {
        std::vector<std::vector<bool>> meshCell3DsFacesOrientation;
        mesh_utilities.ImportOpenVolumeMesh(file_path, mesh, meshCell3DsFacesOrientation);
    }
    break;
    case MeshGenerator_Types_3D::VtkImporter: {
        mesh_utilities.ImportVtkMesh3D(file_path, mesh);
    }
    default:
        throw std::runtime_error("MeshGenerator_Types_3D " + std::to_string((unsigned int)mesh_type) + " not supported");
    }
}

inline Gedim::MeshUtilities::MeshGeometricData1D compute_mesh_1D_geometry_data(const Gedim::GeometryUtilities &geometry_utilities,
                                                                               const Gedim::MeshUtilities &mesh_utilities,
                                                                               const Gedim::MeshMatricesDAO &mesh)
{
    return mesh_utilities.FillMesh1DGeometricData(geometry_utilities, mesh);
}

inline Gedim::MeshUtilities::MeshGeometricData2D compute_mesh_2D_geometry_data(const Gedim::GeometryUtilities &geometry_utilities,
                                                                               const Gedim::MeshUtilities &mesh_utilities,
                                                                               const Gedim::MeshMatricesDAO &mesh)
{
    std::vector<Gedim::GeometryUtilities::PolygonTypes> cell2Ds_types(mesh.Cell2DTotalNumber(),
                                                                      Gedim::GeometryUtilities::PolygonTypes::Generic_Concave);
    return mesh_utilities.FillMesh2DGeometricData(geometry_utilities, mesh, cell2Ds_types);
}

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
