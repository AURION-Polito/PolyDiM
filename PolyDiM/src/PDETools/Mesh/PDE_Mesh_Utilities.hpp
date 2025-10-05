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
struct PDE_Mesh_Utilities
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
        Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D::Domain_Shape_Types shape_type;
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
        Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_3D::Domain_Shape_Types shape_type;
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
        Squared = 5,     ///< squared mesh
        RandomDistorted = 6
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

    static inline void create_mesh_1D(const Gedim::GeometryUtilities &geometry_utilities,
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

            mesh_utilities.FillMesh1D(geometry_utilities, segment_origin, segment_tangent, {0.0, 1.0}, mesh);
        }
        break;
        case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_1D::Equispaced: {
            const Eigen::Vector3d segment_origin = pde_domain.vertices.col(0);
            const Eigen::Vector3d segment_tangent = pde_domain.vertices.col(1) - segment_origin;

            mesh_utilities.FillMesh1D(geometry_utilities,
                                      segment_origin,
                                      segment_tangent,
                                      geometry_utilities.EquispaceCoordinates(max_relative_length, true),
                                      mesh);
        }
        break;
        default:
            throw std::runtime_error("MeshGenerator_Types_1D " + std::to_string((unsigned int)mesh_type) + " not supported");
        }
    }

    static void create_mesh_2D(const Gedim::GeometryUtilities &geometry_utilities,
                               const Gedim::MeshUtilities &mesh_utilities,
                               const Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D &mesh_type,
                               const Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D &pde_domain,
                               const double &max_relative_area,
                               Gedim::MeshMatricesDAO &mesh);

    static void create_mesh_3D(const Gedim::GeometryUtilities &geometry_utilities,
                               const Gedim::MeshUtilities &mesh_utilities,
                               const Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_3D &mesh_type,
                               const Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_3D &pde_domain,
                               const double &max_relative_volume,
                               Gedim::MeshMatricesDAO &mesh);

    static inline void import_mesh_1D(const Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_1D &mesh_type,
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

    static inline void import_mesh_2D(const Gedim::MeshUtilities &mesh_utilities,
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
        default:
            throw std::runtime_error("MeshGenerator_Types_2D " + std::to_string((unsigned int)mesh_type) + " not supported");
        }
    }

    static inline void import_mesh_3D(const Gedim::MeshUtilities &mesh_utilities,
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

    static inline Gedim::MeshUtilities::MeshGeometricData1D compute_mesh_1D_geometry_data(const Gedim::GeometryUtilities &geometry_utilities,
                                                                                          const Gedim::MeshUtilities &mesh_utilities,
                                                                                          const Gedim::MeshMatricesDAO &mesh)
    {
        return mesh_utilities.FillMesh1DGeometricData(geometry_utilities, mesh);
    }

    static inline Gedim::MeshUtilities::MeshGeometricData2D compute_mesh_2D_geometry_data(const Gedim::GeometryUtilities &geometry_utilities,
                                                                                          const Gedim::MeshUtilities &mesh_utilities,
                                                                                          const Gedim::MeshMatricesDAO &mesh)
    {
        std::vector<Gedim::GeometryUtilities::PolygonTypes> cell2Ds_types(mesh.Cell2DTotalNumber(),
                                                                          Gedim::GeometryUtilities::PolygonTypes::Generic_Concave);
        return mesh_utilities.FillMesh2DGeometricData(geometry_utilities, mesh, cell2Ds_types);
    }

    static inline Gedim::MeshUtilities::MeshGeometricData3D compute_mesh_3D_geometry_data(const Gedim::GeometryUtilities &geometry_utilities,
                                                                                          Gedim::MeshMatricesDAO &mesh)
    {
        Gedim::MeshUtilities mesh_utilities;
        mesh_utilities.ComputeCell2DCell3DNeighbours(mesh);
        return mesh_utilities.FillMesh3DGeometricData(geometry_utilities, mesh);
    }
};
} // namespace Mesh
} // namespace PDETools
} // namespace Polydim

#endif
