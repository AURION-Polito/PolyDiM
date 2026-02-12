#ifndef __TEST_ZFEM_PCC_2D_LocalSpace_H
#define __TEST_ZFEM_PCC_2D_LocalSpace_H

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "GeometryUtilities.hpp"
#include "MeshUtilities.hpp"
#include "SphereMeshUtilities.hpp"
#include "VTKUtilities.hpp"
#include "ZFEM_PCC_2D_LocalSpace.hpp"
#include "ZFEM_PCC_2D_ReferenceElement.hpp"
#include "ZFEM_PCC_PerformanceAnalysis.hpp"

namespace Polydim
{
namespace UnitTesting
{
//****************************************************************************
TEST(TEST_ZFEM_PCC_2D_LocalSpace, TEST_ZFEM_PCC_2D_LocalSpace_1_2_3_4_5_6)
{
    const std::string exportFolder = "Export/ZFEM/PCC/"
                                     "TEST_ZFEM_PCC_2D_LocalSpace_1_2_3_4_5_6";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = 1.0e-12;
    geometry_utilities_config.Tolerance1D = 1.0e-14;
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    std::vector<Eigen::MatrixXd> polygons(4);

    polygons[0] = Eigen::MatrixXd::Zero(3, 3); // triangle
    polygons[0] << 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0;

    polygons[1] = Eigen::MatrixXd::Zero(3, 4); // square
    polygons[1] << 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0;

    polygons[2] = Eigen::MatrixXd::Zero(3, 4); // concave
    polygons[2] << 0.0, -1.0, 0.0, 1.0, 0.0, 2.0, -1.0, 2.0, 0.0, 0.0, 0.0, 0.0;

    polygons[3] = Eigen::MatrixXd::Zero(3, 8); // hanging nodes
    polygons[3] << 0.0, 0.2, 0.4, 0.6, 0.8, 0.7, 0.5, -0.2, 0.0, 0.0, 0.0, 0.0, 0.3, 0.4, 0.6, 0.5, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0;

    ZFEM::PCC::ZFEM_PCC_2D_ReferenceElement reference_element;
    ZFEM::PCC::ZFEM_PCC_2D_LocalSpace local_space;
    const ZFEM::PCC::ZFEM_PCC_2D_ReferenceElement_Data reference_element_data_1 = reference_element.Create(1);

    std::vector<unsigned int> num_kernel_points = {3, 4, 4, 5, 3, 3, 5};
    for (unsigned int p = 0; p < polygons.size(); p++)
    {
        Eigen::MatrixXd polygon_vertices = polygons[p];

        Gedim::MeshUtilities mesh_utilities;

        std::vector<unsigned int> vertexMarkers;
        std::vector<unsigned int> edgeMarkers;

        Gedim::MeshMatrices mesh_data;
        Gedim::MeshMatricesDAO mesh(mesh_data);
        mesh_utilities.Mesh2DFromPolygon(polygon_vertices, vertexMarkers, edgeMarkers, mesh);
        std::vector<Gedim::GeometryUtilities::PolygonTypes> cell2Ds_types(mesh.Cell2DTotalNumber(),
                                                                          Gedim::GeometryUtilities::PolygonTypes::Generic_Concave);

        const auto mesh_geometric_data =
            mesh_utilities.FillMesh2DGeometricData(geometry_utilities, mesh, cell2Ds_types, reference_element_data_1.mesh_geometric_data_config);

        for (unsigned int order = 1; order < 7; order++)
        {

            const ZFEM::PCC::ZFEM_PCC_2D_ReferenceElement_Data reference_element_data = reference_element.Create(order);

            {
                Gedim::VTKUtilities exporter;
                exporter.AddPoints(reference_element_data.fem_reference_element_data.DofPositions);
                exporter.Export(exportFolder + "/FEMNodes_order_" + std::to_string(order) + ".vtu");
            }

            unsigned int c = 0;

            const ZFEM::PCC::ZFEM_PCC_2D_Polygon_Geometry polygon = {
                geometry_utilities_config.Tolerance1D,
                geometry_utilities_config.Tolerance2D,
                mesh_geometric_data.Cell2DsVertices.at(c),
                mesh_geometric_data.Cell2DsAreas.at(c),
                mesh_geometric_data.Cell2DsDiameters.at(c),
                mesh_geometric_data.Cell2DsEdgeLengths.at(c),
                mesh_geometric_data.Cell2DsEdgeDirections.at(c),
                mesh_geometric_data.Cell2DsEdgeTangents.at(c),
                mesh_geometric_data.Cell2DsEdgeNormals.at(c),
                mesh_geometric_data.Cell2DsChebyshevCenter.at(c),
                mesh_geometric_data.Cell2DsInRadius.at(c),
                mesh_geometric_data.Cell2DsTriangulationsByChebyshevCenter.at(c)};

            const auto local_space_data = local_space.CreateLocalSpace(reference_element_data, polygon);

            ZFEM::PCC::ZFEM_PCC_PerformanceAnalysis performanceAnalysis;
            const auto monomials = Polydim::Utilities::Monomials_2D();
            ZFEM::PCC::ZFEM_PCC_PerformanceAnalysis_Data analysis_data =
                performanceAnalysis.Compute2D(monomials, reference_element_data.monomials_data, local_space, local_space_data);

            const double consistency_error = analysis_data.ConsistencyError.array().sqrt().matrix().maxCoeff();
            const double consistency_error_derivatives =
                (analysis_data.DerivativesConsistencyError[0] + analysis_data.DerivativesConsistencyError[1])
                    .array()
                    .sqrt()
                    .matrix()
                    .maxCoeff();

            unsigned int el = 231 + p;

            const auto jacobi_svd = local_space_data.Dmatrix.jacobiSvd();
            const auto singular_values = jacobi_svd.singularValues();
            ASSERT_TRUE(jacobi_svd.rank() == std::min(local_space_data.Dmatrix.rows(), local_space_data.Dmatrix.cols()));
            ASSERT_TRUE(singular_values.size() == std::min(local_space_data.Dmatrix.rows(), local_space_data.Dmatrix.cols()));

            ASSERT_TRUE(consistency_error <= 1.0e-8);
            ASSERT_TRUE(consistency_error_derivatives <= 1.0e-8);

#ifdef PRINT_ERROR
            std::cout.precision(4);
            std::cout << std::scientific;
            std::cout << el << "\t" << consistency_error << "\t" << consistency_error_derivatives << std::endl;
#endif

            {
                Gedim::VTKUtilities exporter;
                exporter.AddPoints(local_space_data.VirtualNodes);
                exporter.Export(exportFolder + "/VirtualNodes_" + std::to_string(el) + "_order_" + std::to_string(order) + ".vtu");
            }

            {
                Gedim::VTKUtilities exporter;
                exporter.AddPoints(local_space_data.DOFsCoordinates);
                exporter.Export(exportFolder + "/CoarseNodes_" + std::to_string(el) + "_order_" + std::to_string(order) + ".vtu");
            }

            {
                Gedim::VTKUtilities exporter;
                exporter.AddPoints(polygon.Vertices);
                exporter.AddPolygon(polygon.Vertices);
                exporter.Export(exportFolder + "/Polygon_" + std::to_string(el) + ".vtu");
            }

            {
                Gedim::VTKUtilities exporter;
                exporter.AddPoint(local_space_data.KernelIncenter);

                Gedim::SphereMeshUtilities sphere_utilities(geometry_utilities, mesh_utilities);
                Eigen::MatrixXd points =
                    sphere_utilities.circle(local_space_data.KernelIncenter, local_space_data.KernelInRadius, 100);

                exporter.AddPoints(points);
                exporter.AddPolygon(points);
                exporter.Export(exportFolder + "/InCircle_" + std::to_string(el) + ".vtu");
            }

            for (unsigned int t = 0; t < polygon.Vertices.cols(); t++)
            {
                // Export domain
                {
                    Gedim::VTKUtilities vtkUtilities;
                    vtkUtilities.AddPolygon(local_space_data.fem_geometry[t].Vertices);
                    vtkUtilities.Export(exportFolder + "/Domain_" + std::to_string(el) + "_Triangle_" + std::to_string(t) + ".vtu");
                }
            }

            ASSERT_TRUE((local_space_data.VirtualWeights.colwise().sum().transpose() -
                         Eigen::VectorXd::Ones(local_space_data.VirtualNodes.cols()))
                            .norm() < 1.0e-3);
            ASSERT_TRUE(
                (local_space_data.VirtualNodes - local_space_data.DOFsCoordinates * local_space_data.VirtualWeights).norm() < 1.0e-3);
        }
    }
}
// ****************************************************************************
TEST(Test_ZFEM_PCC, Test_ZFEM_PCC_2D_1_to_6)
{
    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = std::numeric_limits<double>::epsilon();
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);
    Gedim::MeshUtilities mesh_utilities;

    ZFEM::PCC::ZFEM_PCC_2D_ReferenceElement reference_element;
    ZFEM::PCC::ZFEM_PCC_2D_LocalSpace local_space;

    const std::string exportTestFolderParaview = "Export/ZFEM/PCC/"
                                                 "Test_ZFEM_PCC_2D_1_to_6";
    Gedim::Output::CreateFolder(exportTestFolderParaview);

    unsigned int num_vertices = 8;
    Eigen::MatrixXd cell2DVertices = Eigen::MatrixXd::Zero(3, num_vertices);
    cell2DVertices.row(0) << 2.8000000000000004e-03, 3.7200000000000004e-02, 3.6900000000000002e-02, 3.3300000000000003e-02,
        2.5600000000000001e-02, 1.4400000000000000e-02, 6.7000000000000002e-03, 3.0999999999999999e-03;
    cell2DVertices.row(1) << 2.1499999999999998e-02, 2.1499999999999998e-02, 2.9500000000000002e-02, 4.0599999999999997e-02,
        5.0700000000000002e-02, 5.0700000000000002e-02, 4.0599999999999997e-02, 2.9500000000000002e-02;

    unsigned int c = 0;
    Gedim::MeshMatrices mesh_data;
    Gedim::MeshMatricesDAO mesh(mesh_data);
    mesh_utilities.Mesh2DFromPolygon(cell2DVertices,
                                     std::vector<unsigned int>(num_vertices, 1),
                                     std::vector<unsigned int>(num_vertices, 1),
                                     mesh);

    const ZFEM::PCC::ZFEM_PCC_2D_ReferenceElement_Data reference_element_data_1 = reference_element.Create(1);
    Gedim::MeshUtilities::MeshGeometricData2D mesh_geometric_data =
        mesh_utilities.FillMesh2DGeometricData(geometry_utilities,
                                               mesh,
                                               {Gedim::GeometryUtilities::PolygonTypes::Generic_Concave},
                                               reference_element_data_1.mesh_geometric_data_config);

    const ZFEM::PCC::ZFEM_PCC_2D_Polygon_Geometry polygon = {geometry_utilities_config.Tolerance1D,
                                                             geometry_utilities_config.Tolerance2D,
                                                             mesh_geometric_data.Cell2DsVertices.at(c),
                                                             mesh_geometric_data.Cell2DsAreas.at(c),
                                                             mesh_geometric_data.Cell2DsDiameters.at(c),
                                                             mesh_geometric_data.Cell2DsEdgeLengths.at(c),
                                                             mesh_geometric_data.Cell2DsEdgeDirections.at(c),
                                                             mesh_geometric_data.Cell2DsEdgeTangents.at(c),
                                                             mesh_geometric_data.Cell2DsEdgeNormals.at(c),
                                                             mesh_geometric_data.Cell2DsChebyshevCenter.at(c),
                                                             mesh_geometric_data.Cell2DsInRadius.at(c),
                                                             mesh_geometric_data.Cell2DsTriangulationsByChebyshevCenter.at(c)};

    for (unsigned int k = 1; k < 7; k++)
    {

        const ZFEM::PCC::ZFEM_PCC_2D_ReferenceElement_Data reference_element_data = reference_element.Create(k);

        const auto local_space_data = local_space.CreateLocalSpace(reference_element_data, polygon);

        ZFEM::PCC::ZFEM_PCC_PerformanceAnalysis performanceAnalysis;
        const auto monomials = Polydim::Utilities::Monomials_2D();
        ZFEM::PCC::ZFEM_PCC_PerformanceAnalysis_Data analysis_data =
            performanceAnalysis.Compute2D(monomials, reference_element_data.monomials_data, local_space, local_space_data);

        const double consistency_error = analysis_data.ConsistencyError.array().sqrt().matrix().maxCoeff();
        const double consistency_error_derivatives =
            (analysis_data.DerivativesConsistencyError[0] + analysis_data.DerivativesConsistencyError[1])
                .array()
                .sqrt()
                .matrix()
                .maxCoeff();

        ASSERT_TRUE(consistency_error <= 1.0e-10);
        ASSERT_TRUE(consistency_error_derivatives <= 1.0e-10);

#ifdef PRINT_ERROR
        std::cout.precision(4);
        std::cout << std::scientific;
        std::cout << c << "\t" << consistency_error << "\t" << consistency_error_derivatives << std::endl;
#endif

        const auto jacobi_svd = local_space_data.Dmatrix.jacobiSvd();
        const auto singular_values = jacobi_svd.singularValues();
        ASSERT_TRUE(jacobi_svd.rank() == std::min(local_space_data.Dmatrix.rows(), local_space_data.Dmatrix.cols()));
        ASSERT_TRUE(singular_values.size() == std::min(local_space_data.Dmatrix.rows(), local_space_data.Dmatrix.cols()));

        {
            Gedim::VTKUtilities exporter;
            exporter.AddPoints(local_space_data.VirtualNodes);
            exporter.Export(exportTestFolderParaview + "/VirtualNodes_" + std::to_string(c) + "_order_" + std::to_string(k) + ".vtu");
        }

        {
            Gedim::VTKUtilities exporter;
            exporter.AddPoints(local_space_data.DOFsCoordinates);
            exporter.Export(exportTestFolderParaview + "/CoarseNodes_" + std::to_string(c) + "_order_" + std::to_string(k) + ".vtu");
        }

        {
            Gedim::VTKUtilities exporter;
            exporter.AddPolygon(polygon.Vertices);
            exporter.Export(exportTestFolderParaview + "/Polygon_" + std::to_string(c) + ".vtu");
        }

        {
            Gedim::VTKUtilities exporter;
            exporter.AddPoint(local_space_data.KernelIncenter);

            Gedim::SphereMeshUtilities sphere_utilities(geometry_utilities, mesh_utilities);
            Eigen::MatrixXd points =
                sphere_utilities.circle(local_space_data.KernelIncenter, local_space_data.KernelInRadius, 100);

            exporter.AddPoints(points);
            exporter.AddPolygon(points);
            exporter.Export(exportTestFolderParaview + "/InCircle_" + std::to_string(c) + ".vtu");
        }

        for (unsigned int t = 0; t < polygon.Vertices.cols(); t++)
        {
            // Export domain
            {
                Gedim::VTKUtilities vtkUtilities;
                vtkUtilities.AddPolygon(local_space_data.fem_geometry[t].Vertices);
                vtkUtilities.Export(exportTestFolderParaview + "/Domain_" + std::to_string(c) + "_Triangle_" +
                                    std::to_string(t) + ".vtu");
            }
        }
    }
}
//****************************************************************************
TEST(Test_ZFEM_PCC, Test_ZFEM_PCC_2D_Degenerate_1)
{

    /// Create folders
    const std::string exportTestFolderParaview = "Export/ZFEM/PCC/"
                                                 "Test_ZFEM_PCC_2D_Degenerate_1";
    Gedim::Output::CreateFolder(exportTestFolderParaview);

    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = std::numeric_limits<double>::epsilon();
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);
    Gedim::MeshUtilities mesh_utilities;

    ZFEM::PCC::ZFEM_PCC_2D_ReferenceElement reference_element;
    ZFEM::PCC::ZFEM_PCC_2D_LocalSpace local_space;

    unsigned int num_vertices = 5;
    Eigen::MatrixXd cell2DVertices = Eigen::MatrixXd::Zero(3, num_vertices);
    cell2DVertices << 0.0, -1.0, -0.5, 0.0, 1.0, 0.0, 2.0, 0.5, -1.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0;

    Gedim::MeshMatrices mesh_data;
    Gedim::MeshMatricesDAO mesh(mesh_data);
    mesh_utilities.Mesh2DFromPolygon(cell2DVertices,
                                     std::vector<unsigned int>(num_vertices, 1),
                                     std::vector<unsigned int>(num_vertices, 1),
                                     mesh);

    const ZFEM::PCC::ZFEM_PCC_2D_ReferenceElement_Data reference_element_data_1 = reference_element.Create(1);
    Gedim::MeshUtilities::MeshGeometricData2D mesh_geometric_data =
        mesh_utilities.FillMesh2DGeometricData(geometry_utilities,
                                               mesh,
                                               {Gedim::GeometryUtilities::PolygonTypes::Generic_Concave},
                                               reference_element_data_1.mesh_geometric_data_config);

    unsigned int c = 0;
    unsigned int el = 66;
    const ZFEM::PCC::ZFEM_PCC_2D_Polygon_Geometry polygon = {geometry_utilities_config.Tolerance1D,
                                                             geometry_utilities_config.Tolerance2D,
                                                             mesh_geometric_data.Cell2DsVertices.at(c),
                                                             mesh_geometric_data.Cell2DsAreas.at(c),
                                                             mesh_geometric_data.Cell2DsDiameters.at(c),
                                                             mesh_geometric_data.Cell2DsEdgeLengths.at(c),
                                                             mesh_geometric_data.Cell2DsEdgeDirections.at(c),
                                                             mesh_geometric_data.Cell2DsEdgeTangents.at(c),
                                                             mesh_geometric_data.Cell2DsEdgeNormals.at(c),
                                                             mesh_geometric_data.Cell2DsChebyshevCenter.at(c),
                                                             mesh_geometric_data.Cell2DsInRadius.at(c),
                                                             mesh_geometric_data.Cell2DsTriangulationsByChebyshevCenter.at(c)};

    for (unsigned int k = 1; k < 7; k++)
    {

        const ZFEM::PCC::ZFEM_PCC_2D_ReferenceElement_Data reference_element_data = reference_element.Create(k);
        const auto local_space_data = local_space.CreateLocalSpace(reference_element_data, polygon);

        ZFEM::PCC::ZFEM_PCC_PerformanceAnalysis performanceAnalysis;
        const auto monomials = Polydim::Utilities::Monomials_2D();
        ZFEM::PCC::ZFEM_PCC_PerformanceAnalysis_Data analysis_data =
            performanceAnalysis.Compute2D(monomials, reference_element_data.monomials_data, local_space, local_space_data);

        const double consistency_error = analysis_data.ConsistencyError.array().sqrt().matrix().maxCoeff();
        const double consistency_error_derivatives =
            (analysis_data.DerivativesConsistencyError[0] + analysis_data.DerivativesConsistencyError[1])
                .array()
                .sqrt()
                .matrix()
                .maxCoeff();

        ASSERT_TRUE(consistency_error <= 1.0e-10);
        ASSERT_TRUE(consistency_error_derivatives <= 1.0e-8);

#ifdef PRINT_ERROR
        std::cout.precision(4);
        std::cout << std::scientific;
        std::cout << el << "\t" << consistency_error << "\t" << consistency_error_derivatives << std::endl;
#endif

        const auto jacobi_svd = local_space_data.Dmatrix.jacobiSvd();
        const auto singular_values = jacobi_svd.singularValues();
        ASSERT_TRUE(jacobi_svd.rank() == std::min(local_space_data.Dmatrix.rows(), local_space_data.Dmatrix.cols()));
        ASSERT_TRUE(singular_values.size() == std::min(local_space_data.Dmatrix.rows(), local_space_data.Dmatrix.cols()));

        {
            Gedim::VTKUtilities exporter;
            exporter.AddPoints(local_space_data.VirtualNodes);
            exporter.Export(exportTestFolderParaview + "/VirtualNodes_" + std::to_string(el) + "_order_" + std::to_string(k) + ".vtu");
        }

        {
            Gedim::VTKUtilities exporter;
            exporter.AddPoints(local_space_data.DOFsCoordinates);
            exporter.Export(exportTestFolderParaview + "/CoarseNodes_" + std::to_string(el) + "_order_" + std::to_string(k) + ".vtu");
        }

        {
            Gedim::VTKUtilities exporter;
            exporter.AddPolygon(polygon.Vertices);
            exporter.Export(exportTestFolderParaview + "/Polygon_" + std::to_string(el) + ".vtu");
        }

        {
            Gedim::VTKUtilities exporter;
            exporter.AddPoint(local_space_data.KernelIncenter);

            Gedim::SphereMeshUtilities sphere_utilities(geometry_utilities, mesh_utilities);
            Eigen::MatrixXd points =
                sphere_utilities.circle(local_space_data.KernelIncenter, local_space_data.KernelInRadius, 100);

            exporter.AddPoints(points);
            exporter.AddPolygon(points);
            exporter.Export(exportTestFolderParaview + "/InCircle_" + std::to_string(el) + ".vtu");
        }

        for (unsigned int t = 0; t < polygon.Vertices.cols(); t++)
        {
            // Export domain
            {
                Gedim::VTKUtilities vtkUtilities;
                vtkUtilities.AddPolygon(local_space_data.fem_geometry[t].Vertices);
                vtkUtilities.Export(exportTestFolderParaview + "/Domain_" + std::to_string(el) + "_Triangle_" +
                                    std::to_string(t) + ".vtu");
            }
        }
    }
}
//****************************************************************************
TEST(Test_ZFEM_PCC, Test_ZFEM_PCC_2D_Degenerate_2)
{

    /// Create folders
    const std::string exportTestFolderParaview = "Export/ZFEM/PCC/"
                                                 "Test_ZFEM_PCC_2D_Degenerate_2";
    Gedim::Output::CreateFolder(exportTestFolderParaview);

    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = std::numeric_limits<double>::epsilon();
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);
    Gedim::MeshUtilities mesh_utilities;

    ZFEM::PCC::ZFEM_PCC_2D_ReferenceElement reference_element;
    ZFEM::PCC::ZFEM_PCC_2D_LocalSpace local_space;

    unsigned int num_vertices = 5;
    Eigen::MatrixXd cell2DVertices = Eigen::MatrixXd::Zero(3, num_vertices);
    cell2DVertices << -1.7573830359442127e-01, -8.8849979536405832e-02, -1.4347600452097295e-01,
        -2.5885943953862300e-01, 6.6692372082946627e-01, 2.0642475138832439e-01, 1.1236184401433286e-01,
        -1.6545532276760880e-01, -2.6569309905462624e-01, 1.1236184401433121e-01, 0.0, 0.0, 0.0, 0.0, 0.0;

    Gedim::MeshMatrices mesh_data;
    Gedim::MeshMatricesDAO mesh(mesh_data);
    mesh_utilities.Mesh2DFromPolygon(cell2DVertices,
                                     std::vector<unsigned int>(num_vertices, 1),
                                     std::vector<unsigned int>(num_vertices, 1),
                                     mesh);

    const ZFEM::PCC::ZFEM_PCC_2D_ReferenceElement_Data reference_element_data_1 = reference_element.Create(1);
    Gedim::MeshUtilities::MeshGeometricData2D mesh_geometric_data =
        mesh_utilities.FillMesh2DGeometricData(geometry_utilities,
                                               mesh,
                                               {Gedim::GeometryUtilities::PolygonTypes::Generic_Concave},
                                               reference_element_data_1.mesh_geometric_data_config);

    unsigned int c = 0;
    unsigned int el = 67;
    const ZFEM::PCC::ZFEM_PCC_2D_Polygon_Geometry polygon = {geometry_utilities_config.Tolerance1D,
                                                             geometry_utilities_config.Tolerance2D,
                                                             mesh_geometric_data.Cell2DsVertices.at(c),
                                                             mesh_geometric_data.Cell2DsAreas.at(c),
                                                             mesh_geometric_data.Cell2DsDiameters.at(c),
                                                             mesh_geometric_data.Cell2DsEdgeLengths.at(c),
                                                             mesh_geometric_data.Cell2DsEdgeDirections.at(c),
                                                             mesh_geometric_data.Cell2DsEdgeTangents.at(c),
                                                             mesh_geometric_data.Cell2DsEdgeNormals.at(c),
                                                             mesh_geometric_data.Cell2DsChebyshevCenter.at(c),
                                                             mesh_geometric_data.Cell2DsInRadius.at(c),
                                                             mesh_geometric_data.Cell2DsTriangulationsByChebyshevCenter.at(c)};

    for (unsigned int k = 1; k < 7; k++)
    {

        const ZFEM::PCC::ZFEM_PCC_2D_ReferenceElement_Data reference_element_data = reference_element.Create(k);

        const auto local_space_data = local_space.CreateLocalSpace(reference_element_data, polygon);

        ZFEM::PCC::ZFEM_PCC_PerformanceAnalysis performanceAnalysis;
        const auto monomials = Polydim::Utilities::Monomials_2D();
        ZFEM::PCC::ZFEM_PCC_PerformanceAnalysis_Data analysis_data =
            performanceAnalysis.Compute2D(monomials, reference_element_data.monomials_data, local_space, local_space_data);

        const double consistency_error = analysis_data.ConsistencyError.array().sqrt().matrix().maxCoeff();
        const double consistency_error_derivatives =
            (analysis_data.DerivativesConsistencyError[0] + analysis_data.DerivativesConsistencyError[1])
                .array()
                .sqrt()
                .matrix()
                .maxCoeff();

        ASSERT_TRUE(consistency_error <= 1.0e-10);
        ASSERT_TRUE(consistency_error_derivatives <= 1.0e-8);

#ifdef PRINT_ERROR
        std::cout.precision(4);
        std::cout << std::scientific;
        std::cout << el << "\t" << consistency_error << "\t" << consistency_error_derivatives << std::endl;
#endif

        const auto jacobi_svd = local_space_data.Dmatrix.jacobiSvd();
        const auto singular_values = jacobi_svd.singularValues();
        ASSERT_TRUE(jacobi_svd.rank() == std::min(local_space_data.Dmatrix.rows(), local_space_data.Dmatrix.cols()));
        ASSERT_TRUE(singular_values.size() == std::min(local_space_data.Dmatrix.rows(), local_space_data.Dmatrix.cols()));

        {
            Gedim::VTKUtilities exporter;
            exporter.AddPoints(local_space_data.VirtualNodes);
            exporter.Export(exportTestFolderParaview + "/VirtualNodes_" + std::to_string(el) + "_order_" + std::to_string(k) + ".vtu");
        }

        {
            Gedim::VTKUtilities exporter;
            exporter.AddPoints(local_space_data.DOFsCoordinates);
            exporter.Export(exportTestFolderParaview + "/CoarseNodes_" + std::to_string(el) + "_order_" + std::to_string(k) + ".vtu");
        }

        {
            Gedim::VTKUtilities exporter;
            exporter.AddPolygon(polygon.Vertices);
            exporter.Export(exportTestFolderParaview + "/Polygon_" + std::to_string(el) + ".vtu");
        }

        {
            Gedim::VTKUtilities exporter;
            exporter.AddPoint(local_space_data.KernelIncenter);

            Gedim::SphereMeshUtilities sphere_utilities(geometry_utilities, mesh_utilities);
            Eigen::MatrixXd points =
                sphere_utilities.circle(local_space_data.KernelIncenter, local_space_data.KernelInRadius, 100);

            exporter.AddPoints(points);
            exporter.AddPolygon(points);
            exporter.Export(exportTestFolderParaview + "/InCircle_" + std::to_string(el) + ".vtu");
        }

        for (unsigned int t = 0; t < polygon.Vertices.cols(); t++)
        {
            // Export domain
            {
                Gedim::VTKUtilities vtkUtilities;
                vtkUtilities.AddPolygon(local_space_data.fem_geometry[t].Vertices);
                vtkUtilities.Export(exportTestFolderParaview + "/Domain_" + std::to_string(el) + "_Triangle_" +
                                    std::to_string(t) + ".vtu");
            }
        }
    }
}

} // namespace UnitTesting
} // namespace Polydim

#endif
