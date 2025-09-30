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

#ifndef __TEST_VEM_MCC_2D_LocalSpace_H
#define __TEST_VEM_MCC_2D_LocalSpace_H

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <numbers>

#include "GeometryUtilities.hpp"
#include "VEM_MCC_2D_EdgeOrtho_ReferenceElement.hpp"
#include "VEM_MCC_2D_EdgeOrtho_Velocity_LocalSpace.hpp"
#include "VEM_MCC_2D_Ortho_EdgeOrtho_Velocity_LocalSpace.hpp"
#include "VEM_MCC_2D_Ortho_Velocity_LocalSpace.hpp"
#include "VEM_MCC_2D_Partial_Velocity_LocalSpace.hpp"
#include "VEM_MCC_2D_Velocity_LocalSpace.hpp"
#include "VEM_MCC_PerformanceAnalysis.hpp"
#include "VTKUtilities.hpp"

namespace Polydim
{
namespace UnitTesting
{
struct Test_VEM_MCC_2D_Polygon_Geometry final
{
    Eigen::MatrixXd Vertices;
    Eigen::Vector3d Centroid;
    double Measure;
    double Diameter;
    std::vector<Eigen::Matrix3d> TriangulationVertices;
    Eigen::VectorXd EdgesLength;
    std::vector<bool> EdgesDirection;
    Eigen::MatrixXd EdgesTangent;
    Eigen::MatrixXd EdgesNormal;
};

void Test_VEM_MCC_2D_Velocity_Export_Dofs(const Polydim::VEM::MCC::VEM_MCC_2D_Polygon_Geometry &polygon,
                                          const Polydim::VEM::MCC::VEM_MCC_2D_Velocity_ReferenceElement_Data &vem_reference_element,
                                          const std::string &exportVtuFolder)
{
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = polygon.Tolerance1D;
    geometryUtilitiesConfig.Tolerance2D = polygon.Tolerance2D;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    const unsigned int num_vertices = polygon.Vertices.cols();
    const unsigned int num_dofs = num_vertices * vem_reference_element.NumDofs0D +
                                  num_vertices * vem_reference_element.NumDofs1D + vem_reference_element.NumDofs2D;

    Eigen::MatrixXd dofs_coordinate = Eigen::MatrixXd::Zero(3, num_dofs);
    std::vector<double> dof_global_index_values(num_dofs);
    std::vector<double> dof_cell_index_values(num_dofs);
    std::vector<double> dof_dimension_values(num_dofs);

    unsigned int id_dofs = 0;
    if (vem_reference_element.NumDofs0D != 0)
    {
        for (unsigned int c = 0; c < num_vertices; ++c)
        {
            for (unsigned int loc_i = 0; loc_i < vem_reference_element.NumDofs0D; ++loc_i)
            {
                dof_cell_index_values[id_dofs] = c;
                dof_dimension_values[id_dofs] = 0;
                dofs_coordinate.col(id_dofs) = polygon.Vertices.col(c);
                dof_global_index_values[id_dofs] = id_dofs;
                id_dofs++;
            }
        }
    }

    if (vem_reference_element.NumDofs1D != 0)
    {
        for (unsigned int c = 0; c < num_vertices; ++c)
        {
            const std::vector<double> local_edge_coordinates =
                geometryUtilities.EquispaceCoordinates(vem_reference_element.NumDofs1D, 0.0, 1.0, false);
            const Eigen::Vector3d edge_origin = polygon.Vertices.col(c);
            const Eigen::Vector3d edge_tangent = polygon.Vertices.col((c + 1) % num_vertices) - edge_origin;

            for (unsigned int loc_i = 0; loc_i < vem_reference_element.NumDofs1D; ++loc_i)
            {
                dof_cell_index_values[id_dofs] = c;
                dof_dimension_values[id_dofs] = 1;
                dofs_coordinate.col(id_dofs) = edge_origin + local_edge_coordinates[loc_i] * edge_tangent;
                dof_global_index_values[id_dofs] = id_dofs;
                id_dofs++;
            }
        }
    }

    if (vem_reference_element.NumDofs2D != 0)
    {

        const auto local_polygon_coordinates =
            geometryUtilities.EquispaceCoordinates(vem_reference_element.NumDofs2D + 1, 0.0, 1.0, true);
        const Eigen::Vector3d polygon_centroid = polygon.Centroid;
        const auto polygonCentroidEdgesDistance =
            geometryUtilities.PolygonCentroidEdgesDistance(polygon.Vertices, polygon.Centroid, polygon.EdgesNormal);

        double circle_diameter = 0.0;
        if (vem_reference_element.NumDofs2D > 1)
            circle_diameter = 0.5 * geometryUtilities.PolygonInRadius(polygonCentroidEdgesDistance);

        for (unsigned int loc_i = 0; loc_i < vem_reference_element.NumDofs2D; ++loc_i)
        {
            dof_cell_index_values[id_dofs] = 0;
            dof_dimension_values[id_dofs] = 2;
            dofs_coordinate.col(id_dofs) =
                polygon_centroid +
                circle_diameter * Eigen::Vector3d(cos(2.0 * std::numbers::pi * local_polygon_coordinates.at(loc_i)),
                                                  sin(2.0 * std::numbers::pi * local_polygon_coordinates.at(loc_i)),
                                                  0.0);
            dof_global_index_values[id_dofs] = id_dofs;
            id_dofs++;
        }
    }

    {
        Gedim::VTKUtilities exporter;
        exporter.AddPoints(dofs_coordinate,
                           {{"cell_dimension",
                             Gedim::VTPProperty::Formats::Points,
                             static_cast<unsigned int>(dof_dimension_values.size()),
                             dof_dimension_values.data()},
                            {"cell_index",
                             Gedim::VTPProperty::Formats::Points,
                             static_cast<unsigned int>(dof_cell_index_values.size()),
                             dof_cell_index_values.data()},
                            {"dof_global_index",
                             Gedim::VTPProperty::Formats::Points,
                             static_cast<unsigned int>(dof_global_index_values.size()),
                             dof_global_index_values.data()}});

        exporter.Export(exportVtuFolder + "/velocity_dofs_" + std::to_string(vem_reference_element.Order) + ".vtu");
    }
}

void Test_VEM_MCC_2D_Pressure_Export_Dofs(const Polydim::VEM::MCC::VEM_MCC_2D_Polygon_Geometry &polygon,
                                          const Polydim::VEM::MCC::VEM_MCC_2D_Pressure_ReferenceElement_Data &vem_reference_element,
                                          const std::string &exportVtuFolder)
{
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = polygon.Tolerance1D;
    geometryUtilitiesConfig.Tolerance2D = polygon.Tolerance2D;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    const unsigned int num_vertices = polygon.Vertices.cols();
    const unsigned int num_dofs = num_vertices * vem_reference_element.NumDofs0D +
                                  num_vertices * vem_reference_element.NumDofs1D + vem_reference_element.NumDofs2D;

    Eigen::MatrixXd dofs_coordinate = Eigen::MatrixXd::Zero(3, num_dofs);
    std::vector<double> dof_global_index_values(num_dofs);
    std::vector<double> dof_cell_index_values(num_dofs);
    std::vector<double> dof_dimension_values(num_dofs);

    unsigned int id_dofs = 0;
    if (vem_reference_element.NumDofs0D != 0)
    {
        for (unsigned int c = 0; c < num_vertices; ++c)
        {
            for (unsigned int loc_i = 0; loc_i < vem_reference_element.NumDofs0D; ++loc_i)
            {
                dof_cell_index_values[id_dofs] = c;
                dof_dimension_values[id_dofs] = 0;
                dofs_coordinate.col(id_dofs) = polygon.Vertices.col(c);
                dof_global_index_values[id_dofs] = id_dofs;
                id_dofs++;
            }
        }
    }

    if (vem_reference_element.NumDofs1D != 0)
    {
        for (unsigned int c = 0; c < num_vertices; ++c)
        {
            const std::vector<double> local_edge_coordinates =
                geometryUtilities.EquispaceCoordinates(vem_reference_element.NumDofs1D, 0.0, 1.0, false);
            const Eigen::Vector3d edge_origin = polygon.Vertices.col(c);
            const Eigen::Vector3d edge_tangent = polygon.Vertices.col((c + 1) % num_vertices) - edge_origin;

            for (unsigned int loc_i = 0; loc_i < vem_reference_element.NumDofs1D; ++loc_i)
            {
                dof_cell_index_values[id_dofs] = c;
                dof_dimension_values[id_dofs] = 1;
                dofs_coordinate.col(id_dofs) = edge_origin + local_edge_coordinates[loc_i] * edge_tangent;
                dof_global_index_values[id_dofs] = id_dofs;
                id_dofs++;
            }
        }
    }

    if (vem_reference_element.NumDofs2D != 0)
    {

        const auto local_polygon_coordinates =
            geometryUtilities.EquispaceCoordinates(vem_reference_element.NumDofs2D + 1, 0.0, 1.0, true);
        const Eigen::Vector3d polygon_centroid = polygon.Centroid;
        const auto polygonCentroidEdgesDistance =
            geometryUtilities.PolygonCentroidEdgesDistance(polygon.Vertices, polygon.Centroid, polygon.EdgesNormal);

        double circle_diameter = 0.0;
        if (vem_reference_element.NumDofs2D > 1)
            circle_diameter = 0.2 * geometryUtilities.PolygonInRadius(polygonCentroidEdgesDistance);

        for (unsigned int loc_i = 0; loc_i < vem_reference_element.NumDofs2D; ++loc_i)
        {
            dof_cell_index_values[id_dofs] = 0;
            dof_dimension_values[id_dofs] = 2;
            dofs_coordinate.col(id_dofs) =
                polygon_centroid +
                circle_diameter * Eigen::Vector3d(cos(2.0 * std::numbers::pi * local_polygon_coordinates.at(loc_i)),
                                                  sin(2.0 * std::numbers::pi * local_polygon_coordinates.at(loc_i)),
                                                  0.0);
            dof_global_index_values[id_dofs] = id_dofs;
            id_dofs++;
        }
    }

    {
        Gedim::VTKUtilities exporter;
        exporter.AddPoints(dofs_coordinate,
                           {{"cell_dimension",
                             Gedim::VTPProperty::Formats::Points,
                             static_cast<unsigned int>(dof_dimension_values.size()),
                             dof_dimension_values.data()},
                            {"cell_index",
                             Gedim::VTPProperty::Formats::Points,
                             static_cast<unsigned int>(dof_cell_index_values.size()),
                             dof_cell_index_values.data()},
                            {"dof_global_index",
                             Gedim::VTPProperty::Formats::Points,
                             static_cast<unsigned int>(dof_global_index_values.size()),
                             dof_global_index_values.data()}});

        exporter.Export(exportVtuFolder + "/pressure_dofs_" + std::to_string(vem_reference_element.Order) + ".vtu");
    }
}

Test_VEM_MCC_2D_Polygon_Geometry Test_VEM_MCC_2D_Geometry(const Gedim::GeometryUtilities &geometry_utilities)
{
    Test_VEM_MCC_2D_Polygon_Geometry result;

    Eigen::MatrixXd &polygonVertices = result.Vertices;
    polygonVertices.setZero(3, 8);
    polygonVertices.row(0) << 2.8000000000000004e-03, 3.7200000000000004e-02, 3.6900000000000002e-02, 3.3300000000000003e-02,
        2.5600000000000001e-02, 1.4400000000000000e-02, 6.7000000000000002e-03, 3.0999999999999999e-03;
    polygonVertices.row(1) << 2.1499999999999998e-02, 2.1499999999999998e-02, 2.9500000000000002e-02, 4.0599999999999997e-02,
        5.0700000000000002e-02, 5.0700000000000002e-02, 4.0599999999999997e-02, 2.9500000000000002e-02;

    result.Measure = 7.9890999999999990e-04;
    result.Diameter = 3.7046997179258676e-02;
    result.Centroid = Eigen::Vector3d(2.0000000000000004e-02, 3.4061346918509809e-02, 0.0);
    result.EdgesLength.resize(8);
    result.EdgesLength << 3.4400000000000000e-02, 8.0056230238501787e-03, 1.1669190203266030e-02, 1.2700393694685222e-02,
        1.1200000000000002e-02, 1.2700393694685220e-02, 1.1669190203266030e-02, 8.0056230238501769e-03;

    result.EdgesDirection.resize(8, true);

    Eigen::MatrixXd &edgeTangents = result.EdgesTangent;
    edgeTangents.setZero(3, 8);
    edgeTangents.row(0) << 3.4400000000000000e-02, -3.0000000000000165e-04, -3.5999999999999990e-03, -7.7000000000000020e-03,
        -1.1200000000000002e-02, -7.6999999999999994e-03, -3.6000000000000003e-03, -2.9999999999999949e-04;
    edgeTangents.row(1) << 0.0000000000000000e+00, 8.0000000000000036e-03, 1.1099999999999995e-02, 1.0100000000000005e-02,
        0.0000000000000000e+00, -1.0100000000000005e-02, -1.1099999999999995e-02, -8.0000000000000036e-03;

    Eigen::MatrixXd &edgeNormals = result.EdgesNormal;
    edgeNormals.setZero(3, 8);
    edgeNormals.row(0) << 0.0000000000000000e+00, 9.9929761570918052e-01, 9.5122281894876248e-01, 7.9525093810490199e-01,
        0.0000000000000000e+00, -7.9525093810490211e-01, -9.5122281894876248e-01, -9.9929761570918074e-01;
    edgeNormals.row(1) << -1.0000000000000000e+00, 3.7473660589094460e-02, 3.0850469803743652e-01, 6.0628041815918254e-01,
        1.0000000000000000e+00, 6.0628041815918243e-01, 3.0850469803743663e-01, 3.7473660589094196e-02;

    // Map triangle points
    std::vector<unsigned int> polygonTriangulation = geometry_utilities.PolygonTriangulationByFirstVertex(polygonVertices);
    result.TriangulationVertices = geometry_utilities.ExtractTriangulationPoints(polygonVertices, polygonTriangulation);

    return result;
}

TEST(Test_VEM_MCC, Test_VEM_MCC_2D_O0_O1_O2_O3)
{
    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = std::numeric_limits<double>::epsilon();
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    const std::string exportFolder = "VEM/MCC/Test_VEM_MCC_2D_O0_O1_O2_O3";
    Gedim::Output::CreateFolder(exportFolder);

    const auto polygon_data = Test_VEM_MCC_2D_Geometry(geometry_utilities);

    Polydim::VEM::MCC::VEM_MCC_2D_Polygon_Geometry polygon = {geometry_utilities_config.Tolerance1D,
                                                              geometry_utilities_config.Tolerance2D,
                                                              polygon_data.Vertices,
                                                              polygon_data.Centroid,
                                                              polygon_data.Measure,
                                                              polygon_data.Diameter,
                                                              polygon_data.TriangulationVertices,
                                                              polygon_data.EdgesLength,
                                                              polygon_data.EdgesDirection,
                                                              polygon_data.EdgesTangent,
                                                              polygon_data.EdgesNormal};

    for (unsigned int k = 0; k < 4; k++)
    {
        Polydim::VEM::MCC::VEM_MCC_2D_Velocity_ReferenceElement vem_reference_element;
        Polydim::VEM::MCC::VEM_MCC_2D_Pressure_ReferenceElement vem_pressure_reference_element;
        Polydim::VEM::MCC::VEM_MCC_2D_Velocity_LocalSpace vem_local_space;

        const auto reference_element_data = vem_reference_element.Create(k);
        const auto pressure_reference_element_data = vem_pressure_reference_element.Create(k);
        const auto local_space = vem_local_space.CreateLocalSpace(reference_element_data, polygon);

        // Export domain
        {
            Gedim::VTKUtilities vtkUtilities;
            vtkUtilities.AddPolygon(polygon_data.Vertices);
            vtkUtilities.Export(exportFolder + "/Polygon.vtu");

            Test_VEM_MCC_2D_Velocity_Export_Dofs(polygon, reference_element_data, exportFolder);
            Test_VEM_MCC_2D_Pressure_Export_Dofs(polygon, pressure_reference_element_data, exportFolder);
        }

        // Test VEM performances
        Polydim::VEM::MCC::VEM_MCC_PerformanceAnalysis performanceAnalysis;

        const auto result =
            performanceAnalysis.Compute(Polydim::Utilities::Monomials_2D(), reference_element_data.MonomialsKp1, vem_local_space, local_space);

        ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(1.0e-10, result.ErrorPi0k, geometry_utilities.Tolerance1D()));

        ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(1.0e-10, result.ErrorGBD, geometry_utilities.Tolerance1D()));

        ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(1.0e-10, result.ErrorStabilization, geometry_utilities.Tolerance1D()));
    }
}

TEST(Test_VEM_MCC, Test_VEM_MCC_Partial_2D_O0_O1_O2_O3)
{
    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = std::numeric_limits<double>::epsilon();
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    const auto polygon_data = Test_VEM_MCC_2D_Geometry(geometry_utilities);

    Polydim::VEM::MCC::VEM_MCC_2D_Polygon_Geometry polygon = {geometry_utilities_config.Tolerance1D,
                                                              geometry_utilities_config.Tolerance2D,
                                                              polygon_data.Vertices,
                                                              polygon_data.Centroid,
                                                              polygon_data.Measure,
                                                              polygon_data.Diameter,
                                                              polygon_data.TriangulationVertices,
                                                              polygon_data.EdgesLength,
                                                              polygon_data.EdgesDirection,
                                                              polygon_data.EdgesTangent,
                                                              polygon_data.EdgesNormal};

    for (unsigned int k = 0; k < 4; k++)
    {
        Polydim::VEM::MCC::VEM_MCC_2D_Velocity_ReferenceElement vem_reference_element;
        Polydim::VEM::MCC::VEM_MCC_2D_Partial_Velocity_LocalSpace vem_local_space;

        const auto reference_element_data = vem_reference_element.Create(k);
        const auto local_space = vem_local_space.CreateLocalSpace(reference_element_data, polygon);

        // Test VEM performances
        Polydim::VEM::MCC::VEM_MCC_PerformanceAnalysis performanceAnalysis;

        const auto result =
            performanceAnalysis.Compute(Polydim::Utilities::Monomials_2D(), reference_element_data.MonomialsKp1, vem_local_space, local_space);

        ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(1.0e-10, result.ErrorPi0k, geometry_utilities.Tolerance1D()));

        ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(1.0e-10, result.ErrorGBD, geometry_utilities.Tolerance1D()));

        ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(1.0e-10, result.ErrorStabilization, geometry_utilities.Tolerance1D()));
    }
}

TEST(Test_VEM_MCC, Test_VEM_MCC_Ortho_2D_O0_O1_O2_O3)
{
    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = std::numeric_limits<double>::epsilon();
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    const auto polygon_data = Test_VEM_MCC_2D_Geometry(geometry_utilities);

    Polydim::VEM::MCC::VEM_MCC_2D_Polygon_Geometry polygon = {geometry_utilities_config.Tolerance1D,
                                                              geometry_utilities_config.Tolerance2D,
                                                              polygon_data.Vertices,
                                                              polygon_data.Centroid,
                                                              polygon_data.Measure,
                                                              polygon_data.Diameter,
                                                              polygon_data.TriangulationVertices,
                                                              polygon_data.EdgesLength,
                                                              polygon_data.EdgesDirection,
                                                              polygon_data.EdgesTangent,
                                                              polygon_data.EdgesNormal};

    for (unsigned int k = 0; k < 4; k++)
    {
        Polydim::VEM::MCC::VEM_MCC_2D_Velocity_ReferenceElement vem_reference_element;
        Polydim::VEM::MCC::VEM_MCC_2D_Ortho_Velocity_LocalSpace vem_local_space;

        const auto reference_element_data = vem_reference_element.Create(k);
        const auto local_space = vem_local_space.CreateLocalSpace(reference_element_data, polygon);

        // Test VEM performances
        Polydim::VEM::MCC::VEM_MCC_PerformanceAnalysis performanceAnalysis;

        const auto result =
            performanceAnalysis.Compute(Polydim::Utilities::Monomials_2D(), reference_element_data.MonomialsKp1, vem_local_space, local_space);

        ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(1.0e-10, result.ErrorPi0k, geometry_utilities.Tolerance1D()));

        ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(1.0e-10, result.ErrorGBD, geometry_utilities.Tolerance1D()));

        ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(1.0e-10, result.ErrorStabilization, geometry_utilities.Tolerance1D()));
    }
}

TEST(Test_VEM_MCC, Test_VEM_MCC_EdgeOrtho_2D_O0_O1_O2_O3)
{
    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = std::numeric_limits<double>::epsilon();
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    const auto polygon_data = Test_VEM_MCC_2D_Geometry(geometry_utilities);

    Polydim::VEM::MCC::VEM_MCC_2D_Polygon_Geometry polygon = {geometry_utilities_config.Tolerance1D,
                                                              geometry_utilities_config.Tolerance2D,
                                                              polygon_data.Vertices,
                                                              polygon_data.Centroid,
                                                              polygon_data.Measure,
                                                              polygon_data.Diameter,
                                                              polygon_data.TriangulationVertices,
                                                              polygon_data.EdgesLength,
                                                              polygon_data.EdgesDirection,
                                                              polygon_data.EdgesTangent,
                                                              polygon_data.EdgesNormal};

    for (unsigned int k = 0; k < 4; k++)
    {
        Polydim::VEM::MCC::VEM_MCC_2D_EdgeOrtho_Velocity_ReferenceElement vem_reference_element;
        Polydim::VEM::MCC::VEM_MCC_2D_EdgeOrtho_Velocity_LocalSpace vem_local_space;

        const auto reference_element_data = vem_reference_element.Create(k);
        const auto local_space = vem_local_space.CreateLocalSpace(reference_element_data, polygon);

        // Test VEM performances
        Polydim::VEM::MCC::VEM_MCC_PerformanceAnalysis performanceAnalysis;

        const auto result =
            performanceAnalysis.Compute(Polydim::Utilities::Monomials_2D(), reference_element_data.MonomialsKp1, vem_local_space, local_space);

        ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(1.0e-10, result.ErrorPi0k, geometry_utilities.Tolerance1D()));

        ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(1.0e-10, result.ErrorGBD, geometry_utilities.Tolerance1D()));

        ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(1.0e-10, result.ErrorStabilization, geometry_utilities.Tolerance1D()));
    }
}

TEST(Test_VEM_MCC, Test_VEM_MCC_Ortho_EdgeOrtho_2D_O0_O1_O2_O3)
{
    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = std::numeric_limits<double>::epsilon();
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    const auto polygon_data = Test_VEM_MCC_2D_Geometry(geometry_utilities);

    Polydim::VEM::MCC::VEM_MCC_2D_Polygon_Geometry polygon = {geometry_utilities_config.Tolerance1D,
                                                              geometry_utilities_config.Tolerance2D,
                                                              polygon_data.Vertices,
                                                              polygon_data.Centroid,
                                                              polygon_data.Measure,
                                                              polygon_data.Diameter,
                                                              polygon_data.TriangulationVertices,
                                                              polygon_data.EdgesLength,
                                                              polygon_data.EdgesDirection,
                                                              polygon_data.EdgesTangent,
                                                              polygon_data.EdgesNormal};

    for (unsigned int k = 0; k < 4; k++)
    {
        Polydim::VEM::MCC::VEM_MCC_2D_EdgeOrtho_Velocity_ReferenceElement vem_reference_element;
        Polydim::VEM::MCC::VEM_MCC_2D_Ortho_EdgeOrtho_Velocity_LocalSpace vem_local_space;

        const auto reference_element_data = vem_reference_element.Create(k);
        const auto local_space = vem_local_space.CreateLocalSpace(reference_element_data, polygon);

        // Test VEM performances
        Polydim::VEM::MCC::VEM_MCC_PerformanceAnalysis performanceAnalysis;

        const auto result =
            performanceAnalysis.Compute(Polydim::Utilities::Monomials_2D(), reference_element_data.MonomialsKp1, vem_local_space, local_space);

        ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(1.0e-10, result.ErrorPi0k, geometry_utilities.Tolerance1D()));

        ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(1.0e-10, result.ErrorGBD, geometry_utilities.Tolerance1D()));

        ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(1.0e-10, result.ErrorStabilization, geometry_utilities.Tolerance1D()));
    }
}

} // namespace UnitTesting
} // namespace Polydim

#endif
