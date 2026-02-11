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

#ifndef __TEST_VEM_MCC_3D_LocalSpace_H
#define __TEST_VEM_MCC_3D_LocalSpace_H

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <numbers>

#include "GeometryUtilities.hpp"
#include "VEM_MCC_3D_ReferenceElement.hpp"
#include "VEM_MCC_3D_Velocity_LocalSpace.hpp"
#include "VEM_MCC_PerformanceAnalysis.hpp"
#include "VTKUtilities.hpp"

namespace Polydim
{
namespace UnitTesting
{

struct Test_VEM_MCC_3D_Polyhedron_Geometry final
{
    Eigen::MatrixXd Vertices;
    Eigen::Vector3d Centroid;
    double Measure;
    double Diameter;
    std::vector<Eigen::MatrixXd> TetrahedronVertices;

    std::vector<Eigen::Matrix3d> FacesRotationMatrix;
    std::vector<Eigen::Vector3d> FacesTranslation;
    std::vector<Eigen::Vector3d> FacesNormal;
    std::vector<bool> FacesNormalDirection;
    std::vector<bool> FacesNormalGlobalDirection;

    std::vector<double> FacesMeasure;
    std::vector<Eigen::Vector3d> FacesCentroid2D;
    std::vector<double> FacesDiameter;
    std::vector<std::vector<Eigen::Matrix3d>> FacesTriangulationVertices2D;
};

void Test_VEM_MCC_3D_Velocity_Export_Dofs(const Gedim::GeometryUtilities::Polyhedron &cube,
                                          const Polydim::VEM::MCC::VEM_MCC_3D_Polyhedron_Geometry &polyhedron,
                                          const Polydim::VEM::MCC::VEM_MCC_3D_Velocity_ReferenceElement_Data &vem_reference_element,
                                          const std::string &exportVtuFolder)
{
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = polyhedron.Tolerance1D;
    geometryUtilitiesConfig.Tolerance2D = polyhedron.Tolerance2D;
    geometryUtilitiesConfig.Tolerance3D = polyhedron.Tolerance3D;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    const unsigned int num_vertices = polyhedron.Vertices.cols();
    const unsigned int num_edges = cube.Edges.cols();
    const unsigned int num_faces = cube.Faces.size();
    const unsigned int num_dofs = vem_reference_element.NumDofs0D * num_vertices + vem_reference_element.NumDofs1D * num_edges +
                                  num_faces * vem_reference_element.NumDofs2D + vem_reference_element.NumDofs3D;

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
                dofs_coordinate.col(id_dofs) = polyhedron.Vertices.col(c);
                dof_global_index_values[id_dofs] = id_dofs;
                id_dofs++;
            }
        }
    }

    if (vem_reference_element.NumDofs1D != 0)
    {
        for (unsigned int c = 0; c < num_edges; ++c)
        {
            const std::vector<double> local_edge_coordinates =
                geometryUtilities.EquispaceCoordinates(vem_reference_element.NumDofs1D, 0.0, 1.0, false);
            const Eigen::Vector3d edge_origin = polyhedron.Vertices.col(cube.Edges(0, c));
            const Eigen::Vector3d edge_tangent = polyhedron.Vertices.col(cube.Edges(1, c)) - edge_origin;

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
        const std::vector<Eigen::MatrixXd> faces_vertices_3d = geometryUtilities.PolyhedronFaceVertices(cube.Vertices, cube.Faces);
        const std::vector<Eigen::MatrixXd> faces_vertices_2d =
            geometryUtilities.PolyhedronFaceRotatedVertices(faces_vertices_3d, polyhedron.FacesTranslation, polyhedron.FacesRotationMatrix);
        for (unsigned int c = 0; c < num_faces; c++)
        {
            const auto local_polygon_coordinates =
                geometryUtilities.EquispaceCoordinates(vem_reference_element.NumDofs2D + 1, 0.0, 1.0, true);
            const Eigen::Vector3d polygon_centroid = polyhedron.FacesCentroid2D[c];
            const auto faces_edge_normal = geometryUtilities.PolygonEdgeNormals(faces_vertices_2d[c]);
            const auto polygonCentroidEdgesDistance =
                geometryUtilities.PolygonCentroidEdgesDistance(faces_vertices_2d[c], polygon_centroid, faces_edge_normal);

            double circle_diameter = 0.0;
            if (vem_reference_element.NumDofs2D > 1)
                circle_diameter = 0.5 * geometryUtilities.PolygonInRadius(polygonCentroidEdgesDistance);

            for (unsigned int loc_i = 0; loc_i < vem_reference_element.NumDofs2D; ++loc_i)
            {
                dof_cell_index_values[id_dofs] = 0;
                dof_dimension_values[id_dofs] = 2;
                dofs_coordinate.col(id_dofs) = geometryUtilities.RotatePointsFrom2DTo3D(
                    polygon_centroid +
                        circle_diameter * Eigen::Vector3d(cos(2.0 * std::numbers::pi * local_polygon_coordinates.at(loc_i)),
                                                          sin(2.0 * std::numbers::pi * local_polygon_coordinates.at(loc_i)),
                                                          0.0),
                    polyhedron.FacesRotationMatrix[c],
                    polyhedron.FacesTranslation[c]);
                dof_global_index_values[id_dofs] = id_dofs;
                id_dofs++;
            }
        }
    }

    if (vem_reference_element.NumDofs3D != 0)
    {
        const auto faces3DVertices = geometryUtilities.PolyhedronFaceVertices(polyhedron.Vertices, cube.Faces);

        const auto local_polyhedron_coordinates = geometryUtilities.fibonacci_sphere(vem_reference_element.NumDofs3D);
        const Eigen::Vector3d polyhedron_centroid = polyhedron.Centroid;
        const auto polyhedron_centroid_faces_distance =
            geometryUtilities.PolyhedronCentroidFacesDistance(polyhedron_centroid, polyhedron.FacesNormal, faces3DVertices);
        const double polyhedron_in_radius = geometryUtilities.PolyhedronInRadius(polyhedron_centroid_faces_distance);

        double sphere_diameter = 0.0;
        if (vem_reference_element.NumDofs3D > 0)
            sphere_diameter = 0.5 * polyhedron_in_radius;

        for (unsigned int loc_i = 0; loc_i < vem_reference_element.NumDofs3D; ++loc_i)
        {
            dof_cell_index_values[id_dofs] = 0;
            dof_dimension_values[id_dofs] = 3;
            dofs_coordinate.col(id_dofs) = polyhedron_centroid + sphere_diameter * local_polyhedron_coordinates.col(loc_i);
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

void Test_VEM_MCC_3D_Pressure_Export_Dofs(const Gedim::GeometryUtilities::Polyhedron &cube,
                                          const Polydim::VEM::MCC::VEM_MCC_3D_Polyhedron_Geometry &polyhedron,
                                          const Polydim::VEM::MCC::VEM_MCC_3D_Pressure_ReferenceElement_Data &vem_reference_element,
                                          const std::string &exportVtuFolder)
{
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = polyhedron.Tolerance1D;
    geometryUtilitiesConfig.Tolerance2D = polyhedron.Tolerance2D;
    geometryUtilitiesConfig.Tolerance3D = polyhedron.Tolerance3D;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    const unsigned int num_vertices = polyhedron.Vertices.cols();
    const unsigned int num_edges = cube.Edges.cols();
    const unsigned int num_faces = cube.Faces.size();
    const unsigned int num_dofs = vem_reference_element.NumDofs0D * num_vertices + vem_reference_element.NumDofs1D * num_edges +
                                  num_faces * vem_reference_element.NumDofs2D + vem_reference_element.NumDofs3D;

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
                dofs_coordinate.col(id_dofs) = polyhedron.Vertices.col(c);
                dof_global_index_values[id_dofs] = id_dofs;
                id_dofs++;
            }
        }
    }

    if (vem_reference_element.NumDofs1D != 0)
    {
        for (unsigned int c = 0; c < num_edges; ++c)
        {
            const std::vector<double> local_edge_coordinates =
                geometryUtilities.EquispaceCoordinates(vem_reference_element.NumDofs1D, 0.0, 1.0, false);
            const Eigen::Vector3d edge_origin = polyhedron.Vertices.col(cube.Edges(0, c));
            const Eigen::Vector3d edge_tangent = polyhedron.Vertices.col(cube.Edges(1, c)) - edge_origin;

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
        const std::vector<Eigen::MatrixXd> faces_vertices_3d = geometryUtilities.PolyhedronFaceVertices(cube.Vertices, cube.Faces);
        const std::vector<Eigen::MatrixXd> faces_vertices_2d =
            geometryUtilities.PolyhedronFaceRotatedVertices(faces_vertices_3d, polyhedron.FacesTranslation, polyhedron.FacesRotationMatrix);
        for (unsigned int c = 0; c < num_faces; c++)
        {
            const auto local_polygon_coordinates =
                geometryUtilities.EquispaceCoordinates(vem_reference_element.NumDofs2D + 1, 0.0, 1.0, true);
            const Eigen::Vector3d polygon_centroid = polyhedron.FacesCentroid2D[c];
            const auto faces_edge_normal = geometryUtilities.PolygonEdgeNormals(faces_vertices_2d[c]);
            const auto polygonCentroidEdgesDistance =
                geometryUtilities.PolygonCentroidEdgesDistance(faces_vertices_2d[c], polygon_centroid, faces_edge_normal);

            double circle_diameter = 0.0;
            if (vem_reference_element.NumDofs2D > 1)
                circle_diameter = 0.2 * geometryUtilities.PolygonInRadius(polygonCentroidEdgesDistance);

            for (unsigned int loc_i = 0; loc_i < vem_reference_element.NumDofs2D; ++loc_i)
            {
                dof_cell_index_values[id_dofs] = 0;
                dof_dimension_values[id_dofs] = 2;
                dofs_coordinate.col(id_dofs) = geometryUtilities.RotatePointsFrom2DTo3D(
                    polygon_centroid +
                        circle_diameter * Eigen::Vector3d(cos(2.0 * std::numbers::pi * local_polygon_coordinates.at(loc_i)),
                                                          sin(2.0 * std::numbers::pi * local_polygon_coordinates.at(loc_i)),
                                                          0.0),
                    polyhedron.FacesRotationMatrix[c],
                    polyhedron.FacesTranslation[c]);
                dof_global_index_values[id_dofs] = id_dofs;
                id_dofs++;
            }
        }
    }

    if (vem_reference_element.NumDofs3D != 0)
    {
        const auto faces3DVertices = geometryUtilities.PolyhedronFaceVertices(polyhedron.Vertices, cube.Faces);

        const auto local_polyhedron_coordinates = geometryUtilities.fibonacci_sphere(vem_reference_element.NumDofs3D);
        const Eigen::Vector3d polyhedron_centroid = polyhedron.Centroid;
        const auto polyhedron_centroid_faces_distance =
            geometryUtilities.PolyhedronCentroidFacesDistance(polyhedron_centroid, polyhedron.FacesNormal, faces3DVertices);
        const double polyhedron_in_radius = geometryUtilities.PolyhedronInRadius(polyhedron_centroid_faces_distance);

        double sphere_diameter = 0.0;
        if (vem_reference_element.NumDofs3D > 0)
            sphere_diameter = 0.2 * polyhedron_in_radius;

        for (unsigned int loc_i = 0; loc_i < vem_reference_element.NumDofs3D; ++loc_i)
        {
            dof_cell_index_values[id_dofs] = 0;
            dof_dimension_values[id_dofs] = 3;
            dofs_coordinate.col(id_dofs) = polyhedron_centroid + sphere_diameter * local_polyhedron_coordinates.col(loc_i);
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

Test_VEM_MCC_3D_Polyhedron_Geometry Test_VEM_MCC_3D_Geometry(const Gedim::GeometryUtilities &geometry_utilities,
                                                             Gedim::GeometryUtilities::Polyhedron &cube)
{
    Test_VEM_MCC_3D_Polyhedron_Geometry result;

    const Eigen::Vector3d origin = Eigen::Vector3d(0.0, 0.0, 0.0);
    const Eigen::Vector3d length = Eigen::Vector3d(1.0, 0.0, 0.0);
    const Eigen::Vector3d height = Eigen::Vector3d(0.0, 0.0, 1.0);
    const Eigen::Vector3d width = Eigen::Vector3d(0.0, 1.0, 0.0);
    cube = geometry_utilities.CreateParallelepipedWithOrigin(origin, length, height, width);

    result.Vertices = cube.Vertices;

    result.Measure = 1.0e+00;
    result.Diameter = sqrt(3.0);
    result.Centroid << 0.5, 0.5, 0.5;

    const std::vector<Eigen::MatrixXd> facesVertices = geometry_utilities.PolyhedronFaceVertices(result.Vertices, cube.Faces);
    result.FacesTranslation = geometry_utilities.PolyhedronFaceTranslations(facesVertices);
    result.FacesNormal = geometry_utilities.PolyhedronFaceNormals(facesVertices);
    result.FacesRotationMatrix =
        geometry_utilities.PolyhedronFaceRotationMatrices(facesVertices, result.FacesNormal, result.FacesTranslation);
    result.FacesNormalDirection =
        geometry_utilities.PolyhedronFaceNormalDirections(facesVertices, result.Centroid, result.FacesNormal);

    const unsigned int numFaces = cube.Faces.size();
    result.FacesTriangulationVertices2D.resize(numFaces);
    result.FacesMeasure.resize(numFaces);
    result.FacesDiameter.resize(numFaces);
    result.FacesCentroid2D.resize(numFaces);
    result.FacesNormalGlobalDirection.resize(numFaces, true);

    for (unsigned int f = 0; f < numFaces; f++)
    {
        const Eigen::MatrixXd facesVertices2D =
            geometry_utilities.RotatePointsFrom3DTo2D(facesVertices[f],
                                                      result.FacesRotationMatrix[f].transpose(),
                                                      result.FacesTranslation[f]);
        std::vector<unsigned int> facesTriangulation = geometry_utilities.PolygonTriangulationByFirstVertex(facesVertices2D);
        result.FacesTriangulationVertices2D[f] = geometry_utilities.ExtractTriangulationPoints(facesVertices2D, facesTriangulation);
        result.FacesMeasure[f] = geometry_utilities.PolygonArea(facesVertices2D);
        result.FacesDiameter[f] = geometry_utilities.PolygonDiameter(facesVertices2D);
        result.FacesCentroid2D[f] = geometry_utilities.PolygonCentroid(facesVertices2D, result.FacesMeasure[f]);
    }

    std::vector<std::vector<unsigned int>> facesTriangulation3D =
        geometry_utilities.PolyhedronFaceTriangulationsByFirstVertex(cube.Faces, facesVertices);

    std::vector<unsigned int> polyhedronTetrahedrons =
        geometry_utilities.PolyhedronTetrahedronsByFaceTriangulations(result.Vertices,
                                                                      cube.Faces,
                                                                      facesTriangulation3D,
                                                                      result.Centroid);

    result.TetrahedronVertices =
        geometry_utilities.ExtractTetrahedronPoints(result.Vertices, result.Centroid, polyhedronTetrahedrons);

    return result;
}

TEST(Test_VEM_MCC, Test_VEM_MCC_3D_O0_O1_O2_O3)
{
    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = std::numeric_limits<double>::epsilon();
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    const std::string exportFolder = "Export/VEM/MCC/Test_VEM_MCC_3D_O0_O1_O2_O3";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilities::Polyhedron cube;
    const Test_VEM_MCC_3D_Polyhedron_Geometry polyhedron_data = Test_VEM_MCC_3D_Geometry(geometry_utilities, cube);

    Polydim::VEM::MCC::VEM_MCC_3D_Polyhedron_Geometry polyhedron = {geometry_utilities_config.Tolerance1D,
                                                                    geometry_utilities_config.Tolerance2D,
                                                                    geometry_utilities_config.Tolerance3D,

                                                                    polyhedron_data.Vertices,
                                                                    polyhedron_data.Centroid,
                                                                    polyhedron_data.Measure,
                                                                    polyhedron_data.Diameter,
                                                                    polyhedron_data.TetrahedronVertices,

                                                                    polyhedron_data.FacesRotationMatrix,
                                                                    polyhedron_data.FacesTranslation,
                                                                    polyhedron_data.FacesNormal,
                                                                    polyhedron_data.FacesNormalDirection,
                                                                    polyhedron_data.FacesNormalGlobalDirection,

                                                                    polyhedron_data.FacesMeasure,
                                                                    polyhedron_data.FacesCentroid2D,
                                                                    polyhedron_data.FacesDiameter,
                                                                    polyhedron_data.FacesTriangulationVertices2D};

    for (unsigned int k = 0; k < 4; k++)
    {
        Polydim::VEM::MCC::VEM_MCC_3D_Velocity_ReferenceElement vem_reference_element;
        Polydim::VEM::MCC::VEM_MCC_3D_Pressure_ReferenceElement vem_pressure_reference_element;
        Polydim::VEM::MCC::VEM_MCC_3D_Velocity_LocalSpace vem_local_space;

        const auto reference_element_data = vem_reference_element.Create(k);
        const auto pressure_reference_element_data = vem_pressure_reference_element.Create(k);
        const auto local_space = vem_local_space.CreateLocalSpace(reference_element_data, polyhedron);

        // Export domain
        {
            Gedim::VTKUtilities vtkUtilities;
            vtkUtilities.AddPolyhedron(cube.Vertices, cube.Edges, cube.Faces);
            vtkUtilities.Export(exportFolder + "/Polyhedron.vtu");

            Test_VEM_MCC_3D_Velocity_Export_Dofs(cube, polyhedron, reference_element_data, exportFolder);
            Test_VEM_MCC_3D_Pressure_Export_Dofs(cube, polyhedron, pressure_reference_element_data, exportFolder);
        }

        // Test VEM performances
        Polydim::VEM::MCC::VEM_MCC_PerformanceAnalysis performanceAnalysis;

        const auto result =
            performanceAnalysis.Compute(Polydim::Utilities::Monomials_3D(), reference_element_data.MonomialsKp1, vem_local_space, local_space);

        ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(1.0e-10, result.ErrorPi0k, geometry_utilities.Tolerance1D()));

        ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(1.0e-10, result.ErrorGBD, geometry_utilities.Tolerance1D()));

        ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(1.0e-10, result.ErrorStabilization, geometry_utilities.Tolerance1D()));
    }
}

} // namespace UnitTesting
} // namespace Polydim

#endif
