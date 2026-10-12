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

#ifndef __TEST_VEM_DF_PCC_3D_LocalSpace_H
#define __TEST_VEM_DF_PCC_3D_LocalSpace_H

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <numbers>

#include "GeometryUtilities.hpp"
#include "VEM_DF_PCC_3D_ReferenceElement.hpp"
#include "VEM_DF_PCC_3D_Velocity_LocalSpace.hpp"
#include "VEM_DF_PCC_PerformanceAnalysis.hpp"
#include "VEM_PCC_2D_ReferenceElement.hpp"
#include "VTKUtilities.hpp"

namespace Polydim
{
namespace UnitTesting
{
struct Test_VEM_DF_PCC_3D_PolygonalFaces_Geometry final
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

struct Test_VEM_DF_PCC_3D_Polyhedron_Geometry final
{
    std::vector<Test_VEM_DF_PCC_3D_PolygonalFaces_Geometry> PolygonalFaces;

    Eigen::MatrixXd Vertices;
    Eigen::MatrixXi Edges;
    std::vector<Eigen::MatrixXi> Faces;
    Eigen::Vector3d Centroid;
    double Measure;
    double Diameter;
    std::vector<Eigen::MatrixXd> TetrahedronVertices;

    std::vector<double> FacesMeasure;
    std::vector<Eigen::Matrix3d> FacesRotationMatrix;
    std::vector<Eigen::Vector3d> FacesTranslation;
    std::vector<Eigen::Vector3d> FacesNormals;
    std::vector<bool> FacesNormalDirection;
    std::vector<bool> FacesNormalGlobalDirection;
    std::vector<std::array<Eigen::Vector3d, 2>> FacesTangents;
    std::vector<std::array<bool, 2>> FacesTangentsGlobalDirection;

    std::vector<bool> EdgesDirection;
    Eigen::MatrixXd EdgesTangent;
};

void Test_VEM_DF_PCC_3D_Velocity_Export_Dofs(const Polydim::VEM::DF_PCC::VEM_DF_PCC_3D_Polyhedron_Geometry &polyhedron,
                                             const std::vector<Polydim::VEM::PCC::VEM_PCC_2D_Polygon_Geometry> &polygonalFaces,
                                             const Polydim::VEM::DF_PCC::VEM_DF_PCC_3D_Velocity_ReferenceElement_Data &vem_reference_element,
                                             const std::string &exportVtuFolder)
{
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = polyhedron.Tolerance1D;
    geometryUtilitiesConfig.Tolerance2D = polyhedron.Tolerance2D;
    geometryUtilitiesConfig.Tolerance3D = polyhedron.Tolerance3D;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    const unsigned int num_vertices = polyhedron.Vertices.cols();
    const unsigned int num_edges = polyhedron.Edges.cols();
    const unsigned int num_faces = polyhedron.Faces.size();
    const unsigned int num_dofs = vem_reference_element.NumDofs0D * num_vertices +
                                  vem_reference_element.NumDofs1D * num_edges + num_faces * vem_reference_element.NumDofs2D +
                                  vem_reference_element.NumDofs3D_BigOPlus + vem_reference_element.NumDofs3D_Divergence;
    const unsigned int num_dofs_3d = vem_reference_element.NumDofs3D_BigOPlus + vem_reference_element.NumDofs3D_Divergence;

    Eigen::MatrixXd dofs_coordinate = Eigen::MatrixXd::Zero(3, num_dofs);
    std::vector<std::array<double, 3>> dof_global_index_values(num_dofs);
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
                dof_global_index_values[id_dofs][0] = id_dofs;
                dof_global_index_values[id_dofs][1] =
                    (vem_reference_element.NumDofs0D * num_vertices + vem_reference_element.NumDofs1D * num_edges) + id_dofs;
                dof_global_index_values[id_dofs][2] =
                    2 * (vem_reference_element.NumDofs0D * num_vertices + vem_reference_element.NumDofs1D * num_edges) + id_dofs;
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
            const Eigen::Vector3d edge_origin = polyhedron.Vertices.col(polyhedron.Edges(0, c));
            const Eigen::Vector3d edge_tangent = polyhedron.Vertices.col(polyhedron.Edges(1, c)) - edge_origin;

            for (unsigned int loc_i = 0; loc_i < vem_reference_element.NumDofs1D; ++loc_i)
            {
                dof_cell_index_values[id_dofs] = c;
                dof_dimension_values[id_dofs] = 1;
                dofs_coordinate.col(id_dofs) = edge_origin + local_edge_coordinates[loc_i] * edge_tangent;
                dof_global_index_values[id_dofs][0] = id_dofs;
                dof_global_index_values[id_dofs][1] =
                    (vem_reference_element.NumDofs0D * num_vertices + vem_reference_element.NumDofs1D * num_edges) + id_dofs;
                dof_global_index_values[id_dofs][2] =
                    2 * (vem_reference_element.NumDofs0D * num_vertices + vem_reference_element.NumDofs1D * num_edges) + id_dofs;
                id_dofs++;
            }
        }
    }

    unsigned int dof_offset = 3 * (vem_reference_element.NumDofs0D * num_vertices + vem_reference_element.NumDofs1D * num_edges);

    if (vem_reference_element.NumDofs2D != 0)
    {
        for (unsigned int c = 0; c < num_faces; c++)
        {
            const auto local_polygon_coordinates =
                geometryUtilities.EquispaceCoordinates(vem_reference_element.NumDofs2D + 1, 0.0, 1.0, true);
            const Eigen::Vector3d polygon_centroid = polygonalFaces[c].Centroid;
            const auto polygonCentroidEdgesDistance =
                geometryUtilities.PolygonCentroidEdgesDistance(polygonalFaces[c].Vertices,
                                                               polygonalFaces[c].Centroid,
                                                               polygonalFaces[c].EdgesNormal);

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
                dof_global_index_values[id_dofs][0] = dof_offset + id_dofs;
                dof_global_index_values[id_dofs][1] = dof_offset + vem_reference_element.NumDofs2D + id_dofs;
                dof_global_index_values[id_dofs][2] = dof_offset + vem_reference_element.NumDofs2D + id_dofs;
                id_dofs++;
            }
        }
    }

    dof_offset += 3 * vem_reference_element.NumDofs2D;

    if (num_dofs_3d != 0)
    {
        const auto faces3DVertices = geometryUtilities.PolyhedronFaceVertices(polyhedron.Vertices, polyhedron.Faces);

        const auto local_polyhedron_coordinates = geometryUtilities.fibonacci_sphere(num_dofs_3d);
        const Eigen::Vector3d polyhedron_centroid = polyhedron.Centroid;
        const auto polyhedron_centroid_faces_distance =
            geometryUtilities.PolyhedronCentroidFacesDistance(polyhedron_centroid, polyhedron.FacesNormal, faces3DVertices);
        const double polyhedron_in_radius = geometryUtilities.PolyhedronInRadius(polyhedron_centroid_faces_distance);

        double sphere_diameter = 0.0;
        if (num_dofs_3d > 0)
            sphere_diameter = 0.5 * polyhedron_in_radius;

        for (unsigned int loc_i = 0; loc_i < num_dofs_3d; ++loc_i)
        {
            dof_cell_index_values[id_dofs] = 0;
            dof_dimension_values[id_dofs] = 3;
            dofs_coordinate.col(id_dofs) = polyhedron_centroid + sphere_diameter * local_polyhedron_coordinates.col(loc_i);
            dof_global_index_values[id_dofs][0] = dof_offset + id_dofs;
            dof_global_index_values[id_dofs][1] = std::nan("");
            dof_global_index_values[id_dofs][2] = std::nan("");
            id_dofs++;
        }
    }

    Eigen::VectorXd dof_global_index_values_data(3 * num_dofs);
    unsigned int c = 0;
    for (const auto &v : dof_global_index_values)
    {
        for (unsigned int d = 0; d < 3; d++)
            dof_global_index_values_data[3 * c + d] = v[d];
        c++;
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
                             static_cast<unsigned int>(dof_global_index_values_data.size()),
                             dof_global_index_values_data.data()}});

        exporter.Export(exportVtuFolder + "/velocity_dofs_" + std::to_string(vem_reference_element.Order) + ".vtu");
    }
}

void Test_VEM_DF_PCC_3D_Pressure_Export_Dofs(const Polydim::VEM::DF_PCC::VEM_DF_PCC_3D_Polyhedron_Geometry &polyhedron,
                                             const std::vector<Polydim::VEM::PCC::VEM_PCC_2D_Polygon_Geometry> &polygonalFaces,
                                             const Polydim::VEM::DF_PCC::VEM_DF_PCC_3D_Pressure_ReferenceElement_Data &vem_reference_element,
                                             const std::string &exportVtuFolder)
{
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = polyhedron.Tolerance1D;
    geometryUtilitiesConfig.Tolerance2D = polyhedron.Tolerance2D;
    geometryUtilitiesConfig.Tolerance3D = polyhedron.Tolerance3D;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    const unsigned int num_vertices = polyhedron.Vertices.cols();
    const unsigned int num_edges = polyhedron.Edges.cols();
    const unsigned int num_faces = polyhedron.Faces.size();
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
            const Eigen::Vector3d edge_origin = polyhedron.Vertices.col(polyhedron.Edges(0, c));
            const Eigen::Vector3d edge_tangent = polyhedron.Vertices.col(polyhedron.Edges(1, c)) - edge_origin;

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
        for (unsigned int c = 0; c < num_faces; c++)
        {
            const auto local_polygon_coordinates =
                geometryUtilities.EquispaceCoordinates(vem_reference_element.NumDofs2D + 1, 0.0, 1.0, true);
            const Eigen::Vector3d polygon_centroid = polygonalFaces[c].Centroid;
            const auto polygonCentroidEdgesDistance =
                geometryUtilities.PolygonCentroidEdgesDistance(polygonalFaces[c].Vertices,
                                                               polygonalFaces[c].Centroid,
                                                               polygonalFaces[c].EdgesNormal);

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
        const auto faces3DVertices = geometryUtilities.PolyhedronFaceVertices(polyhedron.Vertices, polyhedron.Faces);

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

Test_VEM_DF_PCC_3D_Polyhedron_Geometry Test_VEM_DF_PCC_3D_Geometry(const Gedim::GeometryUtilities &geometry_utilities)
{
    Test_VEM_DF_PCC_3D_Polyhedron_Geometry result;

    const Eigen::Vector3d origin = Eigen::Vector3d(0.0, 0.0, 0.0);
    const Eigen::Vector3d length = Eigen::Vector3d(1.0, 0.0, 0.0);
    const Eigen::Vector3d height = Eigen::Vector3d(0.0, 0.0, 1.0);
    const Eigen::Vector3d width = Eigen::Vector3d(0.0, 1.0, 0.0);
    Gedim::GeometryUtilities::Polyhedron cube = geometry_utilities.CreateParallelepipedWithOrigin(origin, length, height, width);

    result.Vertices = cube.Vertices;
    result.Edges = cube.Edges;
    result.Faces = cube.Faces;

    result.Measure = 1.0e+00;
    result.Diameter = sqrt(3.0);
    result.Centroid << 0.5, 0.5, 0.5;
    result.EdgesTangent = geometry_utilities.PolyhedronEdgeTangents(result.Vertices, result.Edges);
    result.EdgesDirection.resize(12, true);

    const std::vector<Eigen::MatrixXd> facesVertices = geometry_utilities.PolyhedronFaceVertices(result.Vertices, result.Faces);
    result.FacesTranslation = geometry_utilities.PolyhedronFaceTranslations(facesVertices);
    result.FacesNormals = geometry_utilities.PolyhedronFaceNormals(facesVertices);
    result.FacesRotationMatrix =
        geometry_utilities.PolyhedronFaceRotationMatrices(facesVertices, result.FacesNormals, result.FacesTranslation);
    result.FacesNormalDirection =
        geometry_utilities.PolyhedronFaceNormalDirections(facesVertices, result.Centroid, result.FacesNormals);

    const unsigned int numFaces = result.Faces.size();
    result.PolygonalFaces.resize(numFaces);

    std::vector<std::vector<bool>> facesEdgeDirections =
        geometry_utilities.PolyhedronFaceEdgeDirections(result.Vertices, result.Edges, result.Faces);

    for (unsigned int f = 0; f < numFaces; f++)
    {
        result.PolygonalFaces[f].Vertices = geometry_utilities.RotatePointsFrom3DTo2D(facesVertices[f],
                                                                                      result.FacesRotationMatrix[f].transpose(),
                                                                                      result.FacesTranslation[f]);
        std::vector<unsigned int> facesTriangulation =
            geometry_utilities.PolygonTriangulationByFirstVertex(result.PolygonalFaces[f].Vertices);
        result.PolygonalFaces[f].TriangulationVertices =
            geometry_utilities.ExtractTriangulationPoints(result.PolygonalFaces[f].Vertices, facesTriangulation);
        result.PolygonalFaces[f].Measure = geometry_utilities.PolygonArea(result.PolygonalFaces[f].Vertices);
        result.PolygonalFaces[f].Diameter = geometry_utilities.PolygonDiameter(result.PolygonalFaces[f].Vertices);
        result.PolygonalFaces[f].Centroid =
            geometry_utilities.PolygonCentroid(result.PolygonalFaces[f].Vertices, result.PolygonalFaces[f].Measure);

        result.PolygonalFaces[f].EdgesLength = geometry_utilities.PolygonEdgeLengths(result.PolygonalFaces[f].Vertices);
        result.PolygonalFaces[f].EdgesNormal = geometry_utilities.PolygonEdgeNormals(result.PolygonalFaces[f].Vertices);
        result.PolygonalFaces[f].EdgesTangent = geometry_utilities.PolygonEdgeTangents(result.PolygonalFaces[f].Vertices);
        result.PolygonalFaces[f].EdgesDirection = facesEdgeDirections[f];
    }

    std::vector<std::vector<unsigned int>> facesTriangulation3D =
        geometry_utilities.PolyhedronFaceTriangulationsByFirstVertex(result.Faces, facesVertices);

    std::vector<unsigned int> polyhedronTetrahedrons =
        geometry_utilities.PolyhedronTetrahedronsByFaceTriangulations(result.Vertices,
                                                                      result.Faces,
                                                                      facesTriangulation3D,
                                                                      result.Centroid);

    result.TetrahedronVertices =
        geometry_utilities.ExtractTetrahedronPoints(result.Vertices, result.Centroid, polyhedronTetrahedrons);

    result.FacesTangents = geometry_utilities.PolyhedronFaceTangents(facesVertices, result.FacesNormals, result.FacesNormalDirection);

    result.FacesNormalGlobalDirection.resize(numFaces, true);
    result.FacesTangentsGlobalDirection.resize(numFaces);

    for (unsigned int f = 0; f < numFaces; ++f)
    {
        result.FacesTangentsGlobalDirection[f][0] = facesEdgeDirections[f][0];
        result.FacesTangentsGlobalDirection[f][1] =
            result.FacesTangentsGlobalDirection[f][0] == result.FacesNormalGlobalDirection[f];
    }

    return result;
}

TEST(Test_VEM_DF_PCC, Test_VEM_DF_PCC_3D_O2_O3_O4)
{
    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = std::numeric_limits<double>::epsilon();
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    const Test_VEM_DF_PCC_3D_Polyhedron_Geometry polyhedron_data = Test_VEM_DF_PCC_3D_Geometry(geometry_utilities);

    const std::string exportFolder = "Export/VEM/DF_PCC/Test_VEM_DF_PCC_3D_O2_O3_O4";
    Gedim::Output::CreateFolder(exportFolder);

    const unsigned int numFaces = polyhedron_data.PolygonalFaces.size();
    std::vector<Polydim::VEM::PCC::VEM_PCC_2D_Polygon_Geometry> polygonalFaces;

    for (unsigned int f = 0; f < numFaces; f++)
    {
        polygonalFaces.push_back({geometry_utilities_config.Tolerance1D,
                                  geometry_utilities_config.Tolerance2D,
                                  polyhedron_data.PolygonalFaces[f].Vertices,
                                  polyhedron_data.PolygonalFaces[f].Centroid,
                                  polyhedron_data.PolygonalFaces[f].Measure,
                                  polyhedron_data.PolygonalFaces[f].Diameter,
                                  polyhedron_data.PolygonalFaces[f].TriangulationVertices,
                                  polyhedron_data.PolygonalFaces[f].EdgesLength,
                                  polyhedron_data.PolygonalFaces[f].EdgesDirection,
                                  polyhedron_data.PolygonalFaces[f].EdgesTangent,
                                  polyhedron_data.PolygonalFaces[f].EdgesNormal});
    }

    Polydim::VEM::DF_PCC::VEM_DF_PCC_3D_Polyhedron_Geometry polyhedron = {
        geometry_utilities_config.Tolerance1D,
        geometry_utilities_config.Tolerance2D,
        geometry_utilities_config.Tolerance3D,

        polyhedron_data.Vertices,

        polyhedron_data.Edges,
        polyhedron_data.Faces,
        polyhedron_data.Centroid,
        polyhedron_data.Measure,
        polyhedron_data.Diameter,
        polyhedron_data.TetrahedronVertices,

        polyhedron_data.FacesRotationMatrix,
        polyhedron_data.FacesTranslation,
        polyhedron_data.FacesNormals,
        polyhedron_data.FacesNormalDirection,
        polyhedron_data.FacesNormalGlobalDirection,
        polyhedron_data.FacesTangents,
        polyhedron_data.FacesTangentsGlobalDirection,

        polyhedron_data.EdgesDirection,
        polyhedron_data.EdgesTangent,
    };

    for (unsigned int k = 2; k < 5; k++)
    {
        Polydim::VEM::PCC::VEM_PCC_2D_ReferenceElement vem_reference_element_2D;
        Polydim::VEM::DF_PCC::VEM_DF_PCC_3D_Velocity_ReferenceElement vem_reference_element_3D;
        Polydim::VEM::DF_PCC::VEM_DF_PCC_3D_Pressure_ReferenceElement vem_pressure_reference_element_3D;
        Polydim::VEM::DF_PCC::VEM_DF_PCC_3D_Velocity_LocalSpace vem_local_space;

        const auto reference_element_data_2D = vem_reference_element_2D.Create(k);
        const auto reference_element_data_3D = vem_reference_element_3D.Create(k);
        const auto pressure_reference_element_data_3D = vem_pressure_reference_element_3D.Create(k);
        const auto local_space =
            vem_local_space.CreateLocalSpace(reference_element_data_2D, reference_element_data_3D, polygonalFaces, polyhedron);

        // Export domain
        {
            Gedim::VTKUtilities vtkUtilities;
            vtkUtilities.AddPolyhedron(polyhedron.Vertices, polyhedron.Edges, polyhedron.Faces);
            vtkUtilities.Export(exportFolder + "/Polyhedron.vtu");

            Test_VEM_DF_PCC_3D_Velocity_Export_Dofs(polyhedron, polygonalFaces, reference_element_data_3D, exportFolder);
            Test_VEM_DF_PCC_3D_Pressure_Export_Dofs(polyhedron, polygonalFaces, pressure_reference_element_data_3D, exportFolder);
        }

        // Test VEM performances
        Polydim::VEM::DF_PCC::VEM_DF_PCC_PerformanceAnalysis performanceAnalysis;

        const auto result =
            performanceAnalysis.Compute(Polydim::Utilities::Monomials_3D(), reference_element_data_3D.Monomials, vem_local_space, local_space);

        for (unsigned int d1 = 0; d1 < local_space.Dimension; d1++)
        {
            ASSERT_TRUE(result.ErrorPiNabla[d1] < 1.0e-12);
            ASSERT_TRUE(result.ErrorPi0km2[d1] < 1.0e-12);
            ASSERT_TRUE(result.ErrorPi0k[d1] < 1.0e-11);
            ASSERT_EQ(result.ErrorPi0km1Grad.size(), local_space.Dimension * local_space.Dimension);
            for (unsigned int d2 = 0; d2 < local_space.Dimension; ++d2)
                ASSERT_TRUE(result.ErrorPi0km1Grad[d1 * local_space.Dimension + d2] < 1.0e-12);
        }

        for (unsigned int d1 = 0; d1 < local_space.Dimension; d1++)
        {
            ASSERT_TRUE(result.ErrorGBD[d1] < 1.0e-12);
            ASSERT_TRUE(result.ErrorHCD[d1] < 1.0e-12);
            ASSERT_EQ(result.ErrorHED.size(), local_space.Dimension * local_space.Dimension);
            for (unsigned int d2 = 0; d2 < local_space.Dimension; ++d2)
                ASSERT_TRUE(result.ErrorHED[d1 * local_space.Dimension + d2] < 1.0e-12);
        }

        ASSERT_TRUE(result.ErrorStabilization < 1.0e-8);
    }
}

} // namespace UnitTesting
} // namespace Polydim

#endif
