#ifndef __TEST_VEM_PCC_3D_Inertia_LocalSpace_H
#define __TEST_VEM_PCC_3D_Inertia_LocalSpace_H

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "GeometryUtilities.hpp"
#include "VEM_PCC_2D_LocalSpace_Data.hpp"
#include "VEM_PCC_3D_Inertia_LocalSpace.hpp"
#include "VEM_PCC_PerformanceAnalysis.hpp"

namespace Polydim
{
namespace UnitTesting
{
struct Test_VEM_PCC_2D_Inertia_PolygonalFaces_Geometry final
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

struct Test_VEM_PCC_3D_Inertia_Polyhedron_Geometry final
{
    std::vector<Test_VEM_PCC_2D_Inertia_PolygonalFaces_Geometry> PolygonalFaces;

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

    std::vector<bool> EdgesDirection;
    Eigen::MatrixXd EdgesTangent;
};

Test_VEM_PCC_3D_Inertia_Polyhedron_Geometry Test_VEM_PCC_3D_Inertia_Geometry(
    const Gedim::GeometryUtilities &geometry_utilities)
{
    Test_VEM_PCC_3D_Inertia_Polyhedron_Geometry result;

    const Eigen::Vector3d origin = Eigen::Vector3d(0.0, 0.0, 0.0);
    const Eigen::Vector3d length = Eigen::Vector3d(1.0, 0.0, 0.0);
    const Eigen::Vector3d height = Eigen::Vector3d(0.0, 0.0, 1.0);
    const Eigen::Vector3d width = Eigen::Vector3d(0.0, 1.0, 0.0);
    Gedim::GeometryUtilities::Polyhedron cube =
        geometry_utilities.CreateParallelepipedWithOrigin(origin, length, height, width);

    result.Vertices = cube.Vertices;
    result.Edges = cube.Edges;
    result.Faces = cube.Faces;

    result.Measure = 1.0e+00;
    result.Diameter = sqrt(3.0);
    result.Centroid << 0.5, 0.5, 0.5;
    result.EdgesTangent = geometry_utilities.PolyhedronEdgeTangents(result.Vertices, result.Edges);
    result.EdgesDirection.resize(12, true);

    const std::vector<Eigen::MatrixXd> facesVertices =
        geometry_utilities.PolyhedronFaceVertices(result.Vertices, result.Faces);
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
        result.PolygonalFaces[f].Vertices = geometry_utilities.RotatePointsFrom3DTo2D(
            facesVertices[f], result.FacesRotationMatrix[f].transpose(), result.FacesTranslation[f]);
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
        result.PolygonalFaces[f].EdgesTangent =
            geometry_utilities.PolygonEdgeTangents(result.PolygonalFaces[f].Vertices);
        result.PolygonalFaces[f].EdgesDirection = facesEdgeDirections[f];
    }

    std::vector<std::vector<unsigned int>> facesTriangulation3D =
        geometry_utilities.PolyhedronFaceTriangulationsByFirstVertex(result.Faces, facesVertices);

    std::vector<unsigned int> polyhedronTetrahedrons = geometry_utilities.PolyhedronTetrahedronsByFaceTriangulations(
        result.Vertices, result.Faces, facesTriangulation3D, result.Centroid);

    result.TetrahedronVertices =
        geometry_utilities.ExtractTetrahedronPoints(result.Vertices, result.Centroid, polyhedronTetrahedrons);

    return result;
}

TEST(Test_VEM_PCC, Test_VEM_PCC_3D_Inertia_O1_O2_O3)
{
    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = std::numeric_limits<double>::epsilon();
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    const Test_VEM_PCC_3D_Inertia_Polyhedron_Geometry polyhedron_data =
        Test_VEM_PCC_3D_Inertia_Geometry(geometry_utilities);

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

    Polydim::VEM::PCC::VEM_PCC_3D_Polyhedron_Geometry polyhedron = {
        geometry_utilities,

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

        polyhedron_data.EdgesDirection,
        polyhedron_data.EdgesTangent,
    };

    for (unsigned int k = 1; k < 4; k++)
    {
        Polydim::VEM::PCC::VEM_PCC_2D_ReferenceElement vem_reference_element_2D;
        Polydim::VEM::PCC::VEM_PCC_3D_ReferenceElement vem_reference_element_3D;
        Polydim::VEM::PCC::VEM_PCC_3D_Inertia_LocalSpace vem_local_space;

        const auto reference_element_data_2D = vem_reference_element_2D.Create(k);
        const auto reference_element_data_3D = vem_reference_element_3D.Create(k);
        const auto local_space = vem_local_space.CreateLocalSpace(
            reference_element_data_2D, reference_element_data_3D, polygonalFaces, polyhedron);

        // Test VEM performances
        Polydim::VEM::PCC::VEM_PCC_PerformanceAnalysis performanceAnalysis;

        const auto result = performanceAnalysis.Compute(Polydim::VEM::Monomials::VEM_Monomials_3D(),
                                                        reference_element_data_3D.Monomials,
                                                        vem_local_space,
                                                        local_space);

        ASSERT_TRUE(
            geometry_utilities.IsValueGreaterOrEqual(9.0e-14, result.ErrorPiNabla, geometry_utilities.Tolerance1D()));
        ASSERT_TRUE(
            geometry_utilities.IsValueGreaterOrEqual(7.4e-14, result.ErrorPi0km1, geometry_utilities.Tolerance1D()));
        ASSERT_TRUE(
            geometry_utilities.IsValueGreaterOrEqual(3.5e-13, result.ErrorPi0k, geometry_utilities.Tolerance1D()));
        ASSERT_EQ(result.ErrorPi0km1Grad.size(), 3);
        for (unsigned int d = 0; d < 3; ++d)
            ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(
                8e-13, result.ErrorPi0km1Grad[d], geometry_utilities.Tolerance1D()));
        ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(
            6.5e-12, result.ErrorStabilization, geometry_utilities.Tolerance1D()));
    }
}

} // namespace UnitTesting
} // namespace Polydim

#endif
