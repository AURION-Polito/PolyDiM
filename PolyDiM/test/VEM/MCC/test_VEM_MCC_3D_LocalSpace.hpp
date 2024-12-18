#ifndef __TEST_VEM_MCC_3D_LocalSpace_H
#define __TEST_VEM_MCC_3D_LocalSpace_H

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "GeometryUtilities.hpp"
#include "VEM_MCC_3D_Velocity_LocalSpace.hpp"
#include "VEM_MCC_PerformanceAnalysis.hpp"

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

Test_VEM_MCC_3D_Polyhedron_Geometry Test_VEM_MCC_3D_Geometry(const Gedim::GeometryUtilities &geometry_utilities)
{
    Test_VEM_MCC_3D_Polyhedron_Geometry result;

    const Eigen::Vector3d origin = Eigen::Vector3d(0.0, 0.0, 0.0);
    const Eigen::Vector3d length = Eigen::Vector3d(1.0, 0.0, 0.0);
    const Eigen::Vector3d height = Eigen::Vector3d(0.0, 0.0, 1.0);
    const Eigen::Vector3d width = Eigen::Vector3d(0.0, 1.0, 0.0);
    Gedim::GeometryUtilities::Polyhedron cube =
        geometry_utilities.CreateParallelepipedWithOrigin(origin, length, height, width);

    result.Vertices = cube.Vertices;

    result.Measure = 1.0e+00;
    result.Diameter = sqrt(3.0);
    result.Centroid << 0.5, 0.5, 0.5;

    const std::vector<Eigen::MatrixXd> facesVertices =
        geometry_utilities.PolyhedronFaceVertices(result.Vertices, cube.Faces);
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
        const Eigen::MatrixXd facesVertices2D = geometry_utilities.RotatePointsFrom3DTo2D(
            facesVertices[f], result.FacesRotationMatrix[f].transpose(), result.FacesTranslation[f]);
        std::vector<unsigned int> facesTriangulation =
            geometry_utilities.PolygonTriangulationByFirstVertex(facesVertices2D);
        result.FacesTriangulationVertices2D[f] =
            geometry_utilities.ExtractTriangulationPoints(facesVertices2D, facesTriangulation);
        result.FacesMeasure[f] = geometry_utilities.PolygonArea(facesVertices2D);
        result.FacesDiameter[f] = geometry_utilities.PolygonDiameter(facesVertices2D);
        result.FacesCentroid2D[f] = geometry_utilities.PolygonCentroid(facesVertices2D, result.FacesMeasure[f]);
    }

    std::vector<std::vector<unsigned int>> facesTriangulation3D =
        geometry_utilities.PolyhedronFaceTriangulationsByFirstVertex(cube.Faces, facesVertices);

    std::vector<unsigned int> polyhedronTetrahedrons = geometry_utilities.PolyhedronTetrahedronsByFaceTriangulations(
        result.Vertices, cube.Faces, facesTriangulation3D, result.Centroid);

    result.TetrahedronVertices =
        geometry_utilities.ExtractTetrahedronPoints(result.Vertices, result.Centroid, polyhedronTetrahedrons);

    return result;
}

TEST(Test_VEM_MCC, Test_VEM_MCC_3D_O0_O1_O2_O3)
{
    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = std::numeric_limits<double>::epsilon();
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    const Test_VEM_MCC_3D_Polyhedron_Geometry polyhedron_data = Test_VEM_MCC_3D_Geometry(geometry_utilities);

    Polydim::VEM::MCC::VEM_MCC_3D_Polyhedron_Geometry polyhedron = {geometry_utilities,

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
        Polydim::VEM::MCC::VEM_MCC_3D_Velocity_LocalSpace vem_local_space;

        const auto reference_element_data = vem_reference_element.Create(k);
        const auto local_space = vem_local_space.CreateLocalSpace(reference_element_data, polyhedron);

        // Test VEM performances
        Polydim::VEM::MCC::VEM_MCC_PerformanceAnalysis performanceAnalysis;

        const auto result = performanceAnalysis.Compute(polyhedron.Measure,
                                                        polyhedron.Diameter,
                                                        Polydim::VEM::Monomials::VEM_Monomials_3D(),
                                                        reference_element_data.MonomialsKp1,
                                                        vem_local_space,
                                                        local_space);

        ASSERT_TRUE(
            geometry_utilities.IsValueGreaterOrEqual(1.0e-10, result.ErrorPi0k, geometry_utilities.Tolerance1D()));

        ASSERT_TRUE(
            geometry_utilities.IsValueGreaterOrEqual(1.0e-10, result.ErrorGBD, geometry_utilities.Tolerance1D()));

        ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(
            1.0e-10, result.ErrorStabilization, geometry_utilities.Tolerance1D()));
    }
}

} // namespace UnitTesting
} // namespace Polydim

#endif
