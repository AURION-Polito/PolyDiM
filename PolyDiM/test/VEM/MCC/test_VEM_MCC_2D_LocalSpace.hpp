#ifndef __TEST_VEM_MCC_2D_LocalSpace_H
#define __TEST_VEM_MCC_2D_LocalSpace_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "GeometryUtilities.hpp"
#include "VEM_MCC_2D_Velocity_LocalSpace.hpp"
#include "VEM_MCC_2D_Partial_Velocity_LocalSpace.hpp"
#include "VEM_MCC_2D_Ortho_Velocity_LocalSpace.hpp"
#include "VEM_MCC_PerformanceAnalysis.hpp"

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

Test_VEM_MCC_2D_Polygon_Geometry Test_VEM_MCC_2D_Geometry(const Gedim::GeometryUtilities& geometry_utilities)
{
    Test_VEM_MCC_2D_Polygon_Geometry result;

    Eigen::MatrixXd& polygonVertices = result.Vertices;
    polygonVertices.setZero(3, 8);
    polygonVertices.row(0)<< 2.8000000000000004e-03,
        3.7200000000000004e-02,
        3.6900000000000002e-02,
        3.3300000000000003e-02,
        2.5600000000000001e-02,
        1.4400000000000000e-02,
        6.7000000000000002e-03,
        3.0999999999999999e-03;
    polygonVertices.row(1)<< 2.1499999999999998e-02,
        2.1499999999999998e-02,
        2.9500000000000002e-02,
        4.0599999999999997e-02,
        5.0700000000000002e-02,
        5.0700000000000002e-02,
        4.0599999999999997e-02,
        2.9500000000000002e-02;

    result.Measure = 7.9890999999999990e-04;
    result.Diameter = 3.7046997179258676e-02;
    result.Centroid = Eigen::Vector3d(2.0000000000000004e-02, 3.4061346918509809e-02, 0.0);
    result.EdgesLength.resize(8);
    result.EdgesLength<< 3.4400000000000000e-02,
        8.0056230238501787e-03,
        1.1669190203266030e-02,
        1.2700393694685222e-02,
        1.1200000000000002e-02,
        1.2700393694685220e-02,
        1.1669190203266030e-02,
        8.0056230238501769e-03;

    result.EdgesDirection.resize(8, true);

    Eigen::MatrixXd& edgeTangents = result.EdgesTangent;
    edgeTangents.setZero(3, 8);
    edgeTangents.row(0)<< 3.4400000000000000e-02,
        -3.0000000000000165e-04,
        -3.5999999999999990e-03,
        -7.7000000000000020e-03,
        -1.1200000000000002e-02,
        -7.6999999999999994e-03,
        -3.6000000000000003e-03,
        -2.9999999999999949e-04;
    edgeTangents.row(1)<< 0.0000000000000000e+00,
        8.0000000000000036e-03,
        1.1099999999999995e-02,
        1.0100000000000005e-02,
        0.0000000000000000e+00,
        -1.0100000000000005e-02,
        -1.1099999999999995e-02,
        -8.0000000000000036e-03;

    Eigen::MatrixXd& edgeNormals = result.EdgesNormal;
    edgeNormals.setZero(3, 8);
    edgeNormals.row(0)<< 0.0000000000000000e+00,  9.9929761570918052e-01,  9.5122281894876248e-01,  7.9525093810490199e-01,  0.0000000000000000e+00, -7.9525093810490211e-01, -9.5122281894876248e-01, -9.9929761570918074e-01;
    edgeNormals.row(1)<< -1.0000000000000000e+00, 3.7473660589094460e-02, 3.0850469803743652e-01, 6.0628041815918254e-01, 1.0000000000000000e+00, 6.0628041815918243e-01, 3.0850469803743663e-01, 3.7473660589094196e-02;

    // Map triangle points
    std::vector<unsigned int> polygonTriangulation = geometry_utilities.PolygonTriangulationByFirstVertex(polygonVertices);
    result.TriangulationVertices = geometry_utilities.ExtractTriangulationPoints(polygonVertices,
                                                                                 polygonTriangulation);

    return result;
}

TEST(Test_VEM_MCC, Test_VEM_MCC_2D_O0_O1_O2_O3)
{
    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = std::numeric_limits<double>::epsilon();
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    const auto polygon_data = Test_VEM_MCC_2D_Geometry(geometry_utilities);

    Polydim::VEM::MCC::VEM_MCC_2D_Polygon_Geometry polygon =
        {
            polygon_data.Vertices,
            polygon_data.Centroid,
            polygon_data.Measure,
            polygon_data.Diameter,
            polygon_data.TriangulationVertices,
            polygon_data.EdgesLength,
            polygon_data.EdgesDirection,
            polygon_data.EdgesTangent,
            polygon_data.EdgesNormal
        };

    for(unsigned int k = 0; k < 4; k++)
    {
        Polydim::VEM::MCC::VEM_MCC_2D_Velocity_ReferenceElement vem_reference_element;
        Polydim::VEM::MCC::VEM_MCC_2D_Velocity_LocalSpace vem_local_space;

        const auto reference_element_data = vem_reference_element.Create(k);
        const auto local_space = vem_local_space.CreateLocalSpace(reference_element_data,
                                                                  polygon);

        // Test VEM performances
        Polydim::VEM::MCC::VEM_MCC_PerformanceAnalysis performanceAnalysis;

        const auto result = performanceAnalysis.Compute(polygon.Measure,
                                                        polygon.Diameter,
                                                        Polydim::VEM::Monomials::VEM_Monomials_2D(),
                                                        reference_element_data.MonomialsKp1,
                                                        vem_local_space,
                                                        local_space);

        ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(1.0e-10,
                                                             result.ErrorPi0k,
                                                             geometry_utilities.Tolerance1D()));

        ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(1.0e-10,
                                                             result.ErrorGBD,
                                                             geometry_utilities.Tolerance1D()));

        ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(1.0e-10,
                                                             result.ErrorStabilization,
                                                             geometry_utilities.Tolerance1D()));
    }
}

TEST(Test_VEM_MCC, Test_VEM_MCC_Partial_2D_O0_O1_O2_O3)
{
    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = std::numeric_limits<double>::epsilon();
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    const auto polygon_data = Test_VEM_MCC_2D_Geometry(geometry_utilities);

    Polydim::VEM::MCC::VEM_MCC_2D_Polygon_Geometry polygon =
        {
            polygon_data.Vertices,
            polygon_data.Centroid,
            polygon_data.Measure,
            polygon_data.Diameter,
            polygon_data.TriangulationVertices,
            polygon_data.EdgesLength,
            polygon_data.EdgesDirection,
            polygon_data.EdgesTangent,
            polygon_data.EdgesNormal
        };

    for(unsigned int k = 0; k < 4; k++)
    {
        Polydim::VEM::MCC::VEM_MCC_2D_Velocity_ReferenceElement vem_reference_element;
        Polydim::VEM::MCC::VEM_MCC_2D_Partial_Velocity_LocalSpace vem_local_space;

        const auto reference_element_data = vem_reference_element.Create(k);
        const auto local_space = vem_local_space.CreateLocalSpace(reference_element_data,
                                                                  polygon);

        // Test VEM performances
        Polydim::VEM::MCC::VEM_MCC_PerformanceAnalysis performanceAnalysis;

        const auto result = performanceAnalysis.Compute(polygon.Measure,
                                                        polygon.Diameter,
                                                        Polydim::VEM::Monomials::VEM_Monomials_2D(),
                                                        reference_element_data.MonomialsKp1,
                                                        vem_local_space,
                                                        local_space);

        ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(1.0e-10,
                                                             result.ErrorPi0k,
                                                             geometry_utilities.Tolerance1D()));

        ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(1.0e-10,
                                                             result.ErrorGBD,
                                                             geometry_utilities.Tolerance1D()));

        ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(1.0e-10,
                                                             result.ErrorStabilization,
                                                             geometry_utilities.Tolerance1D()));
    }
}

TEST(Test_VEM_MCC, Test_VEM_MCC_Ortho_2D_O0_O1_O2_O3)
{
    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = std::numeric_limits<double>::epsilon();
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    const auto polygon_data = Test_VEM_MCC_2D_Geometry(geometry_utilities);

    Polydim::VEM::MCC::VEM_MCC_2D_Polygon_Geometry polygon =
        {
            polygon_data.Vertices,
            polygon_data.Centroid,
            polygon_data.Measure,
            polygon_data.Diameter,
            polygon_data.TriangulationVertices,
            polygon_data.EdgesLength,
            polygon_data.EdgesDirection,
            polygon_data.EdgesTangent,
            polygon_data.EdgesNormal
        };

    for(unsigned int k = 0; k < 4; k++)
    {
        Polydim::VEM::MCC::VEM_MCC_2D_Velocity_ReferenceElement vem_reference_element;
        Polydim::VEM::MCC::VEM_MCC_2D_Ortho_Velocity_LocalSpace vem_local_space;

        const auto reference_element_data = vem_reference_element.Create(k);
        const auto local_space = vem_local_space.CreateLocalSpace(reference_element_data,
                                                                  polygon);

        // Test VEM performances
        Polydim::VEM::MCC::VEM_MCC_PerformanceAnalysis performanceAnalysis;

        const auto result = performanceAnalysis.Compute(polygon.Measure,
                                                        polygon.Diameter,
                                                        Polydim::VEM::Monomials::VEM_Monomials_2D(),
                                                        reference_element_data.MonomialsKp1,
                                                        vem_local_space,
                                                        local_space);

        ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(1.0e-10,
                                                             result.ErrorPi0k,
                                                             geometry_utilities.Tolerance1D()));

        ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(1.0e-10,
                                                             result.ErrorGBD,
                                                             geometry_utilities.Tolerance1D()));

        ASSERT_TRUE(geometry_utilities.IsValueGreaterOrEqual(1.0e-10,
                                                             result.ErrorStabilization,
                                                             geometry_utilities.Tolerance1D()));
    }
}

}
}

#endif
