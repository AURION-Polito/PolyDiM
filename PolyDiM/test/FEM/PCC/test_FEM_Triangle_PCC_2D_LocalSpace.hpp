#ifndef __TEST_FEM_Triangle_PCC_2D_LocalSpace_H
#define __TEST_FEM_Triangle_PCC_2D_LocalSpace_H

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "FEM_Triangle_PCC_2D_LocalSpace.hpp"
#include "FEM_Triangle_PCC_2D_ReferenceElement.hpp"
#include "GeometryUtilities.hpp"

namespace Polydim
{
namespace UnitTesting
{
TEST(Test_FEM_Triangle_PCC_2D, Test_FEM_Triangle_PCC_2D_Reference_Element)
{
    const Polydim::FEM::PCC::FEM_Triangle_PCC_2D_ReferenceElement reference_element;

    const auto referenceQuadrature = Gedim::Quadrature::Quadrature_Gauss2D_Triangle::FillPointsAndWeights(21);
    const Eigen::MatrixXd &referenceQuadraturePoints = referenceQuadrature.Points;

    for (unsigned int o = 0; o < 5; o++)
    {
        const auto reference_element_data = reference_element.Create(o);

        const Eigen::MatrixXd dofs = reference_element_data.DofPositions;

        Eigen::MatrixXd points(3, dofs.cols() + referenceQuadraturePoints.cols());
        points << dofs, referenceQuadraturePoints;

        const Eigen::MatrixXd basisValues = reference_element.EvaluateBasisFunctions(points, reference_element_data);
        const std::vector<Eigen::MatrixXd> gradBasisValues =
            reference_element.EvaluateBasisFunctionDerivatives(points, reference_element_data);
        const std::array<Eigen::MatrixXd, 4> secondDerBasisValues =
            reference_element.EvaluateBasisFunctionSecondDerivatives(points, reference_element_data);

        const Eigen::VectorXd sumBasisValues = basisValues.rowwise().sum();
        const Eigen::VectorXd sumGradXValues = gradBasisValues[0].rowwise().sum();
        const Eigen::VectorXd sumGradYValues = gradBasisValues[1].rowwise().sum();
        const Eigen::VectorXd sumDerXXValues = secondDerBasisValues[0].rowwise().sum();
        const Eigen::VectorXd sumDerXYValues = secondDerBasisValues[1].rowwise().sum();
        const Eigen::VectorXd sumDerYXValues = secondDerBasisValues[2].rowwise().sum();
        const Eigen::VectorXd sumDerYYValues = secondDerBasisValues[3].rowwise().sum();

        for (unsigned int q = 0; q < points.cols(); q++)
        {
            ASSERT_TRUE(abs(sumBasisValues[q] - 1.0) < 1.0e-14);
            ASSERT_TRUE(abs(sumGradXValues[q]) < 1.0e-14);
            ASSERT_TRUE(abs(sumGradYValues[q]) < 1.0e-14);
            ASSERT_TRUE(abs(sumDerXXValues[q]) < 1.0e-14);
            ASSERT_TRUE(abs(sumDerXYValues[q]) < 1.0e-14);
            ASSERT_TRUE(abs(sumDerYXValues[q]) < 1.0e-14);
            ASSERT_TRUE(abs(sumDerYYValues[q]) < 1.0e-14);
        }
    }
}

TEST(Test_FEM_Triangle_PCC_2D, Test_FEM_Triangle_PCC_2D_Local_Space)
{
    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = 1.0e-8;
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    const auto poligon_vertices = geometry_utilities.CreateTriangle(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                    Eigen::Vector3d(1.0, 0.0, 0.0),
                                                                    Eigen::Vector3d(0.0, 1.0, 0.0));
    const std::vector<bool> polygon_edges_direction(3, true);
    const auto polygon_edges_tangent = geometry_utilities.PolygonEdgeTangents(poligon_vertices);
    const auto polygon_edges_length = geometry_utilities.PolygonEdgeLengths(poligon_vertices);
    const auto polygon_edges_normal = geometry_utilities.PolygonEdgeNormals(poligon_vertices);

    Polydim::FEM::PCC::FEM_Triangle_PCC_2D_Polygon_Geometry polygon_geometry;

    polygon_geometry = {geometry_utilities.Tolerance1D(),
                        geometry_utilities.Tolerance2D(),
                        poligon_vertices,
                        polygon_edges_direction,
                        polygon_edges_tangent,
                        polygon_edges_length};

    const Polydim::FEM::PCC::FEM_Triangle_PCC_2D_ReferenceElement reference_element;
    const Polydim::FEM::PCC::FEM_Triangle_PCC_2D_LocalSpace local_space;

    for (unsigned int o = 0; o < 5; o++)
    {
        const auto reference_element_data = reference_element.Create(o);
        const auto local_space_data = local_space.CreateLocalSpace(reference_element_data, polygon_geometry);

        const auto &internal_quadrature = local_space_data.InternalQuadrature;

        const Eigen::MatrixXd dofs = local_space_data.Dofs;

        Eigen::MatrixXd points(3, dofs.cols() + internal_quadrature.Points.cols());
        points << dofs, internal_quadrature.Points;

        const auto basisValues = local_space.ComputeBasisFunctionsValues(reference_element_data, local_space_data, points);
        const auto gradBasisValues =
            local_space.ComputeBasisFunctionsDerivativeValues(reference_element_data, local_space_data, points);

        const Eigen::VectorXd sumBasisValues = basisValues.rowwise().sum();
        const Eigen::VectorXd sumGradXValues = gradBasisValues[0].rowwise().sum();
        const Eigen::VectorXd sumGradYValues = gradBasisValues[1].rowwise().sum();

        for (unsigned int q = 0; q < points.cols(); q++)
        {
            ASSERT_TRUE(abs(sumBasisValues[q] - 1.0) < 1.0e-14);
            ASSERT_TRUE(abs(sumGradXValues[q]) < 1.0e-14);
            ASSERT_TRUE(abs(sumGradYValues[q]) < 1.0e-14);
        }

        const auto &derivative_values = local_space.ComputeBasisFunctionsDerivativeValues(reference_element_data,
                                                                                          local_space_data,
                                                                                          internal_quadrature.Points);

        Eigen::VectorXd internal_integral = Eigen::VectorXd::Zero(reference_element_data.NumBasisFunctions);
        for (unsigned int dim = 0; dim < reference_element_data.Dimension; ++dim)
            internal_integral += derivative_values[dim].transpose() * internal_quadrature.Weights;

        Eigen::VectorXd boundary_integral = Eigen::VectorXd::Zero(reference_element_data.NumBasisFunctions);
        for (unsigned int b = 0; b < polygon_edges_normal.cols(); ++b)
        {
            const Eigen::Vector3d boundary_normal = polygon_edges_normal.col(b);
            const auto &boundary_quadrature = local_space_data.BoundaryQuadrature[b];
            const auto boundary_values = reference_element.EvaluateBasisFunctions(boundary_quadrature.Points, reference_element_data);
            boundary_integral += boundary_values.transpose() * boundary_quadrature.Weights * boundary_normal.sum();
        }
        ASSERT_TRUE((internal_integral - boundary_integral).norm() < 1.0e-14 * std::max(1.0, boundary_integral.norm()));
    }
}

} // namespace UnitTesting
} // namespace Polydim

#endif
