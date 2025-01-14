#ifndef __TEST_FEM_PCC_1D_LocalSpace_H
#define __TEST_FEM_PCC_1D_LocalSpace_H

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "FEM_PCC_1D_ReferenceElement.hpp"

namespace Polydim
{
namespace UnitTesting
{
TEST(Test_FEM_PCC_1D, Test_FEM_PCC_1D_Reference_Element)
{
    const Polydim::FEM::PCC::FEM_PCC_1D_ReferenceElement reference_element;

    const auto referenceQuadrature = Gedim::Quadrature::Quadrature_Gauss1D::FillPointsAndWeights(10);
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

        const Eigen::VectorXd sumBasisValues = basisValues.rowwise().sum();
        const Eigen::VectorXd sumGradXValues = gradBasisValues[0].rowwise().sum();
        for (unsigned int q = 0; q < points.cols(); q++)
        {
            ASSERT_TRUE(abs(sumBasisValues[q] - 1.0) < 1.0e-14);
            ASSERT_TRUE(abs(sumGradXValues[q]) < 1.0e-14);
        }
    }
}
} // namespace UnitTesting
} // namespace Polydim

#endif
