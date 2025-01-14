#ifndef __TEST_FEM_Tetrahedron_PCC_3D_LocalSpace_H
#define __TEST_FEM_Tetrahedron_PCC_3D_LocalSpace_H

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "FEM_Tetrahedron_PCC_3D_ReferenceElement.hpp"

namespace Polydim
{
namespace UnitTesting
{
TEST(Test_FEM_Tetrahedron_PCC_3D, Test_FEM_Tetrahedron_PCC_3D_Reference_Element)
{
    const Polydim::FEM::PCC::FEM_Tetrahedron_PCC_3D_ReferenceElement reference_element;

    const auto referenceQuadrature = Gedim::Quadrature::Quadrature_Gauss3D_Tetrahedron_PositiveWeights::FillPointsAndWeights(10);
    const Eigen::MatrixXd &referenceQuadraturePoints = referenceQuadrature.Points;

    for (unsigned int o = 1; o < 3; o++)
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
        const Eigen::VectorXd sumGradYValues = gradBasisValues[1].rowwise().sum();
        const Eigen::VectorXd sumGradZValues = gradBasisValues[2].rowwise().sum();
        for (unsigned int q = 0; q < points.cols(); q++)
        {
            ASSERT_TRUE((basisValues.topRows(reference_element_data.NumBasisFunctions) -
                         Eigen::MatrixXd::Identity(reference_element_data.NumBasisFunctions, reference_element_data.NumBasisFunctions))
                            .norm() < 1.0e-14);
            ASSERT_TRUE(abs(sumBasisValues[q] - 1.0) < 1.0e-14);
            ASSERT_TRUE(abs(sumGradXValues[q]) < 1.0e-14);
            ASSERT_TRUE(abs(sumGradYValues[q]) < 1.0e-14);
            ASSERT_TRUE(abs(sumGradZValues[q]) < 1.0e-14);
        }
    }
}
} // namespace UnitTesting
} // namespace Polydim

#endif
