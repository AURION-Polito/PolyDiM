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

    std::array<Eigen::Vector3d, 2> boundaries_normal;
    boundaries_normal[0] = Eigen::Vector3d(-1.0, 0.0, 0.0).normalized();
    boundaries_normal[1] = Eigen::Vector3d(+1.0, 0.0, 0.0).normalized();
    std::array<Gedim::Quadrature::QuadratureData, 2> boundaries_quadrature;
    boundaries_quadrature[0].Points = Eigen::Vector3d(0.0, 0.0, 0.0);
    boundaries_quadrature[0].Weights = Eigen::VectorXd::Constant(1, 1.0);
    boundaries_quadrature[1].Points = Eigen::Vector3d(1.0, 0.0, 0.0);
    boundaries_quadrature[1].Weights = Eigen::VectorXd::Constant(1, 1.0);

    for (unsigned int o = 0; o < 5; o++)
    {
        const auto reference_element_data = reference_element.Create(o);

        const Eigen::MatrixXd dofs = reference_element_data.DofPositions;

        Eigen::MatrixXd points(3, dofs.cols() + referenceQuadraturePoints.cols());
        points << dofs, referenceQuadraturePoints;

        const Eigen::MatrixXd basisValues = reference_element.EvaluateBasisFunctions(points, reference_element_data);
        ASSERT_TRUE((basisValues.topRows(reference_element_data.NumBasisFunctions) -
                     Eigen::MatrixXd::Identity(reference_element_data.NumBasisFunctions, reference_element_data.NumBasisFunctions))
                        .norm() < 1.0e-13);

        const std::vector<Eigen::MatrixXd> gradBasisValues =
            reference_element.EvaluateBasisFunctionDerivatives(points, reference_element_data);

        const Eigen::VectorXd sumBasisValues = basisValues.rowwise().sum();
        const Eigen::VectorXd sumGradXValues = gradBasisValues[0].rowwise().sum();
        for (unsigned int q = 0; q < points.cols(); q++)
        {
            ASSERT_TRUE(abs(sumBasisValues[q] - 1.0) < 1.0e-14);
            ASSERT_TRUE(abs(sumGradXValues[q]) < 1.0e-14);
        }

        const auto &internal_quadrature = reference_element_data.ReferenceSegmentQuadrature;
        const auto &derivative_values = reference_element_data.ReferenceBasisFunctionDerivativeValues;

        Eigen::VectorXd internal_integral = Eigen::VectorXd::Zero(reference_element_data.NumBasisFunctions);
        for (unsigned int dim = 0; dim < reference_element_data.Dimension; ++dim)
            internal_integral += derivative_values[dim].transpose() * internal_quadrature.Weights;

        Eigen::VectorXd boundary_integral = Eigen::VectorXd::Zero(reference_element_data.NumBasisFunctions);
        for (unsigned int b = 0; b < boundaries_normal.size(); ++b)
        {
            const Eigen::Vector3d boundary_normal = boundaries_normal[b];
            const auto &boundary_quadrature = boundaries_quadrature[b];
            const auto boundary_values = reference_element.EvaluateBasisFunctions(boundary_quadrature.Points, reference_element_data);
            boundary_integral += boundary_values.transpose() * boundary_quadrature.Weights * boundary_normal.sum();
        }
        ASSERT_TRUE((internal_integral - boundary_integral).norm() < 1.0e-14 * std::max(1.0, boundary_integral.norm()));
    }
}
} // namespace UnitTesting
} // namespace Polydim

#endif
