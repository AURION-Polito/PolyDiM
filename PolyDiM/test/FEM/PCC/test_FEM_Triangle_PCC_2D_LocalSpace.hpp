#ifndef __TEST_FEM_Triangle_PCC_2D_LocalSpace_H
#define __TEST_FEM_Triangle_PCC_2D_LocalSpace_H

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "FEM_Triangle_PCC_2D_ReferenceElement.hpp"

namespace Polydim
{
  namespace UnitTesting
  {
    TEST(Test_FEM_Triangle_PCC, Test_FEM_Triangle_PCC_2D_Reference_Element)
    {
      const Polydim::FEM::PCC::FEM_RefElement_Langrange_PCC_Triangle_2D reference_element;

      const auto referenceQuadrature = Gedim::Quadrature::Quadrature_Gauss2D_Triangle::FillPointsAndWeights(21);
      const Eigen::MatrixXd& referenceQuadraturePoints = referenceQuadrature.Points;

      for (unsigned int o = 0; o < 5; o++)
      {
        const auto reference_element_data = reference_element.Create(o);

        const Eigen::MatrixXd dofs = reference_element_data.DofPositions;

        Eigen::MatrixXd points(3,
                               dofs.cols() +
                               referenceQuadraturePoints.cols());
        points << dofs, referenceQuadraturePoints;

        const Eigen::MatrixXd basisValues = reference_element.EvaluateBasisFunctions(points,
                                                                                     reference_element_data);
        const std::vector<Eigen::MatrixXd> gradBasisValues = reference_element.EvaluateBasisFunctionDerivatives(points,
                                                                                                                reference_element_data);

        const Eigen::VectorXd sumBasisValues = basisValues.rowwise().sum();
        const Eigen::VectorXd sumGradXValues = gradBasisValues[0].rowwise().sum();
        const Eigen::VectorXd sumGradYValues = gradBasisValues[0].rowwise().sum();
        for (unsigned int q = 0; q < points.cols(); q++)
        {
          ASSERT_TRUE(abs(sumBasisValues[q] - 1.0) < 1.0e-14);
          ASSERT_TRUE(abs(sumGradXValues[q]) < 1.0e-14);
          ASSERT_TRUE(abs(sumGradYValues[q]) < 1.0e-14);
        }
      }
    }
  } // namespace UnitTesting
} // namespace Polydim

#endif
