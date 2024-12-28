#ifndef __FEM_Tetrahedron_PCC_3D_ReferenceElement_H
#define __FEM_Tetrahedron_PCC_3D_ReferenceElement_H

#include "Eigen/Eigen"
#include "QuadratureData.hpp"
#include "Quadrature_Gauss1D.hpp"
#include "Quadrature_Gauss2D_Triangle.hpp"
#include "Quadrature_Gauss3D_Tetrahedron_PositiveWeights.hpp"

namespace Polydim
{
  namespace FEM
  {
    namespace PCC
    {
      struct FEM_Tetrahedron_PCC_3D_ReferenceElement_Data final
      {
          unsigned int Dimension; ///< Geometric dimension
          unsigned int Order;     ///< Order of the method
          unsigned int NumDofs0D; ///< Number of dofs for each vertex.
          unsigned int NumDofs1D; ///< Number of dofs internal to each edge.
          unsigned int NumDofs2D; ///< Number of dofs internal to each polygon.
          unsigned int NumDofs3D; ///< Number of dofs internal to each polyhedron.

          unsigned int NumBasisFunctions; ///< Number of total basis functions
          Eigen::MatrixXd DofPositions;   ///< reference element dof points

          Gedim::Quadrature::QuadratureData ReferenceSegmentQuadrature;
          Gedim::Quadrature::QuadratureData ReferenceTriangleQuadrature;
          Gedim::Quadrature::QuadratureData ReferenceTetrahedronQuadrature;

          Eigen::MatrixXd ReferenceBasisFunctionValues;
          std::vector<Eigen::MatrixXd> ReferenceBasisFunctionDerivativeValues;
      };

      class FEM_RefElement_Langrange_PCC_Tetrahedron_3D final
      {
        public:
          FEM_Tetrahedron_PCC_3D_ReferenceElement_Data Create(const unsigned int order) const
          {
            FEM_Tetrahedron_PCC_3D_ReferenceElement_Data result;

            result.Order = order;
            result.Dimension = 3;

            switch (order)
            {
              case 1:
              {
                result.NumDofs0D = 1;
                result.NumDofs1D = 0;
                result.NumDofs2D = 0;
                result.NumDofs3D = 0;
                result.NumBasisFunctions = 4;

                result.DofPositions.resize(result.Dimension,
                                           result.NumBasisFunctions);
                result.DofPositions.col(0)<< 0.0, 0.0, 0.0;
                result.DofPositions.col(1)<< 1.0, 0.0, 0.0;
                result.DofPositions.col(2)<< 0.0, 1.0, 0.0;
                result.DofPositions.col(3)<< 0.0, 0.0, 1.0;
              }
                break;
              default:
                throw std::runtime_error("order " + std::to_string(order) + "not supported yet");
            }

            result.ReferenceTetrahedronQuadrature = Gedim::Quadrature::Quadrature_Gauss3D_Tetrahedron_PositiveWeights::FillPointsAndWeights(2 * order);
            result.ReferenceTriangleQuadrature = Gedim::Quadrature::Quadrature_Gauss2D_Triangle::FillPointsAndWeights(2 * order);
            result.ReferenceSegmentQuadrature = Gedim::Quadrature::Quadrature_Gauss1D::FillPointsAndWeights(2 * order);

            result.ReferenceBasisFunctionValues = EvaluateBasisFunctions(result.ReferenceTetrahedronQuadrature.Points, result);
            result.ReferenceBasisFunctionDerivativeValues =
                EvaluateBasisFunctionDerivatives(result.ReferenceTetrahedronQuadrature.Points, result);

            return result;
          }

          Eigen::MatrixXd EvaluateBasisFunctions(const Eigen::MatrixXd &points,
                                                 const FEM_Tetrahedron_PCC_3D_ReferenceElement_Data &reference_element_data) const
          {
            switch (reference_element_data.Order)
            {
              case 0:
                return Eigen::VectorXd::Constant(points.cols(), 1.0);
              default: {
                const double h = 1.0 / reference_element_data.Order;
                const Eigen::ArrayXd x = points.row(0).transpose().array();
                const Eigen::ArrayXd y = points.row(1).transpose().array();
                Eigen::MatrixXd values = Eigen::MatrixXd::Ones(points.cols(), reference_element_data.NumBasisFunctions);

                for (unsigned int d = 0; d < reference_element_data.NumBasisFunctions; d++)
                {
                  const Eigen::Vector3i &dofType = reference_element_data.DofTypes.col(d);
                  const Eigen::Vector3d &dofPosition = reference_element_data.DofPositions.col(d);

                  // terms of equation 1 - x - y - t * h
                  for (unsigned int t = 0; t < static_cast<unsigned int>(dofType[0]); t++)
                  {
                    values.col(d).array() *= (1.0 - x - y - t * h);
                    values.col(d) /= (1.0 - dofPosition.x() - dofPosition.y() - t * h);
                  }

                  // terms of equation x - t * h
                  for (unsigned int t = 0; t < static_cast<unsigned int>(dofType[1]); t++)
                  {
                    values.col(d).array() *= (x - t * h);
                    values.col(d) /= (dofPosition.x() - t * h);
                  }

                  // terms of equation y - t * h
                  for (unsigned int t = 0; t < static_cast<unsigned int>(dofType[2]); t++)
                  {
                    values.col(d).array() *= (y - t * h);
                    values.col(d) /= (dofPosition.y() - t * h);
                  }
                }
                return values;
              }
            }
          }
          // ***************************************************************************
          std::vector<Eigen::MatrixXd> EvaluateBasisFunctionDerivatives(const Eigen::MatrixXd &points,
                                                                        const FEM_Tetrahedron_PCC_3D_ReferenceElement_Data &reference_element_data) const
          {
            switch (reference_element_data.Order)
            {
              case 0:
                return std::vector<Eigen::MatrixXd>(reference_element_data.Dimension,
                                                    Eigen::MatrixXd::Zero(points.cols(), reference_element_data.NumBasisFunctions));
              default: {
                const double h = 1.0 / reference_element_data.Order;
                const Eigen::ArrayXd x = points.row(0).transpose().array();
                const Eigen::ArrayXd y = points.row(1).transpose().array();
                std::vector<Eigen::MatrixXd> gradValues(reference_element_data.Dimension,
                                                        Eigen::MatrixXd::Zero(points.cols(), reference_element_data.NumBasisFunctions));

                for (unsigned int d = 0; d < reference_element_data.NumBasisFunctions; d++)
                {
                  const Eigen::Vector3i &dofType = reference_element_data.DofTypes.col(d);
                  const Eigen::Vector3d &dofPosition = reference_element_data.DofPositions.col(d);

                  const unsigned int numProds = dofType[0] + dofType[1] + dofType[2];

                  std::vector<Eigen::ArrayXd> prod_terms(numProds);
                  std::vector<Eigen::Array2d> grad_terms(numProds);
                  double denominator = 1.0;

                  unsigned int dt = 0;
                  // terms of equation 1 - x - y - t * h
                  for (unsigned int t = 0; t < static_cast<unsigned int>(dofType[0]); t++)
                  {
                    prod_terms[dt] = (1.0 - x - y - t * h);
                    grad_terms[dt] << -1.0, -1.0;
                    denominator *= (1.0 - dofPosition.x() - dofPosition.y() - t * h);
                    dt++;
                  }

                  // terms of equation x - t * h
                  for (unsigned int t = 0; t < static_cast<unsigned int>(dofType[1]); t++)
                  {
                    prod_terms[dt] = (x - t * h);
                    grad_terms[dt] << 1.0, 0.0;
                    denominator *= (dofPosition.x() - t * h);
                    dt++;
                  }

                  // terms of equation y - t * h
                  for (unsigned int t = 0; t < static_cast<unsigned int>(dofType[2]); t++)
                  {
                    prod_terms[dt] = (y - t * h);
                    grad_terms[dt] << 0.0, 1.0;
                    denominator *= (dofPosition.y() - t * h);
                    dt++;
                  }

                  for (unsigned int i = 0; i < numProds; i++)
                  {
                    Eigen::ArrayXd inner_prod = Eigen::ArrayXd::Ones(points.cols());
                    for (unsigned int j = 0; j < numProds; j++)
                    {
                      if (i != j)
                        inner_prod *= prod_terms[j];
                    }

                    gradValues[0].col(d).array() += inner_prod * grad_terms[i][0];
                    gradValues[1].col(d).array() += inner_prod * grad_terms[i][1];
                  }

                  gradValues[0].col(d) /= denominator;
                  gradValues[1].col(d) /= denominator;
                }

                return gradValues;
              }
            }
          }
      };
    } // namespace PCC
  } // namespace FEM
} // namespace Polydim

#endif
