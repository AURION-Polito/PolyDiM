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
        case 1: {
            result.NumDofs0D = 1;
            result.NumDofs1D = 0;
            result.NumDofs2D = 0;
            result.NumDofs3D = 0;
            result.NumBasisFunctions = 4;

            result.DofPositions.resize(result.Dimension, result.NumBasisFunctions);
            result.DofPositions.col(0) << 0.0, 0.0, 0.0;
            result.DofPositions.col(1) << 1.0, 0.0, 0.0;
            result.DofPositions.col(2) << 0.0, 1.0, 0.0;
            result.DofPositions.col(3) << 0.0, 0.0, 1.0;
        }
        break;
        default:
            throw std::runtime_error("order " + std::to_string(order) + "not supported yet");
        }

        result.ReferenceTetrahedronQuadrature =
            Gedim::Quadrature::Quadrature_Gauss3D_Tetrahedron_PositiveWeights::FillPointsAndWeights(2 * order);
        result.ReferenceTriangleQuadrature = Gedim::Quadrature::Quadrature_Gauss2D_Triangle::FillPointsAndWeights(2 * order);
        result.ReferenceSegmentQuadrature = Gedim::Quadrature::Quadrature_Gauss1D::FillPointsAndWeights(2 * order);

        result.ReferenceBasisFunctionValues = EvaluateBasisFunctions(result.ReferenceTetrahedronQuadrature.Points, result);
        result.ReferenceBasisFunctionDerivativeValues =
            EvaluateBasisFunctionDerivatives(result.ReferenceTetrahedronQuadrature.Points, result);

        return result;
    }
    // ***************************************************************************
    Eigen::MatrixXd EvaluateLambda(const Eigen::MatrixXd &points) const
    {
        Eigen::MatrixXd lambda;
        lambda.setZero(points.cols(), 4);

        lambda.col(0) = 1.0 - points.row(0).array() - points.row(1).array() - points.row(2).array();
        lambda.col(1) = points.row(0);
        lambda.col(2) = points.row(1);
        lambda.col(3) = points.row(2);

        return lambda;
    }
    // ***************************************************************************
    std::array<Eigen::MatrixXd, 3> EvaluateGradLambda(const Eigen::MatrixXd &points) const
    {
        std::array<Eigen::MatrixXd, 3> gradLambda;

        gradLambda[0].setZero(points.cols(), 4);
        gradLambda[1].setZero(points.cols(), 4);
        gradLambda[2].setZero(points.cols(), 4);

        gradLambda[0].col(0).setConstant(-1.0);
        gradLambda[1].col(0).setConstant(-1.0);
        gradLambda[2].col(0).setConstant(-1.0);

        gradLambda[0].col(1).setOnes();
        gradLambda[1].col(1).setZero();
        gradLambda[2].col(1).setZero();

        gradLambda[0].col(2).setZero();
        gradLambda[1].col(2).setOnes();
        gradLambda[2].col(2).setZero();

        gradLambda[0].col(3).setZero();
        gradLambda[1].col(3).setZero();
        gradLambda[2].col(3).setOnes();

        return gradLambda;
    }
    // ***************************************************************************
    Eigen::MatrixXd EvaluateBasisFunctions(const Eigen::MatrixXd &points,
                                           const FEM_Tetrahedron_PCC_3D_ReferenceElement_Data &reference_element_data) const
    {
        switch (reference_element_data.Order)
        {
        case 1:
            return EvaluateLambda(points);
        default:
            throw std::runtime_error("order " + std::to_string(reference_element_data.Order) + "not supported yet");
        }
    }
    // ***************************************************************************
    std::vector<Eigen::MatrixXd> EvaluateBasisFunctionDerivatives(const Eigen::MatrixXd &points,
                                                                  const FEM_Tetrahedron_PCC_3D_ReferenceElement_Data &reference_element_data) const
    {
        std::vector<Eigen::MatrixXd> values(reference_element_data.Dimension);
        for (unsigned int d = 0; d < reference_element_data.Dimension; d++)
            values[d].resize(points.cols(), reference_element_data.NumBasisFunctions);

        const auto grad_lambda = EvaluateGradLambda(points);

        switch (reference_element_data.Order)
        {
        case 1:
            for (unsigned int i = 0; i < reference_element_data.Dimension; i++)
                values[i] = grad_lambda[i];

            return values;
        default:
            throw std::runtime_error("order " + std::to_string(reference_element_data.Order) + "not supported yet");
        }
    }
};
} // namespace PCC
} // namespace FEM
} // namespace Polydim

#endif
