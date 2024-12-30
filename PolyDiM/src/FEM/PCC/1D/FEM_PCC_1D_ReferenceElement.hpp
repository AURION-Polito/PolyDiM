#ifndef __FEM_PCC_1D_ReferenceElement_H
#define __FEM_PCC_1D_ReferenceElement_H

#include "Eigen/Eigen"
#include "QuadratureData.hpp"
#include "Quadrature_Gauss1D.hpp"
#include "Quadrature_Gauss2D_Triangle.hpp"

namespace Polydim
{
namespace FEM
{
namespace PCC
{
/// \brief Base class for storing information related to \ref VEM::PCC::I_VEM_PCC_2D_ReferenceElement
struct FEM_PCC_1D_ReferenceElement_Data final
{
    unsigned int Dimension; ///< Geometric dimension
    unsigned int Order;     ///< Order of the method
    unsigned int NumDofs0D; ///< Number of dofs for each vertex.
    unsigned int NumDofs1D; ///< Number of dofs internal to each edge.

    unsigned int NumBasisFunctions; ///< Number of total basis functions
    Eigen::MatrixXd DofPositions;   ///< reference element dof points

    Gedim::Quadrature::QuadratureData ReferenceSegmentQuadrature;

    Eigen::MatrixXd ReferenceBasisFunctionValues;
    std::vector<Eigen::MatrixXd> ReferenceBasisFunctionDerivativeValues;
};

class FEM_PCC_1D_ReferenceElement final
{
  public:
    FEM_PCC_1D_ReferenceElement_Data Create(const unsigned int order) const
    {
        FEM_PCC_1D_ReferenceElement_Data result;

        result.Dimension = 1;
        result.Order = order;

        if (order == 0)
        {
            result.NumDofs0D = 0;
            result.NumDofs1D = 1;
            result.NumBasisFunctions = 1;
            result.DofPositions = (Eigen::MatrixXd(3, 1) << 0.5, 0.0, 0.0).finished();

            result.ReferenceSegmentQuadrature = Gedim::Quadrature::Quadrature_Gauss1D::FillPointsAndWeights(2 * order);

            result.ReferenceBasisFunctionValues = EvaluateBasisFunctions(result.ReferenceSegmentQuadrature.Points, result);
            result.ReferenceBasisFunctionDerivativeValues =
                EvaluateBasisFunctionDerivatives(result.ReferenceSegmentQuadrature.Points, result);

            return result;
        }

        unsigned int vertexCount = 0;
        std::vector<unsigned int> nodeDofs;            ///< The degrees of freedom of the element on nodes
        std::vector<unsigned int> cellDofs;            ///< The degrees of freedom of the element on cell
        std::vector<unsigned int> nodeDofsIndex(3, 0); ///< Index of DOFS on nodes

        result.NumDofs0D = 0;
        result.NumDofs1D = 0;
        result.NumBasisFunctions = order + 1;

        Eigen::MatrixXd dofPositions = Eigen::MatrixXd::Zero(3, result.NumBasisFunctions);

        for (unsigned int i = 0; i < result.NumBasisFunctions; i++)
        {
            dofPositions(0, vertexCount) = ((double)i) / ((double)order);

            if (dofPositions(0, vertexCount) == 0.0)
            {
                nodeDofs.push_back(vertexCount);
                nodeDofsIndex[1]++;
            }
            else if (dofPositions(0, vertexCount) == 1.0)
            {
                nodeDofs.push_back(vertexCount);
                nodeDofsIndex[2]++;
            }
            else
                cellDofs.push_back(vertexCount);

            vertexCount++;
        }

        for (unsigned int i = 1; i < nodeDofsIndex.size(); i++)
            nodeDofsIndex[i] = nodeDofsIndex[i] + nodeDofsIndex[i - 1];

        if (vertexCount != result.NumBasisFunctions || nodeDofs.size() != nodeDofsIndex.back())
            throw std::runtime_error("Wrong initialization in FE reference element. Number of DOFs found is not "
                                     "correct.");

        /// <li> Reordering Dofs using convention [point, cell]
        result.NumDofs0D = nodeDofs.size() / 2;
        result.NumDofs1D = cellDofs.size();
        result.DofPositions.resize(3, result.NumBasisFunctions);

        unsigned int dofCounter = 0;
        for (unsigned int n = 0; n < 2; n++)
        {
            for (unsigned int i = nodeDofsIndex[n]; i < nodeDofsIndex[n + 1]; i++)
                result.DofPositions.col(dofCounter++) << dofPositions.col(nodeDofs[i]);
        }

        for (unsigned int c = 0; c < cellDofs.size(); c++)
            result.DofPositions.col(dofCounter++) << dofPositions.col(cellDofs[c]);

        result.ReferenceSegmentQuadrature = Gedim::Quadrature::Quadrature_Gauss1D::FillPointsAndWeights(2 * order);

        result.ReferenceBasisFunctionValues = EvaluateBasisFunctions(result.ReferenceSegmentQuadrature.Points, result);
        result.ReferenceBasisFunctionDerivativeValues =
            EvaluateBasisFunctionDerivatives(result.ReferenceSegmentQuadrature.Points, result);

        return result;
    }
    // ***************************************************************************
    Eigen::MatrixXd EvaluateLambda(const Eigen::MatrixXd &points) const
    {
        Eigen::MatrixXd lambda = Eigen::MatrixXd::Zero(points.cols(), 2);

        lambda.setZero(points.cols(), 2);

        lambda.col(0) = 1.0 - points.row(0).array();
        lambda.col(1) = points.row(0);

        return lambda;
    }
    // ***************************************************************************
    std::vector<Eigen::MatrixXd> EvaluateGradLambda(const Eigen::MatrixXd &points) const
    {
        std::vector<Eigen::MatrixXd> gradLambda(1, Eigen::MatrixXd::Zero(points.cols(), 2));

        gradLambda[0].col(0).setConstant(-1.0);
        gradLambda[0].col(1).setOnes();

        return gradLambda;
    }
    // ***************************************************************************
    Eigen::MatrixXd EvaluateBasisFunctions(const Eigen::MatrixXd &points, const FEM_PCC_1D_ReferenceElement_Data &reference_element_data) const
    {
        Eigen::MatrixXd values(points.cols(), reference_element_data.NumBasisFunctions);

        const auto lambda = EvaluateLambda(points);

        switch (reference_element_data.Order)
        {
        case 0:
            values.col(0) << Eigen::VectorXd::Constant(points.cols(), 1.0);
            break;
        case 1:
            values.col(0) << lambda.col(0);
            values.col(1) << lambda.col(1);
            break;
        case 2:
            values.col(0) << 2.0 * lambda.col(0).array() * (lambda.col(0).array() - 0.5);
            values.col(1) << 2.0 * lambda.col(1).array() * (lambda.col(1).array() - 0.5);
            values.col(2) << 4.0 * lambda.col(0).array() * lambda.col(1).array();
            break;
        default:
            throw std::runtime_error("Order not supported yet");
        }

        return values;
    }
    // ***************************************************************************
    std::vector<Eigen::MatrixXd> EvaluateBasisFunctionDerivatives(const Eigen::MatrixXd &points,
                                                                  const FEM_PCC_1D_ReferenceElement_Data &reference_element_data) const
    {
        std::vector<Eigen::MatrixXd> values(reference_element_data.Dimension);

        for (unsigned int d = 0; d < reference_element_data.Dimension; d++)
            values[d].resize(points.cols(), reference_element_data.NumBasisFunctions);

        const auto lambda = EvaluateLambda(points);
        const auto gradLambda = EvaluateGradLambda(points);

        switch (reference_element_data.Order)
        {
        case 0:
            values[0].col(0) = Eigen::VectorXd::Zero(points.cols());
            break;
        case 1:
            values[0].col(0) = gradLambda[0].col(0);
            values[0].col(1) = gradLambda[0].col(1);
            break;
        case 2:
            values[0].col(0) = 2.0 * gradLambda[0].col(0).array() * (lambda.col(0).array() - 0.5) +
                               2.0 * lambda.col(0).array() * gradLambda[0].col(0).array();
            values[0].col(1) = 2.0 * gradLambda[0].col(1).array() * (lambda.col(1).array() - 0.5) +
                               2.0 * lambda.col(1).array() * gradLambda[0].col(1).array();
            values[0].col(2) = 4.0 * gradLambda[0].col(0).array() * lambda.col(1).array() +
                               4.0 * lambda.col(0).array() * gradLambda[0].col(1).array();
            break;
        default:
            throw std::runtime_error("Order not supported yet");
        }

        return values;
    }
    // ***************************************************************************
};
} // namespace PCC
} // namespace FEM
} // namespace Polydim

#endif
