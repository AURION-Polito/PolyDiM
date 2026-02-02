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

#ifndef __FEM_PCC_1D_ReferenceElement_HPP
#define __FEM_PCC_1D_ReferenceElement_HPP

#include "Eigen/Eigen"
#include "IOUtilities.hpp"
#include "QuadratureData.hpp"
#include "Quadrature_Gauss1D.hpp"
#include "Quadrature_GaussLobatto1D.hpp"
#include "lagrange_1D.hpp"

namespace Polydim
{
namespace FEM
{
namespace PCC
{

enum class FEM_PCC_1D_Types
{
    Equispaced = 0,
    GaussLobatto = 1
};

struct FEM_PCC_1D_ReferenceElement_Data final
{
    unsigned int Dimension;
    unsigned int Order;
    unsigned int NumDofs0D;
    unsigned int NumDofs1D;

    unsigned int NumBasisFunctions;
    Eigen::MatrixXd DofPositions;
    Eigen::VectorXd Interpolation_coefficients;

    Gedim::Quadrature::QuadratureData ReferenceSegmentQuadrature;

    Eigen::MatrixXd ReferenceBasisFunctionValues;
    std::vector<Eigen::MatrixXd> ReferenceBasisFunctionDerivativeValues;
};

class FEM_PCC_1D_ReferenceElement final
{
  public:
    Polydim::FEM::PCC::FEM_PCC_1D_ReferenceElement_Data Create(
        const unsigned int order,
        const Polydim::FEM::PCC::FEM_PCC_1D_Types type = Polydim::FEM::PCC::FEM_PCC_1D_Types::Equispaced,
        const unsigned int quadrature_order = 0) const
    {
        Polydim::FEM::PCC::FEM_PCC_1D_ReferenceElement_Data result;

        result.Dimension = 1;
        result.Order = order;

        unsigned int computed_quadrature_order = quadrature_order;
        if (computed_quadrature_order == 0)
          computed_quadrature_order = 2 * order;

        if (order == 0)
        {
            result.NumDofs0D = 0;
            result.NumDofs1D = 1;
            result.NumBasisFunctions = 1;
            result.DofPositions = (Eigen::MatrixXd(3, 1) << 0.5, 0.0, 0.0).finished();

            result.ReferenceSegmentQuadrature = Gedim::Quadrature::Quadrature_Gauss1D::FillPointsAndWeights(computed_quadrature_order);

            result.Interpolation_coefficients = Eigen::VectorXd::Zero(1);

            result.ReferenceBasisFunctionValues = EvaluateBasisFunctions(result.ReferenceSegmentQuadrature.Points, result);
            result.ReferenceBasisFunctionDerivativeValues =
                EvaluateBasisFunctionDerivatives(result.ReferenceSegmentQuadrature.Points, result);

            return result;
        }

        unsigned int vertexCount = 0;
        std::vector<unsigned int> nodeDofs;
        std::vector<unsigned int> cellDofs;
        std::vector<unsigned int> nodeDofsIndex(3, 0);

        result.NumDofs0D = 1;
        result.NumDofs1D = order - 1;
        result.NumBasisFunctions = order + 1;

        Eigen::MatrixXd dofPositions = Eigen::MatrixXd::Zero(3, result.NumBasisFunctions);

        switch (type)
        {
        case Polydim::FEM::PCC::FEM_PCC_1D_Types::Equispaced:
            dofPositions.row(0) = Eigen::VectorXd::LinSpaced(result.NumBasisFunctions, 0.0, 1.0);
            break;
        case Polydim::FEM::PCC::FEM_PCC_1D_Types::GaussLobatto:
            dofPositions.row(0) =
                Gedim::Quadrature::Quadrature_GaussLobatto1D::FillPointsAndWeights(2 * order - 1).Points.row(0);
            break;
        default:
            throw std::runtime_error("not valid fem 1D type");
        }

        Gedim::Output::Assert(dofPositions.row(0).size() == result.NumBasisFunctions);

        for (unsigned int i = 0; i < result.NumBasisFunctions; i++)
        {
            if (i == 0)
            {
                nodeDofs.push_back(vertexCount);
                nodeDofsIndex[1]++;
            }
            else if (i == result.NumBasisFunctions - 1)
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

        result.DofPositions.resize(3, result.NumBasisFunctions);

        unsigned int dofCounter = 0;
        for (unsigned int n = 0; n < 2; n++)
        {
            for (unsigned int i = nodeDofsIndex[n]; i < nodeDofsIndex[n + 1]; i++)
                result.DofPositions.col(dofCounter++) << dofPositions.col(nodeDofs[i]);
        }

        for (unsigned int c = 0; c < cellDofs.size(); c++)
            result.DofPositions.col(dofCounter++) << dofPositions.col(cellDofs[c]);

        result.Interpolation_coefficients =
            Polydim::Interpolation::Lagrange::Lagrange_1D_coefficients(result.DofPositions.row(0).transpose());

        result.ReferenceSegmentQuadrature = Gedim::Quadrature::Quadrature_Gauss1D::FillPointsAndWeights(computed_quadrature_order);

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
    Eigen::MatrixXd EvaluateBasisFunctions(const Eigen::MatrixXd &points,
                                           const Polydim::FEM::PCC::FEM_PCC_1D_ReferenceElement_Data &reference_element_data) const
    {
        return Polydim::Interpolation::Lagrange::Lagrange_1D_values(reference_element_data.DofPositions.row(0).transpose(),
                                                                    reference_element_data.Interpolation_coefficients,
                                                                    points.row(0).transpose());
    }
    // ***************************************************************************
    std::vector<Eigen::MatrixXd> EvaluateBasisFunctionDerivatives(const Eigen::MatrixXd &points,
                                                                  const Polydim::FEM::PCC::FEM_PCC_1D_ReferenceElement_Data &reference_element_data) const
    {
        std::vector<Eigen::MatrixXd> values(reference_element_data.Dimension);

        for (unsigned int d = 0; d < reference_element_data.Dimension; d++)
            values[d].setZero(points.cols(), reference_element_data.NumBasisFunctions);

        values[0] = Polydim::Interpolation::Lagrange::Lagrange_1D_derivative_values(
            reference_element_data.DofPositions.row(0).transpose(),
            reference_element_data.Interpolation_coefficients,
            points.row(0).transpose());

        return values;
    }
    // ***************************************************************************
};
} // namespace PCC
} // namespace FEM
} // namespace Polydim

#endif
