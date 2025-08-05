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

#ifndef __FEM_Quadrilateral_PCC_2D_ReferenceElement_HPP
#define __FEM_Quadrilateral_PCC_2D_ReferenceElement_HPP

#include "Eigen/Eigen"
#include "FEM_PCC_1D_ReferenceElement.hpp"
#include "QuadratureData.hpp"
#include "Quadrature_Gauss1D.hpp"
#include "Quadrature_Gauss2D_Triangle.hpp"

namespace Polydim
{
namespace FEM
{
namespace PCC
{
struct FEM_Quadrilateral_PCC_2D_ReferenceElement_Data final
{
    unsigned int Dimension;
    unsigned int Order;
    unsigned int NumDofs0D;
    unsigned int NumDofs1D;
    unsigned int NumDofs2D;

    unsigned int NumBasisFunctions;
    Eigen::MatrixXd DofPositions;
    Eigen::MatrixXi DofTypes;

    Gedim::Quadrature::QuadratureData ReferenceTriangleQuadrature;
    Gedim::Quadrature::QuadratureData ReferenceSquareQuadrature;

    Eigen::MatrixXd ReferenceBasisFunctionValues;
    std::vector<Eigen::MatrixXd> ReferenceBasisFunctionDerivativeValues;
    std::array<Eigen::MatrixXd, 4> ReferenceBasisFunctionSecondDerivativeValues;

    FEM_PCC_1D_ReferenceElement_Data BoundaryReferenceElement_Data;
};

class FEM_Quadrilateral_PCC_2D_ReferenceElement final
{
  public:
    FEM_Quadrilateral_PCC_2D_ReferenceElement_Data Create(const unsigned int order) const
    {
        FEM_Quadrilateral_PCC_2D_ReferenceElement_Data result;

        result.Dimension = 2;
        result.Order = order;
        result.NumBasisFunctions = (order + 1) * (order + 1);
        result.NumDofs0D = 1;
        result.NumDofs1D = order - 1;
        result.NumDofs2D = (order - 1) * (order - 1);

        std::vector<unsigned int> nodeDofs = {0, order, result.NumBasisFunctions - 1};
        std::list<unsigned int> cellDofs;
        std::vector<std::list<unsigned int>> edgeDofs(3);

        Eigen::MatrixXd localDofPositions = Eigen::MatrixXd::Zero(3, result.NumBasisFunctions);
        Eigen::MatrixXi localDofTypes = Eigen::MatrixXi::Zero(3, result.NumBasisFunctions);

        const double h = 1.0 / order;
        unsigned int dof = 0;
        for (unsigned int i = 0; i < order + 1; i++)
        {
            for (unsigned int j = 0; j < order + 1 - i; j++)
            {
                if (i == 0)
                {
                    if (j > 0 && j < order - i)
                        edgeDofs[0].push_back(dof);
                }
                else if (i < order && (j == 0 || j == order - i))
                {
                    if (j == 0)
                        edgeDofs[2].push_front(dof);
                    else
                        edgeDofs[1].push_back(dof);
                }
                else if (i < order)
                    cellDofs.push_back(dof);

                localDofPositions.col(dof) << (double)j * h, (double)i * h, 0.0;
                localDofTypes.col(dof) << order - i - j, j, i;

                dof++;
            }
        }

        // Reordering Dofs using convention [point, edge, cell]

        result.DofPositions.setZero(3, result.NumBasisFunctions);
        result.DofTypes.setZero(3, result.NumBasisFunctions);

        dof = 0;
        for (const unsigned int dofIndex : nodeDofs)
        {
            result.DofPositions.col(dof) << localDofPositions.col(dofIndex);
            result.DofTypes.col(dof) << localDofTypes.col(dofIndex);
            dof++;
        }
        for (unsigned int e = 0; e < 3; e++)
        {
            for (const unsigned int dofIndex : edgeDofs.at(e))
            {
                result.DofPositions.col(dof) << localDofPositions.col(dofIndex);
                result.DofTypes.col(dof) << localDofTypes.col(dofIndex);
                dof++;
            }
        }
        for (const unsigned int dofIndex : cellDofs)
        {
            result.DofPositions.col(dof) << localDofPositions.col(dofIndex);
            result.DofTypes.col(dof) << localDofTypes.col(dofIndex);
            dof++;
        }

        FEM_PCC_1D_ReferenceElement boundary_reference_element;
        result.BoundaryReferenceElement_Data = boundary_reference_element.Create(order);

        result.ReferenceTriangleQuadrature = Gedim::Quadrature::Quadrature_Gauss2D_Triangle::FillPointsAndWeights(2 * order);

        result.ReferenceBasisFunctionValues = EvaluateBasisFunctions(result.ReferenceTriangleQuadrature.Points, result);
        result.ReferenceBasisFunctionDerivativeValues =
            EvaluateBasisFunctionDerivatives(result.ReferenceTriangleQuadrature.Points, result);
        result.ReferenceBasisFunctionSecondDerivativeValues =
            EvaluateBasisFunctionSecondDerivatives(result.ReferenceTriangleQuadrature.Points, result);

        return result;
    }

    Eigen::MatrixXd EvaluateBasisFunctions(const Eigen::MatrixXd &points,
                                           const FEM_Quadrilateral_PCC_2D_ReferenceElement_Data &reference_element_data) const
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
                                                                  const FEM_Quadrilateral_PCC_2D_ReferenceElement_Data &reference_element_data) const
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
    // ***************************************************************************
    std::array<Eigen::MatrixXd, 4> EvaluateBasisFunctionSecondDerivatives(const Eigen::MatrixXd &points,
                                                                          const FEM_Quadrilateral_PCC_2D_ReferenceElement_Data &reference_element_data) const
    {
        switch (reference_element_data.Order)
        {
        case 0:
        case 1: {
            const Eigen::MatrixXd zero_matrix = Eigen::MatrixXd::Zero(points.cols(), reference_element_data.NumBasisFunctions);
            return {zero_matrix, zero_matrix, zero_matrix, zero_matrix};
        }
        case 2: {
            std::array<Eigen::MatrixXd, 4> constant_laplacian;

            for (unsigned int der = 0; der < 4; ++der)
                constant_laplacian[der] = Eigen::MatrixXd::Zero(points.cols(), reference_element_data.NumBasisFunctions);

            constant_laplacian[0].col(0) = Eigen::VectorXd::Constant(points.cols(), +4.0);
            constant_laplacian[1].col(0) = Eigen::VectorXd::Constant(points.cols(), +4.0);
            constant_laplacian[2].col(0) = Eigen::VectorXd::Constant(points.cols(), +4.0);
            constant_laplacian[3].col(0) = Eigen::VectorXd::Constant(points.cols(), +4.0);

            constant_laplacian[0].col(1) = Eigen::VectorXd::Constant(points.cols(), +4.0);
            constant_laplacian[1].col(1) = Eigen::VectorXd::Constant(points.cols(), +0.0);
            constant_laplacian[2].col(1) = Eigen::VectorXd::Constant(points.cols(), +0.0);
            constant_laplacian[3].col(1) = Eigen::VectorXd::Constant(points.cols(), +0.0);

            constant_laplacian[0].col(2) = Eigen::VectorXd::Constant(points.cols(), +0.0);
            constant_laplacian[1].col(2) = Eigen::VectorXd::Constant(points.cols(), +0.0);
            constant_laplacian[2].col(2) = Eigen::VectorXd::Constant(points.cols(), +0.0);
            constant_laplacian[3].col(2) = Eigen::VectorXd::Constant(points.cols(), +4.0);

            constant_laplacian[0].col(3) = Eigen::VectorXd::Constant(points.cols(), -8.0);
            constant_laplacian[1].col(3) = Eigen::VectorXd::Constant(points.cols(), -4.0);
            constant_laplacian[2].col(3) = Eigen::VectorXd::Constant(points.cols(), -4.0);
            constant_laplacian[3].col(3) = Eigen::VectorXd::Constant(points.cols(), +0.0);

            constant_laplacian[0].col(4) = Eigen::VectorXd::Constant(points.cols(), +0.0);
            constant_laplacian[1].col(4) = Eigen::VectorXd::Constant(points.cols(), +4.0);
            constant_laplacian[2].col(4) = Eigen::VectorXd::Constant(points.cols(), +4.0);
            constant_laplacian[3].col(4) = Eigen::VectorXd::Constant(points.cols(), +0.0);

            constant_laplacian[0].col(5) = Eigen::VectorXd::Constant(points.cols(), +0.0);
            constant_laplacian[1].col(5) = Eigen::VectorXd::Constant(points.cols(), -4.0);
            constant_laplacian[2].col(5) = Eigen::VectorXd::Constant(points.cols(), -4.0);
            constant_laplacian[3].col(5) = Eigen::VectorXd::Constant(points.cols(), -8.0);

            return constant_laplacian;
        }
        default: {
            const Eigen::MatrixXd zero_matrix = Eigen::MatrixXd::Zero(points.cols(), reference_element_data.NumBasisFunctions);
            return {zero_matrix, zero_matrix, zero_matrix, zero_matrix};
        }
        }
    }
};
} // namespace PCC
} // namespace FEM
} // namespace Polydim

#endif
