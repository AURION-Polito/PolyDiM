#ifndef __FEM_Triangle_PCC_2D_ReferenceElement_H
#define __FEM_Triangle_PCC_2D_ReferenceElement_H

#include "Eigen/Eigen"
#include "QuadratureData.hpp"
#include "Quadrature_Gauss1D.hpp"
#include "Quadrature_Gauss2D_Triangle.hpp"
#include "FEM_PCC_1D_ReferenceElement.hpp"

namespace Polydim
{
namespace FEM
{
namespace PCC
{
/// \brief Base class for storing information related to \ref VEM::PCC::I_VEM_PCC_2D_ReferenceElement
struct FEM_Triangle_PCC_2D_ReferenceElement_Data final
{
    unsigned int Dimension; ///< Geometric dimension
    unsigned int Order;     ///< Order of the method
    unsigned int NumDofs0D; ///< Number of dofs for each vertex.
    unsigned int NumDofs1D; ///< Number of dofs internal to each edge.
    unsigned int NumDofs2D; ///< Number of dofs internal to each polygon.

    unsigned int NumBasisFunctions; ///< Number of total basis functions
    Eigen::MatrixXd DofPositions;   ///< reference element dof points
    Eigen::MatrixXi DofTypes;       ///< dof type [num oblique edges, num vertical edges, num horizontal edges]

    Gedim::Quadrature::QuadratureData ReferenceTriangleQuadrature;

    Eigen::MatrixXd ReferenceBasisFunctionValues;
    std::vector<Eigen::MatrixXd> ReferenceBasisFunctionDerivativeValues;
    std::array<Eigen::MatrixXd, 4> ReferenceBasisFunctionSecondDerivativeValues;

    FEM_PCC_1D_ReferenceElement_Data BoundaryReferenceElement_Data;
};

class FEM_Triangle_PCC_2D_ReferenceElement final
{
  public:
    FEM_Triangle_PCC_2D_ReferenceElement_Data Create(const unsigned int order) const
    {
        FEM_Triangle_PCC_2D_ReferenceElement_Data result;

        if (order == 0)
        {
            result.Dimension = 2;
            result.Order = order;
            result.NumDofs0D = 0;
            result.NumDofs1D = 0;
            result.NumDofs2D = 1;
            result.NumBasisFunctions = 1;
            result.DofTypes.setZero(3, result.NumBasisFunctions);
            result.DofPositions.setZero(3, result.NumBasisFunctions);
            result.DofPositions.col(0) << 1.0 / 3.0, 1.0 / 3.0, 0.0;

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

        result.Dimension = 2;
        result.Order = order;
        result.NumBasisFunctions = (order + 1) * (order + 2) / 2;

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
        result.NumDofs0D = nodeDofs.size() / 3;
        result.NumDofs1D = edgeDofs[0].size();
        result.NumDofs2D = cellDofs.size();
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
                                           const FEM_Triangle_PCC_2D_ReferenceElement_Data &reference_element_data) const
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
                                                                  const FEM_Triangle_PCC_2D_ReferenceElement_Data &reference_element_data) const
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
                                                                          const FEM_Triangle_PCC_2D_ReferenceElement_Data &reference_element_data) const
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
