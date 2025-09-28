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

#ifndef __FEM_Tetrahedron_PCC_3D_ReferenceElement_HPP
#define __FEM_Tetrahedron_PCC_3D_ReferenceElement_HPP

#include "Eigen/Eigen"
#include "FEM_Triangle_PCC_2D_ReferenceElement.hpp"
#include "IOStream.hpp"
#include "QuadratureData.hpp"
#include "Quadrature_Gauss3D_Tetrahedron_PositiveWeights.hpp"

namespace Polydim
{
namespace FEM
{
namespace PCC
{
struct FEM_Tetrahedron_PCC_3D_ReferenceElement_Data final
{
    unsigned int Dimension;
    unsigned int Order;
    unsigned int NumDofs0D;
    unsigned int NumDofs1D;
    unsigned int NumDofs2D;
    unsigned int NumDofs3D;

    std::map<std::pair<unsigned int, unsigned int>, std::pair<unsigned int, bool>> Edges_by_vertices;
    std::map<std::pair<unsigned int, unsigned int>, unsigned int> Faces_by_edge_vertex;
    std::map<std::pair<unsigned int, unsigned int>, std::pair<std::array<unsigned int, 2>, bool>> Faces_by_edges;

    unsigned int NumBasisFunctions;
    Eigen::MatrixXd DofPositions;
    std::vector<std::array<unsigned int, 4>> DofTypes;

    Gedim::Quadrature::QuadratureData ReferenceTetrahedronQuadrature;

    Eigen::MatrixXd ReferenceBasisFunctionValues;
    std::vector<Eigen::MatrixXd> ReferenceBasisFunctionDerivativeValues;

    FEM_PCC_1D_ReferenceElement_Data EdgeReferenceElement_Data;
    FEM_Triangle_PCC_2D_ReferenceElement_Data BoundaryReferenceElement_Data;
};

class FEM_Tetrahedron_PCC_3D_ReferenceElement final
{
  public:
    FEM_Tetrahedron_PCC_3D_ReferenceElement_Data Create(const unsigned int order) const
    {
        FEM_Tetrahedron_PCC_3D_ReferenceElement_Data result;

        result.Order = order;
        result.Dimension = 3;

        result.Edges_by_vertices = {{{0, 1}, {0, true}},
                                    {{1, 0}, {0, false}},
                                    {{1, 2}, {1, true}},
                                    {{2, 1}, {1, false}},
                                    {{2, 0}, {2, true}},
                                    {{0, 2}, {2, false}},
                                    {{0, 3}, {3, true}},
                                    {{3, 0}, {3, false}},
                                    {{1, 3}, {4, false}},
                                    {{3, 1}, {4, true}},
                                    {{2, 3}, {5, false}},
                                    {{3, 2}, {5, true}}};

        result.Faces_by_edge_vertex = {{{0, 2}, 0},
                                       {{1, 0}, 0},
                                       {{2, 1}, 0},
                                       {{0, 3}, 1},
                                       {{4, 0}, 1},
                                       {{3, 1}, 1},
                                       {{2, 3}, 2},
                                       {{5, 0}, 2},
                                       {{3, 2}, 2},
                                       {{1, 3}, 3},
                                       {{5, 1}, 3},
                                       {{4, 2}, 3}};

        result.Faces_by_edges = {{{0, 1}, {std::array<unsigned int, 2>({0, 0}), true}},
                                 {{1, 2}, {std::array<unsigned int, 2>({0, 1}), true}},
                                 {{2, 0}, {std::array<unsigned int, 2>({0, 2}), true}},
                                 {{1, 0}, {std::array<unsigned int, 2>({0, 2}), false}},
                                 {{2, 1}, {std::array<unsigned int, 2>({0, 0}), false}},
                                 {{0, 2}, {std::array<unsigned int, 2>({0, 1}), false}},
                                 {{0, 4}, {std::array<unsigned int, 2>({1, 0}), true}},
                                 {{4, 3}, {std::array<unsigned int, 2>({1, 1}), true}},
                                 {{3, 0}, {std::array<unsigned int, 2>({1, 2}), true}},
                                 {{4, 0}, {std::array<unsigned int, 2>({1, 2}), false}},
                                 {{3, 4}, {std::array<unsigned int, 2>({1, 0}), false}},
                                 {{0, 3}, {std::array<unsigned int, 2>({1, 1}), false}},
                                 {{2, 5}, {std::array<unsigned int, 2>({2, 0}), true}},
                                 {{5, 3}, {std::array<unsigned int, 2>({2, 1}), true}},
                                 {{3, 2}, {std::array<unsigned int, 2>({2, 2}), true}},
                                 {{5, 2}, {std::array<unsigned int, 2>({2, 2}), false}},
                                 {{3, 5}, {std::array<unsigned int, 2>({2, 0}), false}},
                                 {{2, 3}, {std::array<unsigned int, 2>({2, 1}), false}},
                                 {{1, 5}, {std::array<unsigned int, 2>({3, 0}), true}},
                                 {{5, 4}, {std::array<unsigned int, 2>({3, 1}), true}},
                                 {{4, 1}, {std::array<unsigned int, 2>({3, 2}), true}},
                                 {{5, 1}, {std::array<unsigned int, 2>({3, 2}), false}},
                                 {{4, 5}, {std::array<unsigned int, 2>({3, 0}), false}},
                                 {{1, 4}, {std::array<unsigned int, 2>({3, 1}), false}}};

        if (order == 0)
        {
            result.NumDofs0D = 0;
            result.NumDofs1D = 0;
            result.NumDofs2D = 0;
            result.NumDofs3D = 1;
            result.NumBasisFunctions = 1;
            result.DofPositions.setZero(3, result.NumBasisFunctions);
            result.DofPositions.col(0) << 1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0;

            FEM_Triangle_PCC_2D_ReferenceElement boundary_reference_element;
            result.BoundaryReferenceElement_Data = boundary_reference_element.Create(order);

            result.ReferenceTetrahedronQuadrature =
                Gedim::Quadrature::Quadrature_Gauss3D_Tetrahedron_PositiveWeights::FillPointsAndWeights(2 * order);

            result.ReferenceBasisFunctionValues = EvaluateBasisFunctions(result.ReferenceTetrahedronQuadrature.Points, result);
            result.ReferenceBasisFunctionDerivativeValues =
                EvaluateBasisFunctionDerivatives(result.ReferenceTetrahedronQuadrature.Points, result);

            return result;
        }

        result.NumDofs0D = 1;
        result.NumDofs1D = order - 1;
        result.NumDofs2D = order > 2 ? (order - 2) * (order - 1) / 2 : 0;
        result.NumDofs3D = order > 3 ? (order - 3) * (order - 2) * (order - 1) / 6 : 0;
        result.NumBasisFunctions = 4 * result.NumDofs0D + 6 * result.NumDofs1D + 4 * result.NumDofs2D + result.NumDofs3D;

        result.DofTypes.resize(result.NumBasisFunctions);
        result.DofPositions.resize(result.Dimension, result.NumBasisFunctions);
        result.DofPositions.col(0) << 0.0, 0.0, 0.0;
        result.DofPositions.col(1) << 1.0, 0.0, 0.0;
        result.DofPositions.col(2) << 0.0, 1.0, 0.0;
        result.DofPositions.col(3) << 0.0, 0.0, 1.0;

        result.DofTypes[0] = {order, 0, 0, 0};
        result.DofTypes[1] = {0, order, 0, 0};
        result.DofTypes[2] = {0, 0, order, 0};
        result.DofTypes[3] = {0, 0, 0, order};

        if (order > 1)
        {
            unsigned int dof = 4;

            // edge x
            for (unsigned int d = 0; d < result.NumDofs1D; d++)
            {
                result.DofPositions.col(dof) << (d + 1.0) / ((double)order), 0.0, 0.0;
                result.DofTypes[dof] = {order - 1 - d, d + 1, 0, 0};
                dof++;
            }

            // edge x - y
            for (unsigned int d = 0; d < result.NumDofs1D; d++)
            {
                result.DofPositions.col(dof) << (result.NumDofs1D - d) / ((double)order), (d + 1.0) / ((double)order), 0.0;
                result.DofTypes[dof] = {0, order - d - 1, d + 1, 0};
                dof++;
            }

            // edge y
            for (unsigned int d = 0; d < result.NumDofs1D; d++)
            {
                result.DofPositions.col(dof) << 0.0, (result.NumDofs1D - d) / ((double)order), 0.0;
                result.DofTypes[dof] = {d + 1, 0, order - d - 1, 0};
                dof++;
            }

            // edge z
            for (unsigned int d = 0; d < result.NumDofs1D; d++)
            {
                result.DofPositions.col(dof) << 0.0, 0.0, (d + 1.0) / ((double)order);
                result.DofTypes[dof] = {order - d - 1, 0, 0, d + 1};
                dof++;
            }

            // edge x - z
            for (unsigned int d = 0; d < result.NumDofs1D; d++)
            {
                result.DofPositions.col(dof) << (d + 1.0) / ((double)order), 0.0, (result.NumDofs1D - d) / ((double)order);
                result.DofTypes[dof] = {0, d + 1, 0, order - d - 1};
                dof++;
            }

            // edge y - z
            for (unsigned int d = 0; d < result.NumDofs1D; d++)
            {
                result.DofPositions.col(dof) << 0.0, (d + 1.0) / ((double)order), (result.NumDofs1D - d) / ((double)order);
                result.DofTypes[dof] = {0, 0, d + 1, order - d - 1};
                dof++;
            }

            if (order > 2)
            {
                using namespace Gedim;
                // face x - y
                std::cout << "0 xy" << std::endl;
                for (unsigned int d1 = 0; d1 < result.NumDofs1D - 1; d1++) // x
                {
                    for (unsigned int d2 = 0; d2 < result.NumDofs1D - 1 - d1; d2++) // y
                    {
                        std::cout << "\t{ " << dof << ", ";
                        result.DofPositions.col(dof) << (d1 + 1.0) / ((double)order), (d2 + 1.0) / ((double)order), 0.0;
                        std::cout << "coord " << std::scientific << result.DofPositions.col(dof).transpose() << ", ";
                        result.DofTypes[dof] = {order - d1 - d2 - 2, d1 + 1, d2 + 1, 0};
                        std::cout << "types " << result.DofTypes[dof] << "} " << std::endl;
                        dof++;
                    }
                }

                // face x - z
                std::cout << "1 xz" << std::endl;
                for (unsigned int d1 = 0; d1 < result.NumDofs1D - 1; d1++) // x
                {
                    for (unsigned int d2 = 0; d2 < result.NumDofs1D - 1 - d1; d2++) // z
                    {
                        std::cout << "\t{ " << dof << ", ";
                        result.DofPositions.col(dof) << (d1 + 1.0) / ((double)order), 0.0, (d2 + 1.0) / ((double)order);
                        std::cout << "coord " << std::scientific << result.DofPositions.col(dof).transpose() << ", ";
                        result.DofTypes[dof] = {order - d1 - d2 - 2, d1 + 1, 0, d2 + 1};
                        std::cout << "types " << result.DofTypes[dof] << "} " << std::endl;
                        dof++;
                    }
                }

                // face y - z
                std::cout << "2 yz" << std::endl;
                for (unsigned int d1 = 0; d1 < result.NumDofs1D - 1; d1++) // y
                {
                    for (unsigned int d2 = 0; d2 < result.NumDofs1D - 1 - d1; d2++) // z
                    {
                        std::cout << "\t{ " << dof << ", ";
                        result.DofPositions.col(dof) << 0.0, (d1 + 1.0) / ((double)order), (d2 + 1.0) / ((double)order);
                        std::cout << "coord " << std::scientific << result.DofPositions.col(dof).transpose() << ", ";
                        result.DofTypes[dof] = {order - d1 - d2 - 2, 0, d1 + 1, d2 + 1};
                        std::cout << "types " << result.DofTypes[dof] << "} " << std::endl;
                        dof++;
                    }
                }

                // face x - y - z
                std::cout << "3 xyz" << std::endl;
                for (unsigned int d1 = 0; d1 < result.NumDofs1D - 1; d1++) // x
                {
                    for (unsigned int d2 = 0; d2 < result.NumDofs1D - 1 - d1; d2++) // y
                    {
                        std::cout << "\t{ " << dof << ", ";
                        result.DofPositions.col(dof) << (d1 + 1.0) / ((double)order), (d2 + 1.0) / ((double)order),
                            1.0 - (d1 + 1.0) / ((double)order) - (d2 + 1.0) / ((double)order);
                        std::cout << "coord " << std::scientific << result.DofPositions.col(dof).transpose() << ", ";
                        result.DofTypes[dof] = {0, d1 + 1, d2 + 1, order - d1 - d2 - 2};
                        std::cout << "types " << result.DofTypes[dof] << "} " << std::endl;
                        dof++;
                    }
                }

                if (order > 3)
                {
                    // internal
                    for (unsigned int d3 = 0; d3 < result.NumDofs1D - 2; d3++)
                    {
                        for (unsigned int d1 = 0; d1 < result.NumDofs1D - 2 - d3; d1++) // x
                        {
                            for (unsigned int d2 = 0; d2 < result.NumDofs1D - 2 - d1 - d3; d2++) // y
                            {
                                result.DofPositions.col(dof) << (d1 + 1.0) / ((double)order),
                                    (d2 + 1.0) / ((double)order), (d3 + 1.0) / ((double)order);

                                result.DofTypes[dof] = {order - d1 - d2 - d3 - 3, d1 + 1, d2 + 1, d3 + 1};
                                dof++;
                            }
                        }
                    }
                }
            }
        }

        FEM_PCC_1D_ReferenceElement edge_reference_element;
        result.EdgeReferenceElement_Data = edge_reference_element.Create(order, FEM_PCC_1D_Types::Equispaced);

        FEM_Triangle_PCC_2D_ReferenceElement boundary_reference_element;
        result.BoundaryReferenceElement_Data = boundary_reference_element.Create(order);

        result.ReferenceTetrahedronQuadrature =
            Gedim::Quadrature::Quadrature_Gauss3D_Tetrahedron_PositiveWeights::FillPointsAndWeights(2 * order);

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
    std::vector<Eigen::MatrixXd> EvaluateGradLambda(const Eigen::MatrixXd &points) const
    {
        std::vector<Eigen::MatrixXd> gradLambda(3);

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
        const Eigen::MatrixXd lambda_functions = EvaluateLambda(points);
        const unsigned int order = reference_element_data.Order;
        if (order == 1)
            return lambda_functions;

        Eigen::MatrixXd basis_functions_values = Eigen::MatrixXd::Ones(points.cols(), reference_element_data.NumBasisFunctions);
        Eigen::VectorXd normalized_factor = Eigen::VectorXd::Ones(reference_element_data.NumBasisFunctions);

        const Eigen::MatrixXd dofs_lambda_functions = EvaluateLambda(reference_element_data.DofPositions);

        for (unsigned int dof = 0; dof < reference_element_data.NumBasisFunctions; dof++)
        {
            for (unsigned int i = 0; i < 4; i++)
            {
                for (unsigned int p = 0; p < reference_element_data.DofTypes[dof][i]; p++)
                {
                    normalized_factor(dof) = normalized_factor(dof) * (dofs_lambda_functions(dof, i) - p / ((double)order));
                    basis_functions_values.col(dof) =
                        basis_functions_values.col(dof).array() * (lambda_functions.col(i).array() - p / ((double)order));
                }
            }

            basis_functions_values.col(dof) /= normalized_factor(dof);
        }

        return basis_functions_values;
    }
    // ***************************************************************************
    std::vector<Eigen::MatrixXd> EvaluateBasisFunctionDerivatives(const Eigen::MatrixXd &points,
                                                                  const FEM_Tetrahedron_PCC_3D_ReferenceElement_Data &reference_element_data) const
    {

        const std::vector<Eigen::MatrixXd> grad_lambda = EvaluateGradLambda(points);
        const unsigned int order = reference_element_data.Order;
        if (order == 1)
            return grad_lambda;

        std::vector<Eigen::MatrixXd> values(reference_element_data.Dimension,
                                            Eigen::MatrixXd::Zero(points.cols(), reference_element_data.NumBasisFunctions));

        const auto basis_function_values = EvaluateBasisFunctions(points, reference_element_data);
        const Eigen::MatrixXd lambda_functions = EvaluateLambda(points);

        for (unsigned int d = 0; d < reference_element_data.Dimension; d++)
        {
            for (unsigned int dof = 0; dof < reference_element_data.NumBasisFunctions; dof++)
            {
                for (unsigned int i = 0; i < 4; i++)
                {
                    for (unsigned int p = 0; p < reference_element_data.DofTypes[dof][i]; p++)
                    {
                        values[d].col(dof) = values[d].col(dof).array() +
                                             grad_lambda[d].col(i).array() * basis_function_values.col(dof).array() /
                                                 (lambda_functions.col(i).array() - p / ((double)order));
                    }
                }
            }
        }

        return values;
    }
};
} // namespace PCC
} // namespace FEM
} // namespace Polydim

#endif
