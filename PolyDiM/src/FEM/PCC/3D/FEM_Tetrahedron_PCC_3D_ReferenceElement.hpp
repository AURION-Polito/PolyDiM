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
    unsigned int Dimension;
    unsigned int Order;
    unsigned int NumDofs0D;
    unsigned int NumDofs1D;
    unsigned int NumDofs2D;
    unsigned int NumDofs3D;

    std::map<std::pair<unsigned int, unsigned int>, std::pair<unsigned int, bool>> Edges_by_vertices;
    std::map<std::pair<unsigned int, unsigned int>, unsigned int> Faces_by_edge_vertex;

    unsigned int NumBasisFunctions;
    Eigen::MatrixXd DofPositions;

    Gedim::Quadrature::QuadratureData ReferenceTetrahedronQuadrature;

    Eigen::MatrixXd ReferenceBasisFunctionValues;
    std::vector<Eigen::MatrixXd> ReferenceBasisFunctionDerivativeValues;

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
        case 2: {
            result.NumDofs0D = 1;
            result.NumDofs1D = 1;
            result.NumDofs2D = 0;
            result.NumDofs3D = 0;
            result.NumBasisFunctions = 10;

            result.DofPositions.resize(result.Dimension, result.NumBasisFunctions);
            result.DofPositions.col(0) << 0.0, 0.0, 0.0;
            result.DofPositions.col(1) << 1.0, 0.0, 0.0;
            result.DofPositions.col(2) << 0.0, 1.0, 0.0;
            result.DofPositions.col(3) << 0.0, 0.0, 1.0;
            result.DofPositions.col(4) << 0.5, 0.0, 0.0;
            result.DofPositions.col(5) << 0.5, 0.5, 0.0;
            result.DofPositions.col(6) << 0.0, 0.5, 0.0;
            result.DofPositions.col(7) << 0.0, 0.0, 0.5;
            result.DofPositions.col(8) << 0.5, 0.0, 0.5;
            result.DofPositions.col(9) << 0.0, 0.5, 0.5;
        }
        break;
        case 3: {
            result.NumDofs0D = 1;
            result.NumDofs1D = 2;
            result.NumDofs2D = 1;
            result.NumDofs3D = 0;
            result.NumBasisFunctions = 20;

            result.DofPositions.resize(result.Dimension, result.NumBasisFunctions);
            result.DofPositions.col(0) << 0.0, 0.0, 0.0;
            result.DofPositions.col(1) << 1.0, 0.0, 0.0;
            result.DofPositions.col(2) << 0.0, 1.0, 0.0;
            result.DofPositions.col(3) << 0.0, 0.0, 1.0;
            result.DofPositions.col(4) << 1.0 / 3.0, 0.0, 0.0;
            result.DofPositions.col(5) << 2.0 / 3.0, 0.0, 0.0;
            result.DofPositions.col(6) << 2.0 / 3.0, 1.0 / 3.0, 0.0;
            result.DofPositions.col(7) << 1.0 / 3.0, 2.0 / 3.0, 0.0;
            result.DofPositions.col(8) << 0.0, 2.0 / 3.0, 0.0;
            result.DofPositions.col(9) << 0.0, 1.0 / 3.0, 0.0;
            result.DofPositions.col(10) << 0.0, 0.0, 1.0 / 3.0;
            result.DofPositions.col(11) << 0.0, 0.0, 2.0 / 3.0;
            result.DofPositions.col(12) << 1.0 / 3.0, 0.0, 2.0 / 3.0;
            result.DofPositions.col(13) << 2.0 / 3.0, 0.0, 1.0 / 3.0;
            result.DofPositions.col(14) << 0.0, 1.0 / 3.0, 2.0 / 3.0;
            result.DofPositions.col(15) << 0.0, 2.0 / 3.0, 1.0 / 3.0;
            result.DofPositions.col(16) << 1.0 / 3.0, 1.0 / 3.0, 0.0;
            result.DofPositions.col(17) << 1.0 / 3.0, 0.0, 1.0 / 3.0;
            result.DofPositions.col(18) << 0.0, 1.0 / 3.0, 1.0 / 3.0;
            result.DofPositions.col(19) << 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0;
        }
        break;
        default:
            throw std::runtime_error("order " + std::to_string(order) + "not supported yet");
        }

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

        const Eigen::MatrixXd lambda_functions = EvaluateLambda(points);
        switch (reference_element_data.Order)
        {
        case 1:
            return lambda_functions;
        case 2: {
            Eigen::MatrixXd basis_functions_values = Eigen::MatrixXd::Zero(points.cols(), reference_element_data.NumBasisFunctions);

            const Eigen::ArrayXd xyz = lambda_functions.col(0);
            const Eigen::ArrayXd x = lambda_functions.col(1);
            const Eigen::ArrayXd y = lambda_functions.col(2);
            const Eigen::ArrayXd z = lambda_functions.col(3);

            basis_functions_values.col(0) = 2.0 * xyz * (xyz - 0.5);
            basis_functions_values.col(1) = 2.0 * x * (x - 0.5);
            basis_functions_values.col(2) = 2.0 * y * (y - 0.5);
            basis_functions_values.col(3) = 2.0 * z * (z - 0.5);
            basis_functions_values.col(4) = 4.0 * xyz * x;
            basis_functions_values.col(5) = 4.0 * x * y;
            basis_functions_values.col(6) = 4.0 * y * xyz;
            basis_functions_values.col(7) = 4.0 * xyz * z;
            basis_functions_values.col(8) = 4.0 * x * z;
            basis_functions_values.col(9) = 4.0 * y * z;
            return basis_functions_values;
        }
        case 3: {
            Eigen::MatrixXd basis_functions_values = Eigen::MatrixXd::Zero(points.cols(), reference_element_data.NumBasisFunctions);

            const Eigen::ArrayXd xyz = lambda_functions.col(0);
            const Eigen::ArrayXd x = lambda_functions.col(1);
            const Eigen::ArrayXd y = lambda_functions.col(2);
            const Eigen::ArrayXd z = lambda_functions.col(3);

            basis_functions_values.col(0) = 9.0 / 2.0 * xyz * (xyz - 1.0 / 3.0) * (xyz - 2.0 / 3.0);
            basis_functions_values.col(1) = 9.0 / 2.0 * x * (x - 1.0 / 3.0) * (x - 2.0 / 3.0);
            basis_functions_values.col(2) = 9.0 / 2.0 * y * (y - 1.0 / 3.0) * (y - 2.0 / 3.0);
            basis_functions_values.col(3) = 9.0 / 2.0 * z * (z - 1.0 / 3.0) * (z - 2.0 / 3.0);
            basis_functions_values.col(4) = 27.0 / 2.0 * x * xyz * (xyz - 1.0 / 3.0);
            basis_functions_values.col(5) = 27.0 / 2.0 * x * xyz * (x - 1.0 / 3.0);
            basis_functions_values.col(6) = 27.0 / 2.0 * y * x * (x - 1.0 / 3.0);
            basis_functions_values.col(7) = 27.0 / 2.0 * x * y * (y - 1.0 / 3.0);
            basis_functions_values.col(8) = 27.0 / 2.0 * y * (y - 1.0 / 3.0) * xyz;
            basis_functions_values.col(9) = 27.0 / 2.0 * y * xyz * (xyz - 1.0 / 3.0);
            basis_functions_values.col(10) = 27.0 / 2.0 * z * (xyz - 1.0 / 3.0) * xyz;
            basis_functions_values.col(11) = 27.0 / 2.0 * z * xyz * (z - 1.0 / 3.0);
            basis_functions_values.col(12) = 27.0 / 2.0 * x * z * (z - 1.0 / 3.0);
            basis_functions_values.col(13) = 27.0 / 2.0 * z * x * (x - 1.0 / 3.0);
            basis_functions_values.col(14) = 27.0 / 2.0 * y * z * (z - 1.0 / 3.0);
            basis_functions_values.col(15) = 27.0 / 2.0 * z * y * (y - 1.0 / 3.0);
            basis_functions_values.col(16) = 27.0 * x * y * xyz;
            basis_functions_values.col(17) = 27.0 * x * z * xyz;
            basis_functions_values.col(18) = 27.0 * y * z * xyz;
            basis_functions_values.col(19) = 27.0 * x * y * z;
            return basis_functions_values;
        }
        default:
            throw std::runtime_error("order " + std::to_string(reference_element_data.Order) + "not supported yet");
        }
    }
    // ***************************************************************************
    std::vector<Eigen::MatrixXd> EvaluateBasisFunctionDerivatives(const Eigen::MatrixXd &points,
                                                                  const FEM_Tetrahedron_PCC_3D_ReferenceElement_Data &reference_element_data) const
    {
        std::vector<Eigen::MatrixXd> values(reference_element_data.Dimension,
                                            Eigen::MatrixXd::Zero(points.cols(), reference_element_data.NumBasisFunctions));
        const auto grad_lambda = EvaluateGradLambda(points);

        switch (reference_element_data.Order)
        {
        case 1: {
            for (unsigned int i = 0; i < reference_element_data.Dimension; i++)
                values[i] = grad_lambda[i];
        }
        break;
        case 2: {
            const Eigen::MatrixXd lambda_functions = EvaluateLambda(points);

            const Eigen::ArrayXd xyz = lambda_functions.col(0);
            const Eigen::ArrayXd x = lambda_functions.col(1);
            const Eigen::ArrayXd y = lambda_functions.col(2);
            const Eigen::ArrayXd z = lambda_functions.col(3);

            values[0].col(0) = -4.0 * xyz + 1.0;
            values[0].col(1) = 4.0 * x - 1.0;
            values[0].col(4) = -4.0 * x + 4.0 * xyz;
            values[0].col(5) = 4.0 * y;
            values[0].col(6) = -4.0 * y;
            values[0].col(7) = -4.0 * z;
            values[0].col(8) = 4.0 * z;

            values[1].col(0) = -4.0 * xyz + 1.0;
            values[1].col(2) = 4.0 * y - 1.0;
            values[1].col(4) = -4.0 * x;
            values[1].col(5) = 4.0 * x;
            values[1].col(6) = -4.0 * y + 4.0 * xyz;
            values[1].col(7) = -4.0 * z;
            values[1].col(9) = 4.0 * z;

            values[2].col(0) = -4.0 * xyz + 1.0;
            values[2].col(3) = 4.0 * z - 1.0;
            values[2].col(4) = -4.0 * x;
            values[2].col(6) = -4.0 * y;
            values[2].col(7) = 4.0 * xyz - 4.0 * z;
            values[2].col(8) = 4.0 * x;
            values[2].col(9) = 4.0 * y;
        }
        break;
        case 3: {
            const Eigen::MatrixXd lambda_functions = EvaluateLambda(points);

            const Eigen::ArrayXd xyz = lambda_functions.col(0);
            const Eigen::ArrayXd x = lambda_functions.col(1);
            const Eigen::ArrayXd y = lambda_functions.col(2);
            const Eigen::ArrayXd z = lambda_functions.col(3);

            values[0].col(0) = -9.0 / 2.0 * (xyz - 1.0 / 3.0) * (xyz - 2.0 / 3.0) -
                               9.0 / 2.0 * xyz * (xyz - 2.0 / 3.0) - 9.0 / 2.0 * xyz * (xyz - 1.0 / 3.0);
            values[0].col(1) = 9.0 / 2.0 * x * (3.0 * x - 2.0) + 1.0;
            values[0].col(4) = 27.0 / 2.0 * xyz * (xyz - 1.0 / 3.0) - 27.0 / 2.0 * x * (xyz - 1.0 / 3.0) - 27.0 / 2.0 * x * xyz;
            values[0].col(5) = 27.0 / 2.0 * xyz * (x - 1.0 / 3.0) + 27.0 / 2.0 * x * xyz - 27.0 / 2.0 * x * (x - 1.0 / 3.0);
            values[0].col(6) = 27.0 / 2.0 * y * (x - 1.0 / 3.0) + 27.0 / 2.0 * y * x;
            values[0].col(7) = 27.0 / 2.0 * y * (y - 1.0 / 3.0);
            values[0].col(8) = -27.0 / 2.0 * y * (y - 1.0 / 3.0);
            values[0].col(9) = -27.0 / 2.0 * y * (xyz - 1.0 / 3.0) - 27.0 / 2.0 * y * xyz;
            values[0].col(10) = -27.0 / 2.0 * z * (xyz - 1.0 / 3.0) - 27.0 / 2.0 * z * xyz;
            values[0].col(11) = -27.0 / 2.0 * z * (z - 1.0 / 3.0);
            values[0].col(12) = 27.0 / 2.0 * z * (z - 1.0 / 3.0);
            values[0].col(13) = 27.0 / 2.0 * z * (x - 1.0 / 3.0) + 27.0 / 2.0 * z * x;
            values[0].col(16) = 27.0 * y * xyz - 27.0 * x * y;
            values[0].col(17) = 27.0 * z * xyz - 27.0 * x * z;
            values[0].col(18) = -27.0 * y * z;
            values[0].col(19) = 27.0 * y * z;

            values[1].col(0) = -9.0 / 2.0 * (xyz - 1.0 / 3.0) * (xyz - 2.0 / 3.0) -
                               9.0 / 2.0 * xyz * (xyz - 2.0 / 3.0) - 9.0 / 2.0 * xyz * (xyz - 1.0 / 3.0);
            values[1].col(2) = 9.0 / 2.0 * y * (3.0 * y - 2.0) + 1.0;
            values[1].col(4) = -27.0 / 2.0 * x * (xyz - 1.0 / 3.0) - 27.0 / 2.0 * x * xyz;
            values[1].col(5) = -27.0 / 2.0 * x * (x - 1.0 / 3.0);
            values[1].col(6) = 27.0 / 2.0 * x * (x - 1.0 / 3.0);
            values[1].col(7) = 27.0 / 2.0 * x * (y - 1.0 / 3.0) + 27.0 / 2.0 * x * y;
            values[1].col(8) = 27.0 / 2.0 * (y - 1.0 / 3.0) * xyz + 27.0 / 2.0 * y * xyz - 27.0 / 2.0 * y * (y - 1.0 / 3.0);
            values[1].col(9) = 27.0 / 2.0 * xyz * (xyz - 1.0 / 3.0) - 27.0 / 2.0 * y * (xyz - 1.0 / 3.0) - 27.0 / 2.0 * y * xyz;
            values[1].col(10) = -27.0 / 2.0 * z * (xyz - 1.0 / 3.0) - 27.0 / 2.0 * z * xyz;
            values[1].col(11) = -27.0 / 2.0 * z * (z - 1.0 / 3.0);
            values[1].col(14) = 27.0 / 2.0 * z * (z - 1.0 / 3.0);
            values[1].col(15) = 27.0 / 2.0 * z * (y - 1.0 / 3.0) + 27.0 / 2.0 * z * y;
            values[1].col(16) = 27.0 * x * xyz - 27.0 * x * y;
            values[1].col(17) = -27.0 * x * z;
            values[1].col(18) = 27.0 * z * xyz - 27.0 * y * z;
            values[1].col(19) = 27.0 * x * z;

            values[2].col(0) = -9.0 / 2.0 * (xyz - 1.0 / 3.0) * (xyz - 2.0 / 3.0) -
                               9.0 / 2.0 * xyz * (xyz - 2.0 / 3.0) - 9.0 / 2.0 * xyz * (xyz - 1.0 / 3.0);
            values[2].col(3) = 9.0 / 2.0 * z * (3.0 * z - 2.0) + 1.0;
            values[2].col(4) = -27.0 / 2.0 * x * (xyz - 1.0 / 3.0) - 27.0 / 2.0 * x * xyz;
            values[2].col(5) = -27.0 / 2.0 * x * (x - 1.0 / 3.0);
            values[2].col(8) = -27.0 / 2.0 * y * (y - 1.0 / 3.0);
            values[2].col(9) = -27.0 / 2.0 * y * (xyz - 1.0 / 3.0) - 27.0 / 2.0 * y * xyz;
            values[2].col(10) = 27.0 / 2.0 * (xyz - 1.0 / 3.0) * xyz - 27.0 / 2.0 * z * xyz - 27.0 / 2.0 * z * (xyz - 1.0 / 3.0);
            values[2].col(11) = 27.0 / 2.0 * xyz * (z - 1.0 / 3.0) + 27.0 / 2.0 * z * xyz - 27.0 / 2.0 * z * (z - 1.0 / 3.0);
            values[2].col(12) = 27.0 / 2.0 * x * (z - 1.0 / 3.0) + 27.0 / 2.0 * x * z;
            values[2].col(13) = 27.0 / 2.0 * x * (x - 1.0 / 3.0);
            values[2].col(14) = 27.0 / 2.0 * y * (z - 1.0 / 3.0) + 27.0 / 2.0 * y * z;
            values[2].col(15) = 27.0 / 2.0 * y * (y - 1.0 / 3.0);
            values[2].col(16) = -27.0 * x * y;
            values[2].col(17) = 27.0 * x * xyz - 27.0 * x * z;
            values[2].col(18) = 27.0 * y * xyz - 27.0 * y * z;
            values[2].col(19) = 27.0 * x * y;
        }
        break;
        default:
            throw std::runtime_error("order " + std::to_string(reference_element_data.Order) + "not supported yet");
        }

        return values;
    }
};
} // namespace PCC
} // namespace FEM
} // namespace Polydim

#endif
