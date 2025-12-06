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

#ifndef __FEM_Triangle_RT_MCC_2D_ReferenceElement_HPP
#define __FEM_Triangle_RT_MCC_2D_ReferenceElement_HPP

#include "Eigen/Eigen"
#include "Monomials_1D.hpp"
#include "Monomials_2D.hpp"
#include "QuadratureData.hpp"
#include "VEM_Quadrature_2D.hpp"

#include <iostream>

namespace Polydim
{
namespace FEM
{
namespace MCC
{

struct FEM_Triangle_RT_MCC_2D_Pressure_ReferenceElement_Data final
{
    unsigned int NumDofs0D;
    unsigned int NumDofs1D;
    unsigned int NumDofs2D;

    unsigned int NumBasisFunctions;

    Eigen::MatrixXd ReferenceBasisFunctionValues;
};

struct FEM_Triangle_RT_MCC_2D_Velocity_ReferenceElement_Data final
{
    struct BasisFunctions
    {
        Eigen::MatrixXd MonomialsCoefficients;

        std::vector<Eigen::MatrixXd> ReferenceBasisFunctionValues;
        Eigen::MatrixXd ReferenceBasisFunctionDivergenceValues;
    };
    unsigned int NumDofs0D;
    unsigned int NumDofs1D;
    unsigned int NumDofs2D;

    unsigned int NumBasisFunctions;

    std::map<std::array<bool, 3>, Polydim::FEM::MCC::FEM_Triangle_RT_MCC_2D_Velocity_ReferenceElement_Data::BasisFunctions> basis_functions;
};

struct FEM_Triangle_RT_MCC_2D_ReferenceElement_Data final
{
    unsigned int Dimension;
    unsigned int Order;
    unsigned int Nk;
    unsigned int Nkm1;

    Polydim::VEM::Quadrature::VEM_QuadratureData_2D Quadrature;
    std::map<std::array<bool, 3>, Polydim::VEM::Quadrature::VEM_Quadrature_2D::Edges_QuadratureData> BoundaryQuadrature;

    Polydim::Utilities::Monomials_Data monomials_2D_data;
    Eigen::Vector3d monomials_2D_center;
    double monomials_2D_scale;

    Polydim::Utilities::Monomials_Data monomials_1D_data;
    Eigen::Vector3d monomials_1D_center;
    double monomials_1D_scale;

    Eigen::MatrixXd VanderBoundary1D;

    Eigen::Matrix3d TriangleVertices = (Eigen::Matrix3d() << 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0).finished();
    Eigen::Vector3d EdgeLengths = (Eigen::Vector3d() << 1.0, sqrt(2.0), 1.0).finished();
    Eigen::Matrix3d EdgeTangents = (Eigen::Matrix3d() << 1.0, -1.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 0.0).finished();
    Eigen::Matrix3d EdgeNormals =
        (Eigen::Matrix3d() << 0.0, 1.0 / sqrt(2.0), -1.0, -1.0, 1.0 / sqrt(2.0), 0.0, 0.0, 0.0, 0.0).finished();

    FEM_Triangle_RT_MCC_2D_Velocity_ReferenceElement_Data reference_element_data_velocity;
    FEM_Triangle_RT_MCC_2D_Pressure_ReferenceElement_Data reference_element_data_pressure;
};

class FEM_Triangle_RT_MCC_2D_ReferenceElement final
{
    Polydim::Utilities::Monomials_1D monomials_1D;
    Polydim::Utilities::Monomials_2D monomials_2D;
    Polydim::VEM::Quadrature::VEM_Quadrature_2D quadrature;

    Eigen::MatrixXd ComputeVelocityMonomialsCoefficient(const std::array<bool, 3> &edge_directions,
                                                        const Polydim::FEM::MCC::FEM_Triangle_RT_MCC_2D_ReferenceElement_Data &reference_element_data) const
    {

        const auto &BoundaryQuadrature = reference_element_data.BoundaryQuadrature.at(edge_directions);

        Eigen::MatrixXd Bmatrix =
            Eigen::MatrixXd::Identity(reference_element_data.reference_element_data_velocity.NumBasisFunctions,
                                      reference_element_data.reference_element_data_velocity.NumBasisFunctions);

        for (unsigned int e = 0; e < 3; e++)
        {
            const double direction = edge_directions[e] ? 1.0 : -1.0;
            Bmatrix.block(e * reference_element_data.reference_element_data_velocity.NumDofs1D,
                          e * reference_element_data.reference_element_data_velocity.NumDofs1D,
                          reference_element_data.reference_element_data_velocity.NumDofs1D,
                          reference_element_data.reference_element_data_velocity.NumDofs1D) *= direction;
        }

        Eigen::MatrixXd Gmatrix =
            Eigen::MatrixXd::Zero(reference_element_data.reference_element_data_velocity.NumBasisFunctions,
                                  reference_element_data.reference_element_data_velocity.NumBasisFunctions);

        const Eigen::MatrixXd VanderBoundary2D = monomials_2D.Vander(reference_element_data.monomials_2D_data,
                                                                     BoundaryQuadrature.Quadrature.Points,
                                                                     reference_element_data.monomials_2D_center,
                                                                     reference_element_data.monomials_2D_scale);

        const Eigen::MatrixXd xy_boundary =
            (1.0 / reference_element_data.monomials_2D_scale) *
            (BoundaryQuadrature.Quadrature.Points.colwise() - reference_element_data.monomials_2D_center);

        // 1 edge
        Gmatrix.block(0, reference_element_data.Nk, reference_element_data.Order + 1, reference_element_data.Nk) =
            reference_element_data.VanderBoundary1D.transpose() *
            BoundaryQuadrature.WeightsTimesNormal[1].segment(0, reference_element_data.Order + 1).asDiagonal() *
            VanderBoundary2D.topRows(reference_element_data.Order + 1);

        Gmatrix.block(0, 2 * reference_element_data.Nk, reference_element_data.Order + 1, reference_element_data.Order + 1) =
            reference_element_data.VanderBoundary1D.transpose() *
            (BoundaryQuadrature.WeightsTimesNormal[1].segment(0, reference_element_data.Order + 1).array() *
             xy_boundary.row(1).segment(0, reference_element_data.Order + 1).transpose().array())
                .matrix()
                .asDiagonal() *
            VanderBoundary2D.topRightCorner(reference_element_data.Order + 1, reference_element_data.Order + 1);

        // 2 edge - obliquo

        Gmatrix.block(1 * (reference_element_data.Order + 1),
                      0,
                      reference_element_data.Order + 1,
                      reference_element_data.Nk) =
            reference_element_data.VanderBoundary1D.transpose() *
            BoundaryQuadrature.WeightsTimesNormal[0]
                .segment(1 * (reference_element_data.Order + 1), reference_element_data.Order + 1)
                .asDiagonal() *
            VanderBoundary2D.middleRows(reference_element_data.Order + 1, reference_element_data.Order + 1);

        Gmatrix.block(1 * (reference_element_data.Order + 1),
                      reference_element_data.Nk,
                      reference_element_data.Order + 1,
                      reference_element_data.Nk) =
            reference_element_data.VanderBoundary1D.transpose() *
            BoundaryQuadrature.WeightsTimesNormal[1]
                .segment(1 * (reference_element_data.Order + 1), reference_element_data.Order + 1)
                .asDiagonal() *
            VanderBoundary2D.middleRows(reference_element_data.Order + 1, reference_element_data.Order + 1);

        Gmatrix.block(1 * (reference_element_data.Order + 1),
                      2 * reference_element_data.Nk,
                      reference_element_data.Order + 1,
                      reference_element_data.Order + 1) =
            reference_element_data.VanderBoundary1D.transpose() *
            (BoundaryQuadrature.WeightsTimesNormal[0]
                     .segment(1 * (reference_element_data.Order + 1), reference_element_data.Order + 1)
                     .array() *
                 xy_boundary.row(0)
                     .segment(1 * (reference_element_data.Order + 1), reference_element_data.Order + 1)
                     .transpose()
                     .array() +
             BoundaryQuadrature.WeightsTimesNormal[1]
                     .segment(1 * (reference_element_data.Order + 1), reference_element_data.Order + 1)
                     .array() *
                 xy_boundary.row(1)
                     .segment(1 * (reference_element_data.Order + 1), reference_element_data.Order + 1)
                     .transpose()
                     .array())
                .matrix()
                .asDiagonal() *
            VanderBoundary2D.block(reference_element_data.Order + 1,
                                   reference_element_data.Nk - (reference_element_data.Order + 1),
                                   reference_element_data.Order + 1,
                                   reference_element_data.Order + 1);

        // 3  edge
        Gmatrix.block(2 * (reference_element_data.reference_element_data_velocity.NumDofs1D),
                      0,
                      reference_element_data.reference_element_data_velocity.NumDofs1D,
                      reference_element_data.Nk) =
            reference_element_data.VanderBoundary1D.transpose() *
            BoundaryQuadrature.WeightsTimesNormal[0]
                .segment(2 * (reference_element_data.Order + 1), reference_element_data.Order + 1)
                .asDiagonal() *
            VanderBoundary2D.bottomRows(reference_element_data.Order + 1);

        Gmatrix.block(2 * (reference_element_data.reference_element_data_velocity.NumDofs1D),
                      2 * reference_element_data.Nk,
                      reference_element_data.reference_element_data_velocity.NumDofs1D,
                      reference_element_data.Order + 1) =
            reference_element_data.VanderBoundary1D.transpose() *
            (BoundaryQuadrature.WeightsTimesNormal[0]
                 .segment(2 * (reference_element_data.Order + 1), reference_element_data.Order + 1)
                 .array() *
             xy_boundary.row(0)
                 .segment(2 * (reference_element_data.Order + 1), reference_element_data.Order + 1)
                 .transpose()
                 .array())
                .matrix()
                .asDiagonal() *
            VanderBoundary2D.bottomRightCorner(reference_element_data.Order + 1, reference_element_data.Order + 1);

        // internal dofs
        if (reference_element_data.Order > 0)
        {
            Gmatrix.block(3 * (reference_element_data.reference_element_data_velocity.NumDofs1D),
                          0,
                          reference_element_data.Nkm1,
                          reference_element_data.Nk) =
                reference_element_data.reference_element_data_pressure.ReferenceBasisFunctionValues
                    .leftCols(reference_element_data.Nkm1)
                    .transpose() *
                reference_element_data.Quadrature.ReferenceTriangleQuadrature.Weights.asDiagonal() *
                reference_element_data.reference_element_data_pressure.ReferenceBasisFunctionValues;

            Gmatrix.block(3 * (reference_element_data.reference_element_data_velocity.NumDofs1D) +
                              reference_element_data.Nkm1,
                          reference_element_data.Nk,
                          reference_element_data.Nkm1,
                          reference_element_data.Nk) =
                reference_element_data.reference_element_data_pressure.ReferenceBasisFunctionValues
                    .leftCols(reference_element_data.Nkm1)
                    .transpose() *
                reference_element_data.Quadrature.ReferenceTriangleQuadrature.Weights.asDiagonal() *
                reference_element_data.reference_element_data_pressure.ReferenceBasisFunctionValues;

            Gmatrix.block(3 * (reference_element_data.reference_element_data_velocity.NumDofs1D),
                          2 * reference_element_data.Nk,
                          reference_element_data.Nkm1,
                          reference_element_data.Order + 1) =
                reference_element_data.reference_element_data_pressure.ReferenceBasisFunctionValues
                    .leftCols(reference_element_data.Nkm1)
                    .transpose() *
                reference_element_data.Quadrature.ReferenceTriangleQuadrature.Weights
                    .cwiseProduct(reference_element_data.reference_element_data_pressure.ReferenceBasisFunctionValues.col(1))
                    .asDiagonal() *
                reference_element_data.reference_element_data_pressure.ReferenceBasisFunctionValues.rightCols(
                    reference_element_data.Order + 1);

            Gmatrix.block(3 * (reference_element_data.reference_element_data_velocity.NumDofs1D) +
                              reference_element_data.Nkm1,
                          2 * reference_element_data.Nk,
                          reference_element_data.Nkm1,
                          reference_element_data.Order + 1) =
                reference_element_data.reference_element_data_pressure.ReferenceBasisFunctionValues
                    .leftCols(reference_element_data.Nkm1)
                    .transpose() *
                reference_element_data.Quadrature.ReferenceTriangleQuadrature.Weights
                    .cwiseProduct(reference_element_data.reference_element_data_pressure.ReferenceBasisFunctionValues.col(2))
                    .asDiagonal() *
                reference_element_data.reference_element_data_pressure.ReferenceBasisFunctionValues.rightCols(
                    reference_element_data.Order + 1);
        }

        return Gmatrix.partialPivLu().solve(Bmatrix);
    }

    void ComputeVelocityMonomialsCoefficient(Polydim::FEM::MCC::FEM_Triangle_RT_MCC_2D_ReferenceElement_Data &reference_element_data) const
    {

        std::vector<std::array<bool, 3>> edge_directions = {{true, true, true},
                                                            {true, true, false},
                                                            {false, true, true},
                                                            {true, false, true},
                                                            {false, false, true},
                                                            {false, true, false},
                                                            {true, false, false},
                                                            {false, false, false}};

        for (unsigned int e = 0; e < edge_directions.size(); e++)
        {
            std::vector<bool> edge_dir(3);
            for (unsigned int d = 0; d < 3; d++)
                edge_dir[d] = edge_directions[e][d];

            reference_element_data.BoundaryQuadrature.insert(
                {edge_directions[e],
                 quadrature.PolygonEdgesQuadrature(reference_element_data.Quadrature.ReferenceSegmentQuadrature,
                                                   reference_element_data.TriangleVertices,
                                                   reference_element_data.EdgeLengths,
                                                   edge_dir,
                                                   reference_element_data.EdgeTangents,
                                                   reference_element_data.EdgeNormals)});

            FEM_Triangle_RT_MCC_2D_Velocity_ReferenceElement_Data::BasisFunctions basis_functions;
            basis_functions.MonomialsCoefficients = ComputeVelocityMonomialsCoefficient(edge_directions[e], reference_element_data);

            basis_functions.ReferenceBasisFunctionValues =
                EvaluateVelocityBasisFunctions(reference_element_data.Quadrature.ReferenceTriangleQuadrature.Points,
                                               basis_functions.MonomialsCoefficients,
                                               reference_element_data);

            basis_functions.ReferenceBasisFunctionDivergenceValues =
                EvaluateVelocityBasisFunctionsDivergence(reference_element_data.Quadrature.ReferenceTriangleQuadrature.Points,
                                                         basis_functions.MonomialsCoefficients,
                                                         reference_element_data);

            reference_element_data.reference_element_data_velocity.basis_functions.insert({edge_directions[e], basis_functions});
        }
    }

  public:
    FEM_Triangle_RT_MCC_2D_ReferenceElement_Data Create(const unsigned int order) const
    {
        Polydim::FEM::MCC::FEM_Triangle_RT_MCC_2D_ReferenceElement_Data result;

        result.Dimension = 2;
        result.Order = order;
        result.Nk = (order + 1) * (order + 2) / 2;
        result.Nkm1 = order * (order + 1) / 2;

        result.reference_element_data_pressure.NumDofs0D = 0;
        result.reference_element_data_pressure.NumDofs1D = 0;
        result.reference_element_data_pressure.NumDofs2D = (order + 1) * (order + 2) / 2;

        result.reference_element_data_pressure.NumBasisFunctions = result.reference_element_data_pressure.NumDofs2D;

        result.reference_element_data_velocity.NumDofs0D = 0;
        result.reference_element_data_velocity.NumDofs1D = order + 1;
        result.reference_element_data_velocity.NumDofs2D = result.Dimension * (order) * (order + 1) / 2;

        result.reference_element_data_velocity.NumBasisFunctions =
            (result.Dimension + 1) * result.reference_element_data_velocity.NumDofs1D +
            result.reference_element_data_velocity.NumDofs2D;

        result.Quadrature = quadrature.Compute_MCC_2D(order);

        result.monomials_2D_data = monomials_2D.Compute(order);
        result.monomials_2D_scale = 1.0;
        result.monomials_2D_center << 1.0 / 3.0, 1.0 / 3.0, 0.0;

        result.monomials_1D_data = monomials_1D.Compute(order);
        result.monomials_1D_scale = 1.0;
        result.monomials_1D_center << 0.5, 0.0, 0.0;

        result.reference_element_data_pressure.ReferenceBasisFunctionValues =
            monomials_2D.Vander(result.monomials_2D_data,
                                result.Quadrature.ReferenceTriangleQuadrature.Points,
                                result.monomials_2D_center,
                                result.monomials_2D_scale);

        result.VanderBoundary1D = monomials_1D.Vander(result.monomials_1D_data,
                                                      result.Quadrature.ReferenceSegmentQuadrature.Points,
                                                      result.monomials_1D_center,
                                                      result.monomials_1D_scale);

        ComputeVelocityMonomialsCoefficient(result);

        return result;
    }
    // ***************************************************************************
    std::vector<Eigen::MatrixXd> EvaluateVelocityBasisFunctions(const Eigen::MatrixXd &points,
                                                                const Eigen::MatrixXd &MonomialsCoefficients,
                                                                const Polydim::FEM::MCC::FEM_Triangle_RT_MCC_2D_ReferenceElement_Data &reference_element_data) const
    {

        const Eigen::MatrixXd xy_internal = (1.0 / reference_element_data.monomials_2D_scale) *
                                            (points.colwise() - reference_element_data.monomials_2D_center);

        std::vector<Eigen::MatrixXd> BasisFunctionValues(
            2,
            Eigen::MatrixXd::Zero(points.cols(), reference_element_data.reference_element_data_velocity.NumBasisFunctions));

        const auto Vander = monomials_2D.Vander(reference_element_data.monomials_2D_data,
                                                points,
                                                reference_element_data.monomials_2D_center,
                                                reference_element_data.monomials_2D_scale);

        BasisFunctionValues[0] = Vander * MonomialsCoefficients.topRows(reference_element_data.Nk) +
                                 xy_internal.row(0).asDiagonal() * Vander.rightCols(reference_element_data.Order + 1) *
                                     MonomialsCoefficients.bottomRows(reference_element_data.Order + 1);

        BasisFunctionValues[1] =
            Vander * MonomialsCoefficients.middleRows(reference_element_data.Nk, reference_element_data.Nk) +
            xy_internal.row(1).asDiagonal() * Vander.rightCols(reference_element_data.Order + 1) *
                MonomialsCoefficients.bottomRows(reference_element_data.Order + 1);

        return BasisFunctionValues;
    }
    // ***************************************************************************
    Eigen::MatrixXd EvaluatePressureBasisFunctions(const Eigen::MatrixXd &points,
                                                   const Polydim::FEM::MCC::FEM_Triangle_RT_MCC_2D_ReferenceElement_Data &reference_element_data) const
    {
        return monomials_2D.Vander(reference_element_data.monomials_2D_data,
                                   points,
                                   reference_element_data.monomials_2D_center,
                                   reference_element_data.monomials_2D_scale);
    }
    // ***************************************************************************
    Eigen::MatrixXd EvaluateVelocityBasisFunctionsDivergence(const Eigen::MatrixXd &points,
                                                             const Eigen::MatrixXd &MonomialsCoefficients,
                                                             const Polydim::FEM::MCC::FEM_Triangle_RT_MCC_2D_ReferenceElement_Data &reference_element_data) const
    {

        const Eigen::MatrixXd xy_internal = (1.0 / reference_element_data.monomials_2D_scale) *
                                            (points.colwise() - reference_element_data.monomials_2D_center);

        const auto Vander = monomials_2D.Vander(reference_element_data.monomials_2D_data,
                                                points,
                                                reference_element_data.monomials_2D_center,
                                                reference_element_data.monomials_2D_scale);

        const auto VanderDerivatives =
            monomials_2D.VanderDerivatives(reference_element_data.monomials_2D_data, Vander, reference_element_data.monomials_2D_scale);

        Eigen::MatrixXd divergence_values =
            VanderDerivatives[0] * MonomialsCoefficients.topRows(reference_element_data.Nk) +
            VanderDerivatives[1] * MonomialsCoefficients.middleRows(reference_element_data.Nk, reference_element_data.Nk) +
            2.0 * (1.0 / reference_element_data.monomials_2D_scale) * Vander.rightCols(reference_element_data.Order + 1) *
                MonomialsCoefficients.bottomRows(reference_element_data.Order + 1) +
            xy_internal.row(0).asDiagonal() * VanderDerivatives[0].rightCols(reference_element_data.Order + 1) *
                MonomialsCoefficients.bottomRows(reference_element_data.Order + 1) +
            xy_internal.row(1).asDiagonal() * VanderDerivatives[1].rightCols(reference_element_data.Order + 1) *
                MonomialsCoefficients.bottomRows(reference_element_data.Order + 1);

        return divergence_values;
    }
    // ***************************************************************************
};
} // namespace MCC
} // namespace FEM
} // namespace Polydim

#endif
