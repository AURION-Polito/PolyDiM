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

#ifndef __TEST_FEM_Triangle_RT_MCC_2D_LocalSpace_H
#define __TEST_FEM_Triangle_RT_MCC_2D_LocalSpace_H

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "Assembler_Utilities.hpp"
#include "DOFsManager.hpp"
#include "Eigen_Array.hpp"
#include "EllipticEquation.hpp"
#include "FEM_MCC_2D_LocalSpace.hpp"
#include "FEM_Triangle_RT_MCC_2D_LocalSpace.hpp"
#include "FEM_Triangle_RT_MCC_2D_ReferenceElement.hpp"
#include "GeometryUtilities.hpp"
#include "IMeshDAO.hpp"
#include "LocalSpace_MCC_2D.hpp"
#include "MeshMatricesDAO_mesh_connectivity_data.hpp"
#include "PDE_Mesh_Utilities.hpp"
#include "Quadrature_Gauss2D_Triangle.hpp"

namespace Polydim
{
namespace UnitTesting
{
TEST(Test_FEM_Triangle_RT_MCC_2D, Test_FEM_Triangle_RT_MCC_2D_Reference_Element)
{
    const Polydim::FEM::MCC::FEM_Triangle_RT_MCC_2D_ReferenceElement reference_element;

    const auto referenceQuadrature = Gedim::Quadrature::Quadrature_Gauss2D_Triangle::FillPointsAndWeights(21);
    const Eigen::MatrixXd &referenceQuadraturePoints = referenceQuadrature.Points;

    Polydim::Utilities::Monomials_2D monomials_2D;
    Polydim::Utilities::Monomials_1D monomials_1D;

    for (unsigned int o = 0; o < 5; o++)
    {
        const auto reference_element_data = reference_element.Create(o);

        const Eigen::MatrixXd pressureBasisValues =
            reference_element.EvaluatePressureBasisFunctions(referenceQuadraturePoints, reference_element_data);
        const std::vector<Eigen::MatrixXd> velocityBasisValues =
            reference_element.EvaluateVelocityBasisFunctions(referenceQuadraturePoints, reference_element_data);
        const Eigen::MatrixXd diverVelocityBasisFunctions =
            reference_element.EvaluateVelocityBasisFunctionsDivergence(referenceQuadraturePoints, reference_element_data);

        // Compute dofs of monomials
        Eigen::MatrixXd vander_boundary = monomials_2D.Vander(reference_element_data.monomials_2D_data,
                                                              reference_element_data.BoundaryQuadrature.Quadrature.Points,
                                                              Eigen::Vector3d::Zero(),
                                                              1.0);

        Eigen::MatrixXd vander_internal =
            monomials_2D.Vander(reference_element_data.monomials_2D_data,
                                reference_element_data.Quadrature.ReferenceTriangleQuadrature.Points,
                                Eigen::Vector3d::Zero(),
                                1.0);

        Eigen::MatrixXd Dmatrix = Eigen::MatrixXd::Zero(reference_element_data.reference_element_data_velocity.NumBasisFunctions,
                                                        2 * reference_element_data.Nk);

        for (unsigned int e = 0; e < 3; e++)
        {
            Dmatrix.block(e * (reference_element_data.Order + 1),
                          0,
                          (reference_element_data.Order + 1),
                          reference_element_data.Nk) =
                reference_element_data.VanderBoundary1D.transpose() *
                reference_element_data.BoundaryQuadrature.WeightsTimesNormal[0]
                    .segment(e * (reference_element_data.Order + 1), (reference_element_data.Order + 1))
                    .asDiagonal() *
                vander_boundary.middleRows(e * (reference_element_data.Order + 1), (reference_element_data.Order + 1));

            Dmatrix.block(e * (reference_element_data.Order + 1),
                          reference_element_data.Nk,
                          (reference_element_data.Order + 1),
                          reference_element_data.Nk) =
                reference_element_data.VanderBoundary1D.transpose() *
                reference_element_data.BoundaryQuadrature.WeightsTimesNormal[1]
                    .segment(e * (reference_element_data.Order + 1), (reference_element_data.Order + 1))
                    .asDiagonal() *
                vander_boundary.middleRows(e * (reference_element_data.Order + 1), (reference_element_data.Order + 1));
        }

        if (reference_element_data.Order > 0)
        {
            Dmatrix.block(3 * (reference_element_data.Order + 1), 0, reference_element_data.Nkm1, reference_element_data.Nk) =
                reference_element_data.reference_element_data_pressure.ReferenceBasisFunctionValues
                    .leftCols(reference_element_data.Nkm1)
                    .transpose() *
                reference_element_data.Quadrature.ReferenceTriangleQuadrature.Weights.asDiagonal() * vander_internal;

            Dmatrix.block(3 * (reference_element_data.Order + 1) + reference_element_data.Nkm1,
                          reference_element_data.Nk,
                          reference_element_data.Nkm1,
                          reference_element_data.Nk) =
                reference_element_data.reference_element_data_pressure.ReferenceBasisFunctionValues
                    .leftCols(reference_element_data.Nkm1)
                    .transpose() *
                reference_element_data.Quadrature.ReferenceTriangleQuadrature.Weights.asDiagonal() * vander_internal;
        }

        Eigen::MatrixXd monomials_values =
            monomials_2D.Vander(reference_element_data.monomials_2D_data, referenceQuadrature.Points, Eigen::Vector3d::Zero(), 1.0);

        const Eigen::MatrixXd ConsistencyErrorX =
            monomials_values - velocityBasisValues[0] * Dmatrix.leftCols(reference_element_data.Nk);
        const Eigen::MatrixXd ConsistencyErrorY =
            monomials_values - velocityBasisValues[1] * Dmatrix.rightCols(reference_element_data.Nk);

        const double max_error_x = ConsistencyErrorX.cwiseAbs().maxCoeff();
        const double max_error_y = ConsistencyErrorY.cwiseAbs().maxCoeff();

        ASSERT_TRUE(max_error_x < 1.0e-10);
        ASSERT_TRUE(max_error_y < 1.0e-10);
    }
}

TEST(Test_FEM_Triangle_RT_MCC_2D, Test_FEM_Triangle_RT_MCC_2D_Reference_Element_Patch)
{
    const Polydim::FEM::MCC::FEM_Triangle_RT_MCC_2D_ReferenceElement reference_element;

    const auto referenceQuadrature = Gedim::Quadrature::Quadrature_Gauss2D_Triangle::FillPointsAndWeights(4);
    const Eigen::MatrixXd &referenceQuadraturePoints = referenceQuadrature.Points;

    Polydim::Utilities::Monomials_2D monomials_2D;
    Polydim::Utilities::Monomials_1D monomials_1D;

    for (unsigned int o = 1; o < 2; o++)
    {
        const auto reference_element_data = reference_element.Create(o);

        const Eigen::MatrixXd pressureBasisValues =
            reference_element.EvaluatePressureBasisFunctions(referenceQuadraturePoints, reference_element_data);
        const std::vector<Eigen::MatrixXd> velocityBasisValues =
            reference_element.EvaluateVelocityBasisFunctions(referenceQuadraturePoints, reference_element_data);
        const Eigen::MatrixXd diverVelocityBasisFunctions =
            reference_element.EvaluateVelocityBasisFunctionsDivergence(referenceQuadraturePoints, reference_element_data);

        // Compute dofs of monomials
        Eigen::MatrixXd vander_boundary = monomials_2D.Vander(reference_element_data.monomials_2D_data,
                                                              reference_element_data.BoundaryQuadrature.Quadrature.Points,
                                                              Eigen::Vector3d::Zero(),
                                                              1.0);

        const auto vander_boundary_der = monomials_2D.VanderDerivatives(reference_element_data.monomials_2D_data,
                                                                        vander_boundary,
                                                                        1.0);

        Eigen::MatrixXd vander_internal =
            monomials_2D.Vander(reference_element_data.monomials_2D_data,
                                reference_element_data.Quadrature.ReferenceTriangleQuadrature.Points,
                                Eigen::Vector3d::Zero(),
                                1.0);

        const auto vander_internal_der = monomials_2D.VanderDerivatives(reference_element_data.monomials_2D_data,
                                                                        vander_internal,
                                                                        1.0);

        Eigen::VectorXd pressure_boundary = vander_boundary.col(1) + vander_boundary.col(2) + 0.5 * vander_boundary.col(0);
        std::vector<Eigen::VectorXd> velocity_boundary = {-(vander_boundary_der[0].col(1) + vander_boundary_der[0].col(2)),
                                                          -(vander_boundary_der[1].col(1) + vander_boundary_der[1].col(2))};

        Eigen::VectorXd pressure_internal = vander_internal.col(1) + vander_internal.col(2) + 0.5 * vander_internal.col(0);
        std::vector<Eigen::VectorXd> velocity_internal = {-(vander_internal_der[0].col(1) + vander_internal_der[0].col(2)),
                                                          -(vander_internal_der[1].col(1) + vander_internal_der[1].col(2))};

        Eigen::MatrixXd Dmatrix = Eigen::MatrixXd::Zero(reference_element_data.reference_element_data_velocity.NumBasisFunctions,
                                                        1);

        for (unsigned int e = 0; e < 3; e++)
        {
            Dmatrix.block(e * (reference_element_data.Order + 1),
                          0,
                          (reference_element_data.Order + 1),
                          1) =
                reference_element_data.VanderBoundary1D.transpose() *
                    reference_element_data.BoundaryQuadrature.WeightsTimesNormal[0]
                        .segment(e * (reference_element_data.Order + 1), (reference_element_data.Order + 1))
                        .asDiagonal() *
                    velocity_boundary[0].middleRows(e * (reference_element_data.Order + 1), (reference_element_data.Order + 1)) +
                reference_element_data.VanderBoundary1D.transpose() *
                    reference_element_data.BoundaryQuadrature.WeightsTimesNormal[1]
                        .segment(e * (reference_element_data.Order + 1), (reference_element_data.Order + 1))
                        .asDiagonal() *
                    velocity_boundary[1].middleRows(e * (reference_element_data.Order + 1), (reference_element_data.Order + 1));

        }

        if (reference_element_data.Order > 0)
        {
            Dmatrix.block(3 * (reference_element_data.Order + 1), 0, reference_element_data.Nkm1, 1) =
                reference_element_data.reference_element_data_pressure.ReferenceBasisFunctionValues
                    .leftCols(reference_element_data.Nkm1)
                    .transpose() *
                reference_element_data.Quadrature.ReferenceTriangleQuadrature.Weights.asDiagonal() * velocity_internal[0];

            Dmatrix.block(3 * (reference_element_data.Order + 1) + reference_element_data.Nkm1, 0, reference_element_data.Nkm1, 1) =
                reference_element_data.reference_element_data_pressure.ReferenceBasisFunctionValues
                    .leftCols(reference_element_data.Nkm1)
                    .transpose() *
                reference_element_data.Quadrature.ReferenceTriangleQuadrature.Weights.asDiagonal() * velocity_internal[1];
        }

        const Eigen::VectorXd &weights = reference_element_data.Quadrature.ReferenceTriangleQuadrature.Weights;
        Eigen::MatrixXd monomials_values =
            monomials_2D.Vander(reference_element_data.monomials_2D_data, referenceQuadrature.Points, Eigen::Vector3d::Zero(), 1.0);


        const auto monomials_values_der =
            monomials_2D.VanderDerivatives(reference_element_data.monomials_2D_data, monomials_values, 1.0);


        Eigen::VectorXd pressure = monomials_values.col(1) + monomials_values.col(2) + 0.5 * monomials_values.col(0);
        std::vector<Eigen::VectorXd> velocity = {-(monomials_values_der[0].col(1) + monomials_values_der[0].col(2)),
                                                 -(monomials_values_der[1].col(1) + monomials_values_der[1].col(2))};

        Eigen::VectorXd source_term = pressure;

        const Eigen::VectorXd local_rhs_A  = velocityBasisValues[0].transpose() * weights.asDiagonal() * velocity[0]
                                            + velocityBasisValues[1].transpose() * weights.asDiagonal() * velocity[1];
        const auto global_A = velocityBasisValues[0].transpose() * weights.asDiagonal() * velocityBasisValues[0] +
                              velocityBasisValues[1].transpose() * weights.asDiagonal() * velocityBasisValues[1];

        Eigen::VectorXd Dmatrix_pressure = Eigen::VectorXd::Zero(3);
        Dmatrix_pressure << 0.5 + reference_element_data.monomials_2D_center(0) + reference_element_data.monomials_2D_center(1), 1.0, 1.0;

        Eigen::VectorXd discrete_pressure = pressureBasisValues * Dmatrix_pressure;

        // const Eigen::VectorXd velocity_error = global_A.block(6, 6, 2, 2) * Dmatrix.bottomRows(2) - local_rhs_A.bottomRows(2)  - (- global_A.block(6, 0, 2, 6) * Dmatrix.topRows(6));
        const Eigen::VectorXd value_error = diverVelocityBasisFunctions.rightCols(2).transpose() * weights.asDiagonal()  * pressureBasisValues * Dmatrix_pressure;
        const Eigen::VectorXd velocity_error = global_A.block(6, 6, 2, 2) * Dmatrix.bottomRows(2) - diverVelocityBasisFunctions.rightCols(2).transpose() * weights.asDiagonal()  * pressureBasisValues * Dmatrix_pressure
                                               - (- global_A.block(6, 0, 2, 6) * Dmatrix.topRows(6));
        const Eigen::VectorXd div_velocity_error = pressureBasisValues.transpose() * weights.asDiagonal() * diverVelocityBasisFunctions.rightCols(2) * Dmatrix.bottomRows(2)
                                                   + pressureBasisValues.transpose() * weights.asDiagonal()  * pressureBasisValues * Dmatrix_pressure
                                                   - (pressureBasisValues.transpose() * weights.asDiagonal()  * source_term
                                                      - pressureBasisValues.transpose() * weights.asDiagonal() * diverVelocityBasisFunctions.leftCols(6) * Dmatrix.topRows(6));

        Eigen::MatrixXd global_matrix = Eigen::MatrixXd::Zero(5, 5);
        Eigen::MatrixXd neumann_matrix = Eigen::MatrixXd::Zero(5, 6);
        Eigen::VectorXd global_rhs = Eigen::VectorXd::Zero(5);
        Eigen::VectorXd solution = Eigen::VectorXd::Zero(5);
        Eigen::VectorXd neuman_solution = Eigen::VectorXd::Zero(6);
        neuman_solution = Dmatrix.topRows(6);
        solution << Dmatrix.bottomRows(2), Dmatrix_pressure;
        global_rhs.bottomRows(3) = pressureBasisValues.transpose() * weights.asDiagonal()  * source_term;
        neumann_matrix.block(2, 0, 3, 6) = pressureBasisValues.transpose() * weights.asDiagonal() * diverVelocityBasisFunctions.leftCols(6);

        global_matrix.block(2, 0, 3, 2) = pressureBasisValues.transpose() * weights.asDiagonal() * diverVelocityBasisFunctions.rightCols(2);
        global_matrix.block(2, 2, 3, 3) = pressureBasisValues.transpose() * weights.asDiagonal()  * pressureBasisValues;


        const double max_error_vel = velocity_error.cwiseAbs().maxCoeff();

        const double max_error = (global_matrix * solution - (global_rhs - neumann_matrix * neuman_solution)).cwiseAbs().maxCoeff();
        const double max_error_div = div_velocity_error.cwiseAbs().maxCoeff();

        ASSERT_TRUE(max_error_vel < 1.0e-10);
        ASSERT_TRUE(max_error < 1.0e-10);
        ASSERT_TRUE(max_error_div < 1.0e-10);
    }
}

TEST(Test_FEM_Triangle_RT_MCC_2D, Test_FEM_Triangle_RT_MCC_2D_Local_Space)
{
    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = 1.0e-8;
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    const auto poligon_vertices = geometry_utilities.CreateTriangle(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                    Eigen::Vector3d(1.5, 0.0, 0.0),
                                                                    Eigen::Vector3d(0.0, 2.0, 0.0));
    const std::vector<bool> polygon_edges_direction(3, true);
    const auto polygon_edges_tangent = geometry_utilities.PolygonEdgeTangents(poligon_vertices);
    const auto polygon_edges_length = geometry_utilities.PolygonEdgeLengths(poligon_vertices);
    const auto polygon_edges_normal = geometry_utilities.PolygonEdgeNormals(poligon_vertices);

    Polydim::FEM::MCC::FEM_MCC_2D_Polygon_Geometry polygon_geometry;

    polygon_geometry = {geometry_utilities.Tolerance1D(),
                        geometry_utilities.Tolerance2D(),
                        poligon_vertices,
                        polygon_edges_length,
                        polygon_edges_direction,
                        polygon_edges_tangent,
                        polygon_edges_normal};

    const Polydim::FEM::MCC::FEM_MCC_2D_ReferenceElement reference_element;
    const Polydim::FEM::MCC::FEM_MCC_2D_LocalSpace local_space;

    Polydim::Utilities::Monomials_2D monomials_2D;
    Polydim::Utilities::Monomials_1D monomials_1D;
    Polydim::VEM::Quadrature::VEM_Quadrature_2D quadrature;

    const auto referenceQuadrature = Gedim::Quadrature::Quadrature_Gauss2D_Triangle::FillPointsAndWeights(21);
    const auto quadrature_data = quadrature.PolygonInternalQuadrature(referenceQuadrature, {poligon_vertices});

    for (unsigned int o = 0; o < 5; o++)
    {
        const auto reference_element_data = reference_element.Create(o, Polydim::FEM::MCC::FEM_MCC_Types::RT);
        const auto local_space_data = local_space.CreateLocalSpace(reference_element_data, polygon_geometry);

        const Eigen::MatrixXd pressureBasisValues =
            local_space.ComputePressureBasisFunctionsValues(reference_element_data, local_space_data, quadrature_data.Points);
        const std::vector<Eigen::MatrixXd> velocityBasisValues =
            local_space.ComputeVelocityBasisFunctionsValues(reference_element_data, local_space_data, quadrature_data.Points);
        const Eigen::MatrixXd diverVelocityBasisFunctions =
            local_space.ComputeVelocityBasisFunctionsDivergenceValues(reference_element_data,
                                                                      local_space_data,
                                                                      quadrature_data.Points);

        const auto internal_quadrature = local_space.InternalQuadrature(
            reference_element_data.rt_triangle_reference_element_data.Quadrature.ReferenceTriangleQuadrature,
            local_space_data);
        const auto boundary_quadrature = local_space.BoundaryQuadrature(
            local_space_data,
            reference_element_data.rt_triangle_reference_element_data.Quadrature.ReferenceSegmentQuadrature,
            polygon_geometry);

        // Compute dofs of monomials
        Eigen::MatrixXd vander_internal = monomials_2D.Vander(reference_element_data.rt_triangle_reference_element_data.monomials_2D_data,
                                                              internal_quadrature.Points,
                                                              Eigen::Vector3d::Zero(),
                                                              1.0);

        Eigen::MatrixXd Dmatrix = Eigen::MatrixXd::Zero(local_space_data.NumVelocityBasisFunctions,
                                                        2 * reference_element_data.rt_triangle_reference_element_data.Nk);

        for (unsigned int e = 0; e < 3; e++)
        {
            Eigen::MatrixXd vander_boundary =
                monomials_2D.Vander(reference_element_data.rt_triangle_reference_element_data.monomials_2D_data,
                                    boundary_quadrature[e].Points,
                                    Eigen::Vector3d::Zero(),
                                    1.0);

            Dmatrix.block(e * (reference_element_data.Order + 1),
                          0,
                          (reference_element_data.Order + 1),
                          reference_element_data.rt_triangle_reference_element_data.Nk) =
                reference_element_data.rt_triangle_reference_element_data.VanderBoundary1D.transpose() *
                (polygon_edges_normal(0, e) * boundary_quadrature[e].Weights).asDiagonal() * vander_boundary;

            Dmatrix.block(e * (reference_element_data.Order + 1),
                          reference_element_data.rt_triangle_reference_element_data.Nk,
                          (reference_element_data.Order + 1),
                          reference_element_data.rt_triangle_reference_element_data.Nk) =
                reference_element_data.rt_triangle_reference_element_data.VanderBoundary1D.transpose() *
                (polygon_edges_normal(1, e) * boundary_quadrature[e].Weights).asDiagonal() * vander_boundary;
        }

        if (reference_element_data.Order > 0)
        {
            Dmatrix.block(3 * (reference_element_data.Order + 1),
                          0,
                          reference_element_data.rt_triangle_reference_element_data.Nkm1,
                          reference_element_data.rt_triangle_reference_element_data.Nk) =
                reference_element_data.rt_triangle_reference_element_data.reference_element_data_pressure
                    .ReferenceBasisFunctionValues
                    .leftCols(reference_element_data.rt_triangle_reference_element_data.Nkm1)
                    .transpose() *
                reference_element_data.rt_triangle_reference_element_data.Quadrature.ReferenceTriangleQuadrature.Weights.asDiagonal() *
                local_space.MapInvVelocityValues(
                    local_space_data,
                    {vander_internal, Eigen::MatrixXd::Zero(vander_internal.rows(), vander_internal.cols())})[0];

            Dmatrix.block(3 * (reference_element_data.Order + 1) +
                              reference_element_data.rt_triangle_reference_element_data.Nkm1,
                          reference_element_data.rt_triangle_reference_element_data.Nk,
                          reference_element_data.rt_triangle_reference_element_data.Nkm1,
                          reference_element_data.rt_triangle_reference_element_data.Nk) =
                reference_element_data.rt_triangle_reference_element_data.reference_element_data_pressure
                    .ReferenceBasisFunctionValues
                    .leftCols(reference_element_data.rt_triangle_reference_element_data.Nkm1)
                    .transpose() *
                reference_element_data.rt_triangle_reference_element_data.Quadrature.ReferenceTriangleQuadrature.Weights.asDiagonal() *
                local_space.MapInvVelocityValues(
                    local_space_data,
                    {Eigen::MatrixXd::Zero(vander_internal.rows(), vander_internal.cols()), vander_internal})[1];
        }

        Eigen::MatrixXd monomials_values = monomials_2D.Vander(reference_element_data.rt_triangle_reference_element_data.monomials_2D_data,
                                                               quadrature_data.Points,
                                                               Eigen::Vector3d::Zero(),
                                                               1.0);

        const std::vector<Eigen::MatrixXd> der_mon_values =
            monomials_2D.VanderDerivatives(reference_element_data.rt_triangle_reference_element_data.monomials_2D_data,
                                           monomials_values,
                                           1.0);

        const Eigen::MatrixXd diveregence_values = der_mon_values[0] + der_mon_values[1];

        const Eigen::MatrixXd ConsistencyErrorDiv1 =
            der_mon_values[0] -
            diverVelocityBasisFunctions * Dmatrix.leftCols(reference_element_data.rt_triangle_reference_element_data.Nk);

        const Eigen::MatrixXd ConsistencyErrorDiv2 =
            der_mon_values[1] - diverVelocityBasisFunctions *
                                    Dmatrix.rightCols(reference_element_data.rt_triangle_reference_element_data.Nk);

        const Eigen::MatrixXd ConsistencyErrorX =
            monomials_values -
            velocityBasisValues[0] * Dmatrix.leftCols(reference_element_data.rt_triangle_reference_element_data.Nk);
        const Eigen::MatrixXd ConsistencyErrorY =
            monomials_values -
            velocityBasisValues[1] * Dmatrix.rightCols(reference_element_data.rt_triangle_reference_element_data.Nk);

        const double max_error_x = ConsistencyErrorX.cwiseAbs().maxCoeff();
        const double max_error_y = ConsistencyErrorY.cwiseAbs().maxCoeff();

        const double max_error_div_1 = ConsistencyErrorDiv1.cwiseAbs().maxCoeff();
        const double max_error_div_2 = ConsistencyErrorDiv2.cwiseAbs().maxCoeff();

        ASSERT_TRUE(max_error_x < 1.0e-10);
        ASSERT_TRUE(max_error_y < 1.0e-10);

        ASSERT_TRUE(max_error_div_1 < 1.0e-10);
        ASSERT_TRUE(max_error_div_2 < 1.0e-10);
    }
}


std::array<Eigen::VectorXd, 9> inverse_diffusion_term(const Eigen::MatrixXd &points)
{
    return {Eigen::VectorXd::Constant(points.cols(), 1.0),
            Eigen::VectorXd::Constant(points.cols(), 0.0),
            Eigen::VectorXd::Zero(points.cols()),
            Eigen::VectorXd::Constant(points.cols(), 0.0),
            Eigen::VectorXd::Constant(points.cols(), 1.0),
            Eigen::VectorXd::Zero(points.cols()),
            Eigen::VectorXd::Zero(points.cols()),
            Eigen::VectorXd::Zero(points.cols()),
            Eigen::VectorXd::Constant(points.cols(), 0.0)};

};

Eigen::VectorXd source_term(const Eigen::MatrixXd &points)
{
    unsigned int order = 1;
    Eigen::ArrayXd second_derivatives = Eigen::ArrayXd::Constant(points.cols(), 0.0);
    Eigen::ArrayXd solution = Eigen::ArrayXd::Constant(points.cols(), 1.0);
    const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + 0.5;

    if (order > 1)
    {
        second_derivatives = Eigen::ArrayXd::Constant(points.cols(), 1.0);
        for (int i = 0; i < order - 2; i++)
            second_derivatives = second_derivatives * polynomial;

        solution = second_derivatives * polynomial * polynomial;
        second_derivatives *= order * (order - 1);
    }
    else if (order == 1)
        solution = polynomial;

    return -2.0 * second_derivatives + solution;
};

Eigen::VectorXd reaction_term(const Eigen::MatrixXd &points)
{
    return Eigen::VectorXd::Ones(points.cols());
}

TEST(Test_FEM_Triangle_RT_MCC_2D, Test_FEM_Triangle_RT_MCC_2D_Local_Space_Patch)
{
    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = 1.0e-8;
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);
    Gedim::MeshUtilities meshUtilities;

    unsigned int o = 1;

    const Polydim::FEM::MCC::FEM_MCC_2D_ReferenceElement reference_element;
    const Polydim::FEM::MCC::FEM_MCC_2D_LocalSpace local_space;


    Polydim::Utilities::Monomials_2D monomials_2D;
    Polydim::Utilities::Monomials_1D monomials_1D;


    const auto reference_element_data =
        Polydim::PDETools::LocalSpace_MCC_2D::CreateReferenceElement(Polydim::PDETools::LocalSpace_MCC_2D::MethodTypes::FEM_RT_MCC, 1);

    Gedim::MeshMatrices meshData;
    Gedim::MeshMatricesDAO mesh(meshData);

    Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D domain;

    domain.area = 1.0;

    domain.vertices = Eigen::MatrixXd::Zero(3, 4);
    domain.vertices.row(0) << 0.0, 1.0, 1.0, 0.0;
    domain.vertices.row(1) << 0.0, 0.0, 1.0, 1.0;

    domain.shape_type = Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D::Domain_Shape_Types::Parallelogram;

    Polydim::PDETools::Mesh::PDE_Mesh_Utilities::create_mesh_2D(geometry_utilities,
                                                                meshUtilities,
                                                                Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::Triangular,
                                                                domain,
                                                                0.01,
                                                                mesh);

    const auto mesh_geometric_data = Polydim::PDETools::Mesh::PDE_Mesh_Utilities::compute_mesh_2D_geometry_data(geometry_utilities, meshUtilities, mesh);

    Polydim::PDETools::Mesh::MeshMatricesDAO_mesh_connectivity_data mesh_connectivity_data(mesh);

    const auto reference_element_num_dofs = Polydim::PDETools::LocalSpace_MCC_2D::ReferenceElementNumDOFs(reference_element_data);

    Polydim::PDETools::DOFs::DOFsManager dofManager;
    std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> meshDOFsInfo(2);
    std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> dofs_data(2);

    std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> boundary_info =
        {{0, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
         {1, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
         {2, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
         {3, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
         {4, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
         {5, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
         {6, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 2}},
         {7, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 3}},
         {8, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 4}}};


    for (unsigned int h = 0; h < 2; h++)
    {
        meshDOFsInfo[h] =
            dofManager.Create_Constant_DOFsInfo_2D(mesh_connectivity_data, {reference_element_num_dofs[h], boundary_info});

        dofs_data[h] = dofManager.CreateDOFs_2D(meshDOFsInfo[h], mesh_connectivity_data);
    }

    const auto count_dofs = Polydim::PDETools::Assembler_Utilities::count_dofs(dofs_data);

    Gedim::Eigen_Array<> global_solution;
    global_solution.SetSize(count_dofs.num_total_dofs);
    Gedim::Eigen_Array<> global_neumann_solution;
    global_neumann_solution.SetSize(count_dofs.num_total_strong);
    for(unsigned int c = 0; c < mesh.Cell2DTotalNumber(); c++)
    {

        const auto local_space_data = Polydim::PDETools::LocalSpace_MCC_2D::CreateLocalSpace(geometry_utilities.Tolerance1D(),
                                                                                             geometry_utilities.Tolerance2D(),
                                                                                             mesh_geometric_data,
                                                                                             c,
                                                                                             reference_element_data);

        const Eigen::MatrixXd pressureBasisValues =
            Polydim::PDETools::LocalSpace_MCC_2D::PressureBasisFunctionsValues(reference_element_data, local_space_data);
        const std::vector<Eigen::MatrixXd> velocityBasisValues =
            Polydim::PDETools::LocalSpace_MCC_2D::VelocityBasisFunctionsValues(reference_element_data, local_space_data);
        const Eigen::MatrixXd diverVelocityBasisFunctions =
            Polydim::PDETools::LocalSpace_MCC_2D::VelocityBasisFunctionsDivergenceValues(reference_element_data, local_space_data);

        const auto internal_quadrature = local_space.InternalQuadrature(
            reference_element_data.FEM_ReferenceElement_Data.rt_triangle_reference_element_data.Quadrature.ReferenceTriangleQuadrature,
            local_space_data.FEM_LocalSpace_Data);
        const auto boundary_quadrature = local_space.BoundaryQuadrature(
            local_space_data.FEM_LocalSpace_Data,
            reference_element_data.FEM_ReferenceElement_Data.rt_triangle_reference_element_data.Quadrature.ReferenceSegmentQuadrature,
            local_space_data.FEM_Geometry);

        const Eigen::VectorXd &weights = internal_quadrature.Weights;
        Eigen::MatrixXd monomials_values =
            monomials_2D.Vander(reference_element_data.FEM_ReferenceElement_Data.rt_triangle_reference_element_data.monomials_2D_data, internal_quadrature.Points, Eigen::Vector3d::Zero(), 1.0);


        const auto monomials_values_der =
            monomials_2D.VanderDerivatives(reference_element_data.FEM_ReferenceElement_Data.rt_triangle_reference_element_data.monomials_2D_data, monomials_values, 1.0);


        Eigen::VectorXd pressure = monomials_values.col(1) + monomials_values.col(2) + 0.5 * monomials_values.col(0);
        std::vector<Eigen::VectorXd> velocity = {-(monomials_values_der[0].col(1) + monomials_values_der[0].col(2)),
                                                 -(monomials_values_der[1].col(1) + monomials_values_der[1].col(2))};

        Eigen::VectorXd source_term = pressure;
        Eigen::VectorXd Dmatrix_pressure = (pressureBasisValues.transpose() * weights.asDiagonal() * pressureBasisValues).llt().solve(pressureBasisValues.transpose() * weights.asDiagonal() * source_term);

        for(unsigned int d = 0; d < dofs_data[1].CellsGlobalDOFs[2].at(c).size(); d++)
        {
            const auto global_dof_i = dofs_data[1].CellsGlobalDOFs[2].at(c).at(d);
            const auto local_dof_i =
                dofs_data[1].CellsDOFs.at(global_dof_i.Dimension).at(global_dof_i.CellIndex).at(global_dof_i.DOFIndex);

            global_solution.SetValue(count_dofs.offsets_DOFs[1] + local_dof_i.Global_Index, Dmatrix_pressure(d));
        }

        Eigen::MatrixXd vander_internal =
            monomials_2D.Vander(reference_element_data.FEM_ReferenceElement_Data.rt_triangle_reference_element_data.monomials_2D_data,
                                internal_quadrature.Points,
                                Eigen::Vector3d::Zero(),
                                1.0);

        const auto vander_internal_der = monomials_2D.VanderDerivatives(reference_element_data.FEM_ReferenceElement_Data.rt_triangle_reference_element_data.monomials_2D_data,
                                                                        vander_internal,
                                                                        1.0);

        Eigen::VectorXd pressure_internal = vander_internal.col(1) + vander_internal.col(2) + 0.5 * vander_internal.col(0);
        std::vector<Eigen::VectorXd> velocity_internal = {-(vander_internal_der[0].col(1) + vander_internal_der[0].col(2)),
                                                          -(vander_internal_der[1].col(1) + vander_internal_der[1].col(2))};

        Eigen::VectorXd Dmatrix = Eigen::MatrixXd::Zero(reference_element_data.FEM_ReferenceElement_Data.rt_triangle_reference_element_data.reference_element_data_velocity.NumBasisFunctions,
                                                        1);

        const auto polygon_edges_normal = mesh_geometric_data.Cell2DsEdgeNormals.at(c);
        for (unsigned int e = 0; e < 3; e++)
        {
            const auto direction = mesh_geometric_data.Cell2DsEdgeDirections.at(c)[e] ? 1.0 : -1.0;


            // Compute dofs of monomials
            Eigen::MatrixXd vander_boundary = monomials_2D.Vander(reference_element_data.FEM_ReferenceElement_Data.rt_triangle_reference_element_data.monomials_2D_data,
                                                                  boundary_quadrature[e].Points,
                                                                  Eigen::Vector3d::Zero(),
                                                                  1.0);

            const auto vander_boundary_der = monomials_2D.VanderDerivatives(reference_element_data.FEM_ReferenceElement_Data.rt_triangle_reference_element_data.monomials_2D_data,
                                                                            vander_boundary,
                                                                            1.0);

            Eigen::VectorXd pressure_boundary = vander_boundary.col(1) + vander_boundary.col(2) + 0.5 * vander_boundary.col(0);
            std::vector<Eigen::VectorXd> velocity_boundary = {-(vander_boundary_der[0].col(1) + vander_boundary_der[0].col(2)),
                                                              -(vander_boundary_der[1].col(1) + vander_boundary_der[1].col(2))};

            Dmatrix.block(e * (reference_element_data.Order + 1),
                          0,
                          (reference_element_data.Order + 1),
                          1) =
                direction * reference_element_data.FEM_ReferenceElement_Data.rt_triangle_reference_element_data.VanderBoundary1D.transpose() *
                    (polygon_edges_normal(0, e) * boundary_quadrature[e].Weights).asDiagonal() *
                    velocity_boundary[0] +
                direction * reference_element_data.FEM_ReferenceElement_Data.rt_triangle_reference_element_data.VanderBoundary1D.transpose() *
                    (polygon_edges_normal(1, e) * boundary_quadrature[e].Weights).asDiagonal() *
                    velocity_boundary[1];

        }

        if (reference_element_data.Order > 0)
        {
            Dmatrix.block(3 * (reference_element_data.Order + 1), 0, reference_element_data.FEM_ReferenceElement_Data.rt_triangle_reference_element_data.Nkm1, 1) =
                reference_element_data.FEM_ReferenceElement_Data.rt_triangle_reference_element_data.reference_element_data_pressure.ReferenceBasisFunctionValues
                    .leftCols(reference_element_data.FEM_ReferenceElement_Data.rt_triangle_reference_element_data.Nkm1)
                    .transpose() *
                reference_element_data.FEM_ReferenceElement_Data.rt_triangle_reference_element_data.Quadrature.ReferenceTriangleQuadrature.Weights.asDiagonal() *
                local_space.MapInvVelocityValuesVect(
                    local_space_data.FEM_LocalSpace_Data,
                    velocity_internal)[0];

            Dmatrix.block(3 * (reference_element_data.Order + 1) + reference_element_data.FEM_ReferenceElement_Data.rt_triangle_reference_element_data.Nkm1, 0, reference_element_data.FEM_ReferenceElement_Data.rt_triangle_reference_element_data.Nkm1, 1) =
                reference_element_data.FEM_ReferenceElement_Data.rt_triangle_reference_element_data.reference_element_data_pressure.ReferenceBasisFunctionValues
                    .leftCols(reference_element_data.FEM_ReferenceElement_Data.rt_triangle_reference_element_data.Nkm1)
                    .transpose() *
                reference_element_data.FEM_ReferenceElement_Data.rt_triangle_reference_element_data.Quadrature.ReferenceTriangleQuadrature.Weights.asDiagonal() *
                local_space.MapInvVelocityValuesVect(
                    local_space_data.FEM_LocalSpace_Data,
                    velocity_internal)[1];
        }


        for(unsigned int d = 0; d < dofs_data[0].CellsGlobalDOFs[2].at(c).size(); d++)
        {
            const auto global_dof_i = dofs_data[0].CellsGlobalDOFs[2].at(c).at(d);
            const auto local_dof_i =
                dofs_data[0].CellsDOFs.at(global_dof_i.Dimension).at(global_dof_i.CellIndex).at(global_dof_i.DOFIndex);

            switch (local_dof_i.Type)
            {
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
                global_neumann_solution.SetValue(local_dof_i.Global_Index, Dmatrix(d));
                break;
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                global_solution.SetValue(local_dof_i.Global_Index, Dmatrix(d));
                break;
            default:
                throw std::runtime_error("Unknown DOF Type");
            }
        }
    }

    std::cout << global_neumann_solution << std::endl;


    Polydim::PDETools::Equations::EllipticEquation equation;

    for(unsigned int c = 0; c < mesh.Cell2DTotalNumber(); c++)
    {
        const auto local_space_data = Polydim::PDETools::LocalSpace_MCC_2D::CreateLocalSpace(geometry_utilities.Tolerance1D(),
                                                                                             geometry_utilities.Tolerance2D(),
                                                                                             mesh_geometric_data,
                                                                                             c,
                                                                                             reference_element_data);


        const auto velocity_basis_functions_values =
            Polydim::PDETools::LocalSpace_MCC_2D::VelocityBasisFunctionsValues(reference_element_data, local_space_data);
        const auto velocity_basis_functions_divergence_values =
            Polydim::PDETools::LocalSpace_MCC_2D::VelocityBasisFunctionsDivergenceValues(reference_element_data, local_space_data);
        const auto pressure_basis_functions_values =
            Polydim::PDETools::LocalSpace_MCC_2D::PressureBasisFunctionsValues(reference_element_data, local_space_data);

        const auto cell2D_internal_quadrature =
            Polydim::PDETools::LocalSpace_MCC_2D::InternalQuadrature(reference_element_data, local_space_data);


        const Eigen::VectorXd &weights = cell2D_internal_quadrature.Weights;


        const auto diffusion_term_values = inverse_diffusion_term(cell2D_internal_quadrature.Points);
        const Eigen::VectorXd source_term_values = source_term(cell2D_internal_quadrature.Points);
        const auto reaction_term_values = reaction_term(cell2D_internal_quadrature.Points);

        auto local_A = equation.ComputeCellDiffusionMatrix(diffusion_term_values,
                                                           velocity_basis_functions_values,
                                                           cell2D_internal_quadrature.Weights);

        const auto local_M = equation.ComputeCellReactionMatrix(reaction_term_values,
                                                                pressure_basis_functions_values,
                                                                cell2D_internal_quadrature.Weights);

        const Eigen::MatrixXd local_B = pressure_basis_functions_values.transpose() *
                                        cell2D_internal_quadrature.Weights.asDiagonal() * velocity_basis_functions_divergence_values;


        const auto local_count_dofs = Polydim::PDETools::Assembler_Utilities::local_count_dofs<2>(c, dofs_data);
        const unsigned int num_local_dofs_pressure = dofs_data[1].CellsGlobalDOFs[2].at(c).size();

        const Eigen::VectorXd dofs_values =
            PDETools::Assembler_Utilities::global_solution_to_local_solution<2>(c,
                                                                                dofs_data,
                                                                                local_count_dofs.num_total_dofs,
                                                                                local_count_dofs.offsets_DOFs,
                                                                                count_dofs.offsets_DOFs,
                                                                                count_dofs.offsets_Strongs,
                                                                                global_solution,
                                                                                global_neumann_solution);

        const Eigen::VectorXd Dmatrix =
            dofs_values.segment(0, local_count_dofs.num_total_dofs - num_local_dofs_pressure);
        const Eigen::VectorXd Dmatrix_pressure =
            dofs_values.segment(local_count_dofs.num_total_dofs - num_local_dofs_pressure, num_local_dofs_pressure);


        Eigen::MatrixXd global_matrix = Eigen::MatrixXd::Zero(5, 5);
        Eigen::MatrixXd neumann_matrix = Eigen::MatrixXd::Zero(5, 6);
        Eigen::VectorXd global_rhs = Eigen::VectorXd::Zero(5);
        Eigen::VectorXd local_solution = Eigen::VectorXd::Zero(5);
        Eigen::VectorXd local_neuman_solution = Eigen::VectorXd::Zero(6);

        local_neuman_solution = Dmatrix.topRows(6);
        local_solution << Dmatrix.bottomRows(2), Dmatrix_pressure;

        global_rhs.bottomRows(3) = pressure_basis_functions_values.transpose() * weights.asDiagonal() * source_term_values;


        neumann_matrix.block(0, 0, 2, 6) = local_A.block(6, 0, 2, 6);
        neumann_matrix.block(2, 0, 3, 6) = pressure_basis_functions_values.transpose() * weights.asDiagonal() * velocity_basis_functions_divergence_values.leftCols(6);

        global_matrix.block(0, 0, 2, 2) = local_A.block(6, 6, 2, 2);
        global_matrix.block(0, 2, 2, 3) = - velocity_basis_functions_divergence_values.rightCols(2).transpose() * weights.asDiagonal()  * pressure_basis_functions_values;
        global_matrix.block(2, 0, 3, 2) = pressure_basis_functions_values.transpose() * weights.asDiagonal() * velocity_basis_functions_divergence_values.rightCols(2);
        global_matrix.block(2, 2, 3, 3) = local_M;

        const double max_error = (global_matrix * local_solution - (global_rhs - neumann_matrix * local_neuman_solution)).cwiseAbs().maxCoeff();

        std::cout << std::scientific;
        std::cout.precision(16);
        std::cout << "c: " << c << " " << max_error << std::endl;

        ASSERT_TRUE(max_error < 1.0e-10);
    }

}

} // namespace UnitTesting
} // namespace Polydim

#endif
