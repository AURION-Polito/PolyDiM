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
#include "Eigen_SparseArray.hpp"
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
// ***************************************************************************
std::vector<Eigen::VectorXd> MapInvVelocityValuesVect(
    const Polydim::FEM::MCC::FEM_Triangle_RT_MCC_2D_LocalSpace_Data &local_space,
    const std::vector<Eigen::VectorXd> &values)
{
    std::vector<Eigen::VectorXd> ref_velocity_values(2, Eigen::VectorXd::Zero(values[0].size()));
    ref_velocity_values[0] = local_space.MapData.DetB *
                             (local_space.MapData.BInv(0, 0) * values[0] + local_space.MapData.BInv(0, 1) * values[1]);
    ref_velocity_values[1] = local_space.MapData.DetB *
                             (local_space.MapData.BInv(1, 0) * values[0] + local_space.MapData.BInv(1, 1) * values[1]);

    return ref_velocity_values;
}
// ***************************************************************************
std::vector<Eigen::MatrixXd> MapInvVelocityValues(
    const Polydim::FEM::MCC::FEM_Triangle_RT_MCC_2D_LocalSpace_Data &local_space,
    const std::vector<Eigen::MatrixXd> &values)
{
    std::vector<Eigen::MatrixXd> ref_velocity_values(2, Eigen::MatrixXd::Zero(values[0].rows(), values[0].cols()));
    ref_velocity_values[0] = local_space.MapData.DetB *
                             (local_space.MapData.BInv(0, 0) * values[0] + local_space.MapData.BInv(0, 1) * values[1]);
    ref_velocity_values[1] = local_space.MapData.DetB *
                             (local_space.MapData.BInv(1, 0) * values[0] + local_space.MapData.BInv(1, 1) * values[1]);

    return ref_velocity_values;
}
// ***************************************************************************
TEST(Test_FEM_Triangle_RT_MCC_2D, Test_FEM_Triangle_RT_MCC_2D_Reference_Element)
{
    const Polydim::FEM::MCC::FEM_Triangle_RT_MCC_2D_ReferenceElement reference_element;

    const auto referenceQuadrature = Gedim::Quadrature::Quadrature_Gauss2D_Triangle::FillPointsAndWeights(21);
    const Eigen::MatrixXd &referenceQuadraturePoints = referenceQuadrature.Points;

    Polydim::Utilities::Monomials_2D monomials_2D;
    Polydim::Utilities::Monomials_1D monomials_1D;

    std::vector<std::array<bool, 3>> edge_directions = {{true, true, true},
                                                        {true, true, false},
                                                        {false, true, true},
                                                        {true, false, true},
                                                        {false, false, true},
                                                        {false, true, false},
                                                        {true, false, false},
                                                        {false, false, false}};

    for(unsigned int i = 0; i < edge_directions.size(); i++)
    {
        for (unsigned int o = 0; o < 5; o++)
        {
            const auto reference_element_data = reference_element.Create(o);

            const Eigen::MatrixXd pressureBasisValues =
                reference_element.EvaluatePressureBasisFunctions(referenceQuadraturePoints, reference_element_data);
            const std::vector<Eigen::MatrixXd> velocityBasisValues =
                reference_element.EvaluateVelocityBasisFunctions(referenceQuadraturePoints,
                                                                 reference_element_data.reference_element_data_velocity.basis_functions.at(edge_directions[i]).MonomialsCoefficients,
                                                                 reference_element_data);
            const Eigen::MatrixXd diverVelocityBasisFunctions =
                reference_element.EvaluateVelocityBasisFunctionsDivergence(referenceQuadraturePoints,
                                                                           reference_element_data.reference_element_data_velocity.basis_functions.at(edge_directions[i]).MonomialsCoefficients,
                                                                           reference_element_data);


            // Compute dofs of monomials
            Eigen::MatrixXd vander_boundary = monomials_2D.Vander(reference_element_data.monomials_2D_data,
                                                                  reference_element_data.BoundaryQuadrature.at(edge_directions[i]).Quadrature.Points,
                                                                  Eigen::Vector3d::Zero(),
                                                                  1.0);

            Eigen::MatrixXd vander_internal =
                monomials_2D.Vander(reference_element_data.monomials_2D_data,
                                    reference_element_data.Quadrature.ReferenceTriangleQuadrature.Points,
                                    Eigen::Vector3d::Zero(),
                                    1.0);


            for(unsigned int e = 0; e < 3; e++)
            {
                const double direction = edge_directions[i][e] ? 1.0 : -1.0;

                const std::vector<Eigen::MatrixXd> velocityBasisValues_edge =
                    reference_element.EvaluateVelocityBasisFunctions(reference_element_data.BoundaryQuadrature.at(edge_directions[i]).Quadrature.Points.middleCols(e * reference_element_data.reference_element_data_velocity.NumDofs1D,
                                                                                                                                                                   reference_element_data.reference_element_data_velocity.NumDofs1D),
                                                                     reference_element_data.reference_element_data_velocity.basis_functions.at(edge_directions[i]).MonomialsCoefficients,
                                                                     reference_element_data);

                Eigen::MatrixXd quantity = (velocityBasisValues_edge[0] * reference_element_data.EdgeNormals(0, e)
                                            + velocityBasisValues_edge[1] * reference_element_data.EdgeNormals(1, e)).transpose()
                                           * reference_element_data.BoundaryQuadrature.at(edge_directions[i]).Quadrature.Weights.segment(e * reference_element_data.reference_element_data_velocity.NumDofs1D,
                                                                                                                                         reference_element_data.reference_element_data_velocity.NumDofs1D).asDiagonal()
                                           *  reference_element_data.VanderBoundary1D;
                ASSERT_TRUE((quantity.sum() - direction * (o + 1)) < 1.0e-10);
            }

            Eigen::MatrixXd Dmatrix = Eigen::MatrixXd::Zero(reference_element_data.reference_element_data_velocity.NumBasisFunctions,
                                                            2 * reference_element_data.Nk);

            for (unsigned int e = 0; e < 3; e++)
            {
                const double direction = edge_directions[i][e] ? 1.0 : -1.0;

                Dmatrix.block(e * (reference_element_data.Order + 1),
                              0,
                              (reference_element_data.Order + 1),
                              reference_element_data.Nk) =
                    direction * reference_element_data.VanderBoundary1D.transpose() *
                    reference_element_data.BoundaryQuadrature.at(edge_directions[i]).WeightsTimesNormal[0]
                        .segment(e * (reference_element_data.Order + 1), (reference_element_data.Order + 1))
                        .asDiagonal() *
                    vander_boundary.middleRows(e * (reference_element_data.Order + 1), (reference_element_data.Order + 1));

                Dmatrix.block(e * (reference_element_data.Order + 1),
                              reference_element_data.Nk,
                              (reference_element_data.Order + 1),
                              reference_element_data.Nk) =
                    direction * reference_element_data.VanderBoundary1D.transpose() *
                    reference_element_data.BoundaryQuadrature.at(edge_directions[i]).WeightsTimesNormal[1]
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
}

// TEST(Test_FEM_Triangle_RT_MCC_2D, Test_FEM_Triangle_RT_MCC_2D_Reference_Element_Patch)
// {
//     const Polydim::FEM::MCC::FEM_Triangle_RT_MCC_2D_ReferenceElement reference_element;

//     const auto referenceQuadrature = Gedim::Quadrature::Quadrature_Gauss2D_Triangle::FillPointsAndWeights(4);
//     const Eigen::MatrixXd &referenceQuadraturePoints = referenceQuadrature.Points;

//     Polydim::Utilities::Monomials_2D monomials_2D;
//     Polydim::Utilities::Monomials_1D monomials_1D;

//     for (unsigned int o = 1; o < 2; o++)
//     {
//         const auto reference_element_data = reference_element.Create(o);

//         const Eigen::MatrixXd pressureBasisValues =
//             reference_element.EvaluatePressureBasisFunctions(referenceQuadraturePoints, reference_element_data);
//         const std::vector<Eigen::MatrixXd> velocityBasisValues =
//             reference_element.EvaluateVelocityBasisFunctions(referenceQuadraturePoints, reference_element_data);
//         const Eigen::MatrixXd diverVelocityBasisFunctions =
//             reference_element.EvaluateVelocityBasisFunctionsDivergence(referenceQuadraturePoints, reference_element_data);

//         // Compute dofs of monomials
//         Eigen::MatrixXd vander_boundary = monomials_2D.Vander(reference_element_data.monomials_2D_data,
//                                                               reference_element_data.BoundaryQuadrature.Quadrature.Points,
//                                                               Eigen::Vector3d::Zero(),
//                                                               1.0);

//         const auto vander_boundary_der = monomials_2D.VanderDerivatives(reference_element_data.monomials_2D_data,
//                                                                         vander_boundary,
//                                                                         1.0);

//         Eigen::MatrixXd vander_internal =
//             monomials_2D.Vander(reference_element_data.monomials_2D_data,
//                                 reference_element_data.Quadrature.ReferenceTriangleQuadrature.Points,
//                                 Eigen::Vector3d::Zero(),
//                                 1.0);

//         const auto vander_internal_der = monomials_2D.VanderDerivatives(reference_element_data.monomials_2D_data,
//                                                                         vander_internal,
//                                                                         1.0);

//         Eigen::VectorXd pressure_boundary = vander_boundary.col(1) + vander_boundary.col(2) + 0.5 * vander_boundary.col(0);
//         std::vector<Eigen::VectorXd> velocity_boundary = {-(vander_boundary_der[0].col(1) + vander_boundary_der[0].col(2)),
//                                                           -(vander_boundary_der[1].col(1) + vander_boundary_der[1].col(2))};

//         Eigen::VectorXd pressure_internal = vander_internal.col(1) + vander_internal.col(2) + 0.5 * vander_internal.col(0);
//         std::vector<Eigen::VectorXd> velocity_internal = {-(vander_internal_der[0].col(1) + vander_internal_der[0].col(2)),
//                                                           -(vander_internal_der[1].col(1) + vander_internal_der[1].col(2))};

//         Eigen::MatrixXd Dmatrix = Eigen::MatrixXd::Zero(reference_element_data.reference_element_data_velocity.NumBasisFunctions,
//                                                         1);

//         for (unsigned int e = 0; e < 3; e++)
//         {
//             Dmatrix.block(e * (reference_element_data.Order + 1),
//                           0,
//                           (reference_element_data.Order + 1),
//                           1) =
//                 reference_element_data.VanderBoundary1D.transpose() *
//                     reference_element_data.BoundaryQuadrature.WeightsTimesNormal[0]
//                         .segment(e * (reference_element_data.Order + 1), (reference_element_data.Order + 1))
//                         .asDiagonal() *
//                     velocity_boundary[0].middleRows(e * (reference_element_data.Order + 1), (reference_element_data.Order + 1)) +
//                 reference_element_data.VanderBoundary1D.transpose() *
//                     reference_element_data.BoundaryQuadrature.WeightsTimesNormal[1]
//                         .segment(e * (reference_element_data.Order + 1), (reference_element_data.Order + 1))
//                         .asDiagonal() *
//                     velocity_boundary[1].middleRows(e * (reference_element_data.Order + 1), (reference_element_data.Order + 1));

//         }

//         if (reference_element_data.Order > 0)
//         {
//             Dmatrix.block(3 * (reference_element_data.Order + 1), 0, reference_element_data.Nkm1, 1) =
//                 reference_element_data.reference_element_data_pressure.ReferenceBasisFunctionValues
//                     .leftCols(reference_element_data.Nkm1)
//                     .transpose() *
//                 reference_element_data.Quadrature.ReferenceTriangleQuadrature.Weights.asDiagonal() * velocity_internal[0];

//             Dmatrix.block(3 * (reference_element_data.Order + 1) + reference_element_data.Nkm1, 0, reference_element_data.Nkm1, 1) =
//                 reference_element_data.reference_element_data_pressure.ReferenceBasisFunctionValues
//                     .leftCols(reference_element_data.Nkm1)
//                     .transpose() *
//                 reference_element_data.Quadrature.ReferenceTriangleQuadrature.Weights.asDiagonal() * velocity_internal[1];
//         }

//         const Eigen::VectorXd &weights = reference_element_data.Quadrature.ReferenceTriangleQuadrature.Weights;
//         Eigen::MatrixXd monomials_values =
//             monomials_2D.Vander(reference_element_data.monomials_2D_data, referenceQuadrature.Points, Eigen::Vector3d::Zero(), 1.0);


//         const auto monomials_values_der =
//             monomials_2D.VanderDerivatives(reference_element_data.monomials_2D_data, monomials_values, 1.0);


//         Eigen::VectorXd pressure = monomials_values.col(1) + monomials_values.col(2) + 0.5 * monomials_values.col(0);
//         std::vector<Eigen::VectorXd> velocity = {-(monomials_values_der[0].col(1) + monomials_values_der[0].col(2)),
//                                                  -(monomials_values_der[1].col(1) + monomials_values_der[1].col(2))};

//         Eigen::VectorXd source_term = pressure;

//         const Eigen::VectorXd local_rhs_A  = velocityBasisValues[0].transpose() * weights.asDiagonal() * velocity[0]
//                                             + velocityBasisValues[1].transpose() * weights.asDiagonal() * velocity[1];
//         const auto global_A = velocityBasisValues[0].transpose() * weights.asDiagonal() * velocityBasisValues[0] +
//                               velocityBasisValues[1].transpose() * weights.asDiagonal() * velocityBasisValues[1];

//         Eigen::VectorXd Dmatrix_pressure = Eigen::VectorXd::Zero(3);
//         Dmatrix_pressure << 0.5 + reference_element_data.monomials_2D_center(0) + reference_element_data.monomials_2D_center(1), 1.0, 1.0;

//         Eigen::VectorXd discrete_pressure = pressureBasisValues * Dmatrix_pressure;

//         // const Eigen::VectorXd velocity_error = global_A.block(6, 6, 2, 2) * Dmatrix.bottomRows(2) - local_rhs_A.bottomRows(2)  - (- global_A.block(6, 0, 2, 6) * Dmatrix.topRows(6));
//         const Eigen::VectorXd value_error = diverVelocityBasisFunctions.rightCols(2).transpose() * weights.asDiagonal()  * pressureBasisValues * Dmatrix_pressure;
//         const Eigen::VectorXd velocity_error = global_A.block(6, 6, 2, 2) * Dmatrix.bottomRows(2) - diverVelocityBasisFunctions.rightCols(2).transpose() * weights.asDiagonal()  * pressureBasisValues * Dmatrix_pressure
//                                                - (- global_A.block(6, 0, 2, 6) * Dmatrix.topRows(6));
//         const Eigen::VectorXd div_velocity_error = pressureBasisValues.transpose() * weights.asDiagonal() * diverVelocityBasisFunctions.rightCols(2) * Dmatrix.bottomRows(2)
//                                                    + pressureBasisValues.transpose() * weights.asDiagonal()  * pressureBasisValues * Dmatrix_pressure
//                                                    - (pressureBasisValues.transpose() * weights.asDiagonal()  * source_term
//                                                       - pressureBasisValues.transpose() * weights.asDiagonal() * diverVelocityBasisFunctions.leftCols(6) * Dmatrix.topRows(6));

//         Eigen::MatrixXd global_matrix = Eigen::MatrixXd::Zero(5, 5);
//         Eigen::MatrixXd neumann_matrix = Eigen::MatrixXd::Zero(5, 6);
//         Eigen::VectorXd global_rhs = Eigen::VectorXd::Zero(5);
//         Eigen::VectorXd solution = Eigen::VectorXd::Zero(5);
//         Eigen::VectorXd neuman_solution = Eigen::VectorXd::Zero(6);
//         neuman_solution = Dmatrix.topRows(6);
//         solution << Dmatrix.bottomRows(2), Dmatrix_pressure;
//         global_rhs.bottomRows(3) = pressureBasisValues.transpose() * weights.asDiagonal()  * source_term;
//         neumann_matrix.block(2, 0, 3, 6) = pressureBasisValues.transpose() * weights.asDiagonal() * diverVelocityBasisFunctions.leftCols(6);

//         global_matrix.block(2, 0, 3, 2) = pressureBasisValues.transpose() * weights.asDiagonal() * diverVelocityBasisFunctions.rightCols(2);
//         global_matrix.block(2, 2, 3, 3) = pressureBasisValues.transpose() * weights.asDiagonal()  * pressureBasisValues;


//         const double max_error_vel = velocity_error.cwiseAbs().maxCoeff();

//         const double max_error = (global_matrix * solution - (global_rhs - neumann_matrix * neuman_solution)).cwiseAbs().maxCoeff();
//         const double max_error_div = div_velocity_error.cwiseAbs().maxCoeff();

//         ASSERT_TRUE(max_error_vel < 1.0e-10);
//         ASSERT_TRUE(max_error < 1.0e-10);
//         ASSERT_TRUE(max_error_div < 1.0e-10);
//     }
// }

TEST(Test_FEM_Triangle_RT_MCC_2D, Test_FEM_Triangle_RT_MCC_2D_Local_Space)
{
    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = 1.0e-8;
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    const auto poligon_vertices = geometry_utilities.CreateTriangle(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                    Eigen::Vector3d(0.2, -0.3, 0.0),
                                                                    Eigen::Vector3d(0.4, 0.6, 0.0));

    const auto polygon_edges_tangent = geometry_utilities.PolygonEdgeTangents(poligon_vertices);
    const auto polygon_edges_length = geometry_utilities.PolygonEdgeLengths(poligon_vertices);
    const auto polygon_edges_normal = geometry_utilities.PolygonEdgeNormals(poligon_vertices);

    std::vector<std::array<bool, 3>> edge_directions = {{true, true, true}, // 0
                                                        {true, true, false}, // 1
                                                        {false, true, true}, // 2
                                                        {true, false, true}, // 3
                                                        {false, false, true}, // 4
                                                        {false, true, false}, // 5
                                                        {true, false, false},
                                                        {false, false, false}};

    const Polydim::FEM::MCC::FEM_MCC_2D_ReferenceElement reference_element;
    const Polydim::FEM::MCC::FEM_MCC_2D_LocalSpace local_space;

    Polydim::Utilities::Monomials_2D monomials_2D;
    Polydim::Utilities::Monomials_1D monomials_1D;
    Polydim::VEM::Quadrature::VEM_Quadrature_2D quadrature;

    const auto referenceQuadrature = Gedim::Quadrature::Quadrature_Gauss2D_Triangle::FillPointsAndWeights(21);
    const auto quadrature_data = quadrature.PolygonInternalQuadrature(referenceQuadrature, {poligon_vertices});

    for(unsigned int p = 0; p < edge_directions.size(); p++)
    {
        for (unsigned int o = 1; o < 2; o++)
        {
            std::vector<bool>  polygon_edge_directions(3);
            for(unsigned int d = 0; d < 3; d++)
                polygon_edge_directions[d] = edge_directions[p][d];

            Polydim::FEM::MCC::FEM_MCC_2D_Polygon_Geometry polygon_geometry = {geometry_utilities.Tolerance1D(),
                                                                               geometry_utilities.Tolerance2D(),
                                                                               poligon_vertices,
                                                                               polygon_edges_length,
                                                                               polygon_edge_directions,
                                                                               polygon_edges_tangent,
                                                                               polygon_edges_normal};

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

            const auto internal_quadrature = local_space_data.InternalQuadrature;
            const auto boundary_quadrature = local_space_data.BoundaryQuadrature;


            for(unsigned int e = 0; e < 3; e++)
            {
                const double direction = edge_directions[p][e] ? 1.0 : -1.0;

                const std::vector<Eigen::MatrixXd> velocityBasisValues_edge =
                    local_space.ComputeVelocityBasisFunctionsValues(reference_element_data, local_space_data, boundary_quadrature[e].Points);

                if(p == 0 && e == 1)
                    std::cout << (velocityBasisValues_edge[0] * polygon_edges_normal(0, e) + velocityBasisValues_edge[1] * polygon_edges_normal(1, e)) << std::endl;

                Eigen::MatrixXd quantity = (velocityBasisValues_edge[0] * polygon_edges_normal(0, e) + velocityBasisValues_edge[1] * polygon_edges_normal(1, e)).transpose()
                                           * boundary_quadrature[e].Weights.asDiagonal() *  reference_element_data.rt_triangle_reference_element_data.VanderBoundary1D;

                std::cout << p << " " << e << std::endl;
                std::cout << quantity.transpose() << std::endl;

                ASSERT_TRUE((quantity.sum() - direction * (o + 1)) < 1.0e-10);
            }

            // Compute dofs of monomials
            Eigen::MatrixXd vander_internal = monomials_2D.Vander(reference_element_data.rt_triangle_reference_element_data.monomials_2D_data,
                                                                  internal_quadrature.Points,
                                                                  Eigen::Vector3d::Zero(),
                                                                  1.0);

            Eigen::MatrixXd Dmatrix = Eigen::MatrixXd::Zero(local_space_data.NumVelocityBasisFunctions,
                                                            2 * reference_element_data.rt_triangle_reference_element_data.Nk);

            for (unsigned int e = 0; e < 3; e++)
            {

                const double direction = edge_directions[p][e] ? 1.0 : -1.0;

                Eigen::MatrixXd vander_boundary =
                    monomials_2D.Vander(reference_element_data.rt_triangle_reference_element_data.monomials_2D_data,
                                        boundary_quadrature[e].Points,
                                        Eigen::Vector3d::Zero(),
                                        1.0);

                Dmatrix.block(e * (reference_element_data.Order + 1),
                              0,
                              (reference_element_data.Order + 1),
                              reference_element_data.rt_triangle_reference_element_data.Nk) =
                    direction * reference_element_data.rt_triangle_reference_element_data.VanderBoundary1D.transpose() *
                    (polygon_edges_normal(0, e) * boundary_quadrature[e].Weights).asDiagonal() * vander_boundary;

                Dmatrix.block(e * (reference_element_data.Order + 1),
                              reference_element_data.rt_triangle_reference_element_data.Nk,
                              (reference_element_data.Order + 1),
                              reference_element_data.rt_triangle_reference_element_data.Nk) =
                    direction * reference_element_data.rt_triangle_reference_element_data.VanderBoundary1D.transpose() *
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
                    MapInvVelocityValues(
                        local_space_data.rt_triangle_local_space_data,
                        {vander_internal, Eigen::MatrixXd::Zero(vander_internal.rows(), vander_internal.cols())})[0];


                Dmatrix.block(3 * (reference_element_data.Order + 1) + reference_element_data.rt_triangle_reference_element_data.Nkm1,
                              0,
                              reference_element_data.rt_triangle_reference_element_data.Nkm1,
                              reference_element_data.rt_triangle_reference_element_data.Nk) =
                    reference_element_data.rt_triangle_reference_element_data.reference_element_data_pressure
                        .ReferenceBasisFunctionValues
                        .leftCols(reference_element_data.rt_triangle_reference_element_data.Nkm1)
                        .transpose() *
                    reference_element_data.rt_triangle_reference_element_data.Quadrature.ReferenceTriangleQuadrature.Weights.asDiagonal() *
                    MapInvVelocityValues(
                        local_space_data.rt_triangle_local_space_data,
                        {vander_internal, Eigen::MatrixXd::Zero(vander_internal.rows(), vander_internal.cols())})[1];


                Dmatrix.block(3 * (reference_element_data.Order + 1),
                              reference_element_data.rt_triangle_reference_element_data.Nk,
                              reference_element_data.rt_triangle_reference_element_data.Nkm1,
                              reference_element_data.rt_triangle_reference_element_data.Nk) =
                    reference_element_data.rt_triangle_reference_element_data.reference_element_data_pressure
                        .ReferenceBasisFunctionValues
                        .leftCols(reference_element_data.rt_triangle_reference_element_data.Nkm1)
                        .transpose() *
                    reference_element_data.rt_triangle_reference_element_data.Quadrature.ReferenceTriangleQuadrature.Weights.asDiagonal() *
                    MapInvVelocityValues(
                        local_space_data.rt_triangle_local_space_data,
                        {Eigen::MatrixXd::Zero(vander_internal.rows(), vander_internal.cols()), vander_internal})[0];

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
                    MapInvVelocityValues(
                        local_space_data.rt_triangle_local_space_data,
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
}

TEST(Test_FEM_Triangle_RT_MCC_2D, Test_FEM_Triangle_RT_MCC_2D_Local_Space_2)
{
    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = 1.0e-8;
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    const auto poligon_vertices = geometry_utilities.CreateTriangle(Eigen::Vector3d(0.2, -0.3, 0.0),
                                                                    Eigen::Vector3d(1.0, 0.6, 0.0),
                                                                    Eigen::Vector3d(0.4, 0.6, 0.0));
    std::vector<bool> polygon_edges_direction(3, true);
    polygon_edges_direction[2] = false;
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

    for (unsigned int o = 1; o < 2; o++)
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

        const auto internal_quadrature = local_space_data.InternalQuadrature;
        const auto boundary_quadrature = local_space_data.BoundaryQuadrature;


        for(unsigned int e = 0; e < 3; e++)
        {
            const double direction = polygon_edges_direction[e] ? 1.0 : -1.0;

            const std::vector<Eigen::MatrixXd> velocityBasisValues_edge =
                local_space.ComputeVelocityBasisFunctionsValues(reference_element_data, local_space_data, boundary_quadrature[e].Points);

            if(e == 2)
                std::cout << (velocityBasisValues_edge[0] * polygon_edges_normal(0, e) + velocityBasisValues_edge[1] * polygon_edges_normal(1, e)) << std::endl;

            Eigen::MatrixXd quantity = (velocityBasisValues_edge[0] * polygon_edges_normal(0, e) + velocityBasisValues_edge[1] * polygon_edges_normal(1, e)).transpose()
                                       * boundary_quadrature[e].Weights.asDiagonal() *  reference_element_data.rt_triangle_reference_element_data.VanderBoundary1D;
            ASSERT_TRUE((quantity.sum() - direction * (o + 1)) < 1.0e-10);
        }

        // Compute dofs of monomials
        Eigen::MatrixXd vander_internal = monomials_2D.Vander(reference_element_data.rt_triangle_reference_element_data.monomials_2D_data,
                                                              internal_quadrature.Points,
                                                              Eigen::Vector3d::Zero(),
                                                              1.0);

        Eigen::MatrixXd Dmatrix = Eigen::MatrixXd::Zero(local_space_data.NumVelocityBasisFunctions,
                                                        2 * reference_element_data.rt_triangle_reference_element_data.Nk);

        for (unsigned int e = 0; e < 3; e++)
        {

            double direction = polygon_edges_direction[e] ? 1.0 : -1.0;

            Eigen::MatrixXd vander_boundary =
                monomials_2D.Vander(reference_element_data.rt_triangle_reference_element_data.monomials_2D_data,
                                    boundary_quadrature[e].Points,
                                    Eigen::Vector3d::Zero(),
                                    1.0);

            Dmatrix.block(e * (reference_element_data.Order + 1),
                          0,
                          (reference_element_data.Order + 1),
                          reference_element_data.rt_triangle_reference_element_data.Nk) =
                direction * reference_element_data.rt_triangle_reference_element_data.VanderBoundary1D.transpose() *
                (polygon_edges_normal(0, e) * boundary_quadrature[e].Weights).asDiagonal() * vander_boundary;

            Dmatrix.block(e * (reference_element_data.Order + 1),
                          reference_element_data.rt_triangle_reference_element_data.Nk,
                          (reference_element_data.Order + 1),
                          reference_element_data.rt_triangle_reference_element_data.Nk) =
                direction * reference_element_data.rt_triangle_reference_element_data.VanderBoundary1D.transpose() *
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
                MapInvVelocityValues(
                    local_space_data.rt_triangle_local_space_data,
                    {vander_internal, Eigen::MatrixXd::Zero(vander_internal.rows(), vander_internal.cols())})[0];


            Dmatrix.block(3 * (reference_element_data.Order + 1) + reference_element_data.rt_triangle_reference_element_data.Nkm1,
                          0,
                          reference_element_data.rt_triangle_reference_element_data.Nkm1,
                          reference_element_data.rt_triangle_reference_element_data.Nk) =
                reference_element_data.rt_triangle_reference_element_data.reference_element_data_pressure
                    .ReferenceBasisFunctionValues
                    .leftCols(reference_element_data.rt_triangle_reference_element_data.Nkm1)
                    .transpose() *
                reference_element_data.rt_triangle_reference_element_data.Quadrature.ReferenceTriangleQuadrature.Weights.asDiagonal() *
                MapInvVelocityValues(
                    local_space_data.rt_triangle_local_space_data,
                    {vander_internal, Eigen::MatrixXd::Zero(vander_internal.rows(), vander_internal.cols())})[1];


            Dmatrix.block(3 * (reference_element_data.Order + 1),
                          reference_element_data.rt_triangle_reference_element_data.Nk,
                          reference_element_data.rt_triangle_reference_element_data.Nkm1,
                          reference_element_data.rt_triangle_reference_element_data.Nk) =
                reference_element_data.rt_triangle_reference_element_data.reference_element_data_pressure
                    .ReferenceBasisFunctionValues
                    .leftCols(reference_element_data.rt_triangle_reference_element_data.Nkm1)
                    .transpose() *
                reference_element_data.rt_triangle_reference_element_data.Quadrature.ReferenceTriangleQuadrature.Weights.asDiagonal() *
                MapInvVelocityValues(
                    local_space_data.rt_triangle_local_space_data,
                    {Eigen::MatrixXd::Zero(vander_internal.rows(), vander_internal.cols()), vander_internal})[0];

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
                MapInvVelocityValues(
                    local_space_data.rt_triangle_local_space_data,
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

Eigen::VectorXd exact_pressure(const Eigen::MatrixXd &points)
{
    unsigned int order = 1;
    const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + 0.5;

    Eigen::ArrayXd result = Eigen::ArrayXd::Constant(points.cols(), 1.0);
    for (int i = 0; i < order; i++)
        result = result * polynomial;

    return result;
};

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
                                                                0.25,
                                                                mesh);

    {
        Gedim::MeshUtilities meshUtilities;
        meshUtilities.ExportMeshToVTU(mesh, "./", "Domain_Mesh");
    }

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

        const auto internal_quadrature = local_space_data.FEM_LocalSpace_Data.InternalQuadrature;
        const auto boundary_quadrature = local_space_data.FEM_LocalSpace_Data.BoundaryQuadrature;

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
                MapInvVelocityValuesVect(
                    local_space_data.FEM_LocalSpace_Data.rt_triangle_local_space_data,
                    velocity_internal)[0];

            Dmatrix.block(3 * (reference_element_data.Order + 1) + reference_element_data.FEM_ReferenceElement_Data.rt_triangle_reference_element_data.Nkm1, 0,
                          reference_element_data.FEM_ReferenceElement_Data.rt_triangle_reference_element_data.Nkm1, 1) =
                reference_element_data.FEM_ReferenceElement_Data.rt_triangle_reference_element_data.reference_element_data_pressure.ReferenceBasisFunctionValues
                    .leftCols(reference_element_data.FEM_ReferenceElement_Data.rt_triangle_reference_element_data.Nkm1)
                    .transpose() *
                reference_element_data.FEM_ReferenceElement_Data.rt_triangle_reference_element_data.Quadrature.ReferenceTriangleQuadrature.Weights.asDiagonal() *
                MapInvVelocityValuesVect(
                    local_space_data.FEM_LocalSpace_Data.rt_triangle_local_space_data,
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

    Gedim::Eigen_SparseArray<> globalMatrixA;
    Gedim::Eigen_SparseArray<> neumannMatrixA;
    Gedim::Eigen_Array<> rightHandSide;
    globalMatrixA.SetSize(count_dofs.num_total_dofs, count_dofs.num_total_dofs, Gedim::ISparseArray::SparseArrayTypes::None);
    neumannMatrixA.SetSize(count_dofs.num_total_dofs, count_dofs.num_total_strong);
    rightHandSide.SetSize(count_dofs.num_total_dofs);


    std::vector<Eigen::MatrixXd> local_A(mesh.Cell2DTotalNumber());
    std::vector<std::array<Eigen::MatrixXd, 3>> local_dirichlet(mesh.Cell2DTotalNumber());
    std::vector<Eigen::MatrixXd> local_B(mesh.Cell2DTotalNumber());
    std::vector<Eigen::MatrixXd> local_M(mesh.Cell2DTotalNumber());
    std::vector<Eigen::MatrixXd> local_rhs(mesh.Cell2DTotalNumber());


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


        for(unsigned int e = 0; e < 3; e++)
        {
            const double direction = mesh_geometric_data.Cell2DsEdgeDirections.at(c)[e] ? 1.0 : -1.0;

            const auto exact_sol_values = exact_pressure(local_space_data.FEM_LocalSpace_Data.BoundaryQuadrature[e].Points);

            const std::vector<Eigen::MatrixXd> velocityBasisValues_edge =
                Polydim::PDETools::LocalSpace_MCC_2D::VelocityBasisFunctionsValues(reference_element_data,
                                                                                   local_space_data,
                                                                                   local_space_data.FEM_LocalSpace_Data.BoundaryQuadrature[e].Points);


            std::cout << "c: " << c << " e: " << e << std::endl;
            std::cout << velocityBasisValues_edge[0] * mesh_geometric_data.Cell2DsEdgeNormals.at(c)(0, e)
                             + velocityBasisValues_edge[1] * mesh_geometric_data.Cell2DsEdgeNormals.at(c)(1, e) << std::endl;

            std::cout << "local_space_data.FEM_LocalSpace_Data.BoundaryQuadrature[e].Points" << std::endl;
            std::cout << local_space_data.FEM_LocalSpace_Data.BoundaryQuadrature[e].Points << std::endl;

            std::cout << "exact_sol_values" << std::endl;
            std::cout << exact_sol_values.transpose() << std::endl;


            ASSERT_TRUE((((velocityBasisValues_edge[0] * mesh_geometric_data.Cell2DsEdgeNormals.at(c)(0, e)
                           + velocityBasisValues_edge[1] * mesh_geometric_data.Cell2DsEdgeNormals.at(c)(1, e)).transpose()
                          * local_space_data.FEM_LocalSpace_Data.BoundaryQuadrature[e].Weights.asDiagonal()
                          *  reference_element_data.FEM_ReferenceElement_Data.rt_triangle_reference_element_data.VanderBoundary1D ).sum() - direction * (o + 1)) < 1.0e-12);


            local_dirichlet[c][e] = (velocityBasisValues_edge[0] * mesh_geometric_data.Cell2DsEdgeNormals.at(c)(0, e)
                                     + velocityBasisValues_edge[1] * mesh_geometric_data.Cell2DsEdgeNormals.at(c)(1, e)).transpose()
                                    * local_space_data.FEM_LocalSpace_Data.BoundaryQuadrature[e].Weights.asDiagonal() * exact_sol_values;
        }


        const auto diffusion_term_values = inverse_diffusion_term(cell2D_internal_quadrature.Points);
        const Eigen::VectorXd source_term_values = source_term(cell2D_internal_quadrature.Points);
        const auto reaction_term_values = reaction_term(cell2D_internal_quadrature.Points);

        local_A[c] = equation.ComputeCellDiffusionMatrix(diffusion_term_values,
                                                         velocity_basis_functions_values,
                                                         cell2D_internal_quadrature.Weights);

        local_M[c] = equation.ComputeCellReactionMatrix(reaction_term_values,
                                                        pressure_basis_functions_values,
                                                        cell2D_internal_quadrature.Weights);

        local_B[c] = pressure_basis_functions_values.transpose() *
                     cell2D_internal_quadrature.Weights.asDiagonal() * velocity_basis_functions_divergence_values;


        local_rhs[c] = pressure_basis_functions_values.transpose() * weights.asDiagonal() * source_term_values;
    }

    for(unsigned int c = 0; c < mesh.Cell2DTotalNumber(); c++)
    {
        const auto local_count_dofs = Polydim::PDETools::Assembler_Utilities::local_count_dofs<2>(c, dofs_data);
        const unsigned int num_local_dofs_pressure = dofs_data[1].CellsGlobalDOFs[2].at(c).size();

        Eigen::MatrixXd elemental_matrix = Eigen::MatrixXd::Zero(local_count_dofs.num_total_dofs, local_count_dofs.num_total_dofs);
        Eigen::VectorXd elemental_rhs = Eigen::VectorXd::Zero(local_count_dofs.num_total_dofs);
        elemental_matrix << local_A[c], -(local_B[c]).transpose(), local_B[c], local_M[c];
        elemental_rhs << Eigen::VectorXd::Zero(local_count_dofs.num_total_dofs - num_local_dofs_pressure), local_rhs[c];

        Polydim::PDETools::Assembler_Utilities::local_matrix_to_global_matrix_dofs_data local_matrix_to_global_matrix_dofs_data = {
                                                                                                                                   {std::cref(dofs_data[0]), std::cref(dofs_data[1])},
                                                                                                                                   local_count_dofs.offsets_DOFs,
                                                                                                                                   count_dofs.offsets_DOFs,
                                                                                                                                   count_dofs.offsets_Strongs};

        Polydim::PDETools::Assembler_Utilities::assemble_local_matrix_to_global_matrix<2>(c,
                                                                                          local_matrix_to_global_matrix_dofs_data,
                                                                                          local_matrix_to_global_matrix_dofs_data,
                                                                                          elemental_matrix,
                                                                                          elemental_rhs,
                                                                                          globalMatrixA,
                                                                                          neumannMatrixA,
                                                                                          rightHandSide);
    }


    rightHandSide.Create();
    globalMatrixA.Create();
    neumannMatrixA.Create();

    if (count_dofs.num_total_strong > 0)
        rightHandSide.SubtractionMultiplication(neumannMatrixA, global_neumann_solution);

    rightHandSide.SubtractionMultiplication(globalMatrixA, global_solution);

    std::cout << "rightHandSide.Norm(): " << rightHandSide.Norm() << std::endl;


    {

        const unsigned int num_local_dofs_pressure = dofs_data[1].CellsGlobalDOFs[2].at(0).size();

        const auto local_count_dofs_0 = Polydim::PDETools::Assembler_Utilities::local_count_dofs<2>(0, dofs_data);


        const Eigen::VectorXd dofs_values_0 =
            PDETools::Assembler_Utilities::global_solution_to_local_solution<2>(0,
                                                                                dofs_data,
                                                                                local_count_dofs_0.num_total_dofs,
                                                                                local_count_dofs_0.offsets_DOFs,
                                                                                count_dofs.offsets_DOFs,
                                                                                count_dofs.offsets_Strongs,
                                                                                global_solution,
                                                                                global_neumann_solution);

        const Eigen::VectorXd Dmatrix_0 =
            dofs_values_0.segment(0, local_count_dofs_0.num_total_dofs - num_local_dofs_pressure);
        const Eigen::VectorXd Dmatrix_pressure_0 =
            dofs_values_0.segment(local_count_dofs_0.num_total_dofs - num_local_dofs_pressure, num_local_dofs_pressure);

        const auto local_count_dofs_2 = Polydim::PDETools::Assembler_Utilities::local_count_dofs<2>(2, dofs_data);


        const Eigen::VectorXd dofs_values_2 =
            PDETools::Assembler_Utilities::global_solution_to_local_solution<2>(2,
                                                                                dofs_data,
                                                                                local_count_dofs_2.num_total_dofs,
                                                                                local_count_dofs_2.offsets_DOFs,
                                                                                count_dofs.offsets_DOFs,
                                                                                count_dofs.offsets_Strongs,
                                                                                global_solution,
                                                                                global_neumann_solution);

        const Eigen::VectorXd Dmatrix_2 =
            dofs_values_2.segment(0, local_count_dofs_2.num_total_dofs - num_local_dofs_pressure);
        const Eigen::VectorXd Dmatrix_pressure_2 =
            dofs_values_2.segment(local_count_dofs_2.num_total_dofs - num_local_dofs_pressure, num_local_dofs_pressure);

        Eigen::MatrixXd global_matrix = Eigen::MatrixXd::Zero(36, 36);
        Eigen::MatrixXd neumann_matrix = Eigen::MatrixXd::Zero(36, 0);
        Eigen::VectorXd global_rhs = Eigen::VectorXd::Zero(36);
        Eigen::VectorXd local_solution = Eigen::VectorXd::Zero(36);
        Eigen::VectorXd local_neuman_solution = Eigen::VectorXd::Zero(0);


        // lato 0
        local_solution.segment(0, 2) = Dmatrix_0.segment(0, 2);
        global_rhs.segment(0, 2) = -local_dirichlet[0][0].middleRows(0, 2);

        // lato 1
        local_solution.segment(2, 2) = Dmatrix_0.segment(2, 2);
        global_rhs.segment(2, 2) += -local_dirichlet[0][1].middleRows(2, 2);
        global_rhs.segment(2, 2) += -local_dirichlet[2][0].middleRows(0, 2);

        // lato 2
        local_solution.segment(4, 2) = Dmatrix_0.segment(4, 2);
        global_rhs.segment(4, 2) = -local_dirichlet[0][2].middleRows(4, 2);

        // // lato 3
        // local_solution.segment(6, 2) = Dmatrix_0.segment(0, 2);

        // // lato 4
        // local_solution.segment(8, 2) = Dmatrix_0.segment(0, 2);

        // lato 5
        local_solution.segment(10, 2) = Dmatrix_2.segment(4, 2);
        global_rhs.segment(10, 2) = -local_dirichlet[2][2].middleRows(4, 2);

        // lato 6
        local_solution.segment(12, 2) = Dmatrix_2.segment(2, 2);
        global_rhs.segment(12, 2) = -local_dirichlet[2][1].middleRows(2, 2);

        // // lato 7
        // local_solution.segment(14, 2) = Dmatrix_0.segment(0, 2);

        // interno 0
        local_solution.segment(16, 2) = Dmatrix_0.segment(6, 2);

        // // interno 1
        // local_solution.segment(18, 2) = Dmatrix_1.segment(6, 2);

        // interno 2
        local_solution.segment(20, 2) = Dmatrix_2.segment(6, 2);

        // // interno 3
        // local_solution.segment(22, 2) = Dmatrix_3.segment(6, 2);

        // pressione 0
        local_solution.segment(24, 3) = Dmatrix_pressure_0;
        global_rhs.segment(24, 3) = local_rhs[0];

        // pressione 1
        // local_solution.segment(27, 3) = Dmatrix_pressure_1;
        // global_rhs.segment(27, 3) = local_rhs[1];

        // pressione 2
        local_solution.segment(30, 3) = Dmatrix_pressure_2;
        global_rhs.segment(30, 3) = local_rhs[2];

        // pressione 3
        // local_solution.segment(33, 3) = Dmatrix_pressure_3;
        // global_rhs.segment(33, 3) = local_rhs[3];


        // edge 0,1,2 lati vs edge 0,1,2
        global_matrix.block(0, 0, 6, 6) += local_A[0].block(0, 0, 6, 6);

        // edge 0,1,2 versus interno 0
        global_matrix.block(0, 16, 6, 2) += local_A[0].block(0, 6, 6, 2);

        // interno 0 versus primi edge 0,1,2
        global_matrix.block(16, 0, 2, 6) += local_A[0].block(6, 0, 2, 6);

        // interno 0 versus interno 0
        global_matrix.block(16, 16, 2, 2) += local_A[0].block(6, 6, 2, 2);

        // edge 0,1,2 versus pression 0
        global_matrix.block(0, 24, 6, 3) += - local_B[0].middleCols(0, 6).transpose();

        // interno 0 versus pression 0
        global_matrix.block(16, 24, 2, 3) += - local_B[0].middleCols(6, 2).transpose();

        // pressione  0 versus edge 0,1,2
        global_matrix.block(24, 0, 3, 6) += local_B[0].middleCols(0, 6);

        // pressione  0 versus interno 0
        global_matrix.block(24, 16, 3, 2) += local_B[0].middleCols(6, 2);

        // pressione 0 versus pressione 0
        global_matrix.block(24, 24, 3, 3) += local_M[0];


        // edge 1 lati vs edge 1
        global_matrix.block(2, 2, 2, 2) += local_A[2].block(0, 0, 2, 2);

        // edge 1 lati vs edge 5
        global_matrix.block(2, 10, 2, 2) += local_A[2].block(0, 4, 2, 2);

        // edge 1 lati vs edge 6
        global_matrix.block(2, 12, 2, 2) += local_A[2].block(0, 2, 2, 2);

        // edge 5 lati vs edge 1
        global_matrix.block(10, 2, 2, 2) += local_A[2].block(4, 0, 2, 2);

        // edge 5 lati vs edge 5
        global_matrix.block(10, 10, 2, 2) += local_A[2].block(4, 4, 2, 2);

        // edge 5 lati vs edge 6
        global_matrix.block(10, 12, 2, 2) += local_A[2].block(4, 2, 2, 2);

        // edge 6 lati vs edge 1
        global_matrix.block(12, 2, 2, 2) += local_A[2].block(2, 0, 2, 2);

        // edge 6 lati vs edge 5
        global_matrix.block(12, 10, 2, 2) += local_A[2].block(2, 4, 2, 2);

        // edge 6 lati vs edge 6
        global_matrix.block(12, 12, 2, 2) += local_A[2].block(2, 2, 2, 2);

        // edge 1 versus interno 0
        global_matrix.block(2, 20, 2, 2) += local_A[2].block(0, 6, 2, 2);

        // edge 5 versus interno 0
        global_matrix.block(10, 20, 2, 2) += local_A[2].block(4, 6, 2, 2);

        // edge 6 versus interno 0
        global_matrix.block(12, 20, 2, 2) += local_A[2].block(2, 6, 2, 2);

        // interno 0 versus primi edge 1
        global_matrix.block(20, 2, 2, 2) += local_A[2].block(6, 0, 2, 2);

        // interno 0 versus primi edge 5
        global_matrix.block(20, 10, 2, 2) += local_A[2].block(6, 4, 2, 2);

        // interno 0 versus primi edge 6
        global_matrix.block(20, 12, 2, 2) += local_A[2].block(6, 2, 2, 2);

        // interno 2 versus interno 2
        global_matrix.block(20, 20, 2, 2) += local_A[2].block(6, 6, 2, 2);

        // edge 1 versus pression 2
        global_matrix.block(2, 30, 2, 3) += - local_B[2].middleCols(0, 2).transpose();

        // edge 5 versus pression 2
        global_matrix.block(10, 30, 2, 3) += - local_B[2].middleCols(4, 2).transpose();

        // edge 6 versus pression 2
        global_matrix.block(12, 30, 2, 3) += - local_B[2].middleCols(2, 2).transpose();

        // interno 2 versus pression 2
        global_matrix.block(20, 30, 2, 3) += - local_B[2].middleCols(6, 2).transpose();

        // pressione 2 versus edge 1
        global_matrix.block(30, 2, 3, 2) += local_B[2].middleCols(0, 2);

        // pressione 2 versus edge 5
        global_matrix.block(30, 10, 3, 2) += local_B[2].middleCols(4, 2);

        // pressione 2 versus edge 6
        global_matrix.block(30, 12, 3, 2) += local_B[2].middleCols(2, 2);

        // pressione 2 versus interno 2
        global_matrix.block(30, 20, 3, 2) += local_B[2].middleCols(6, 2);

        // pressione 2 versus pressione 2
        global_matrix.block(30, 30, 3, 3) += local_M[2];



        Eigen::VectorXd error = global_matrix * local_solution - (global_rhs - neumann_matrix * local_neuman_solution);
        const double max_error = error.cwiseAbs().maxCoeff();

        std::cout << std::scientific;
        std::cout.precision(16);
        std::cout << "totale: " << 0 << " " << max_error << std::endl;

        ASSERT_TRUE(max_error < 1.0e-10);
    }

    {

        const unsigned int num_local_dofs_pressure = dofs_data[1].CellsGlobalDOFs[2].at(2).size();

        const auto local_count_dofs_2 = Polydim::PDETools::Assembler_Utilities::local_count_dofs<2>(2, dofs_data);


        const Eigen::VectorXd dofs_values_2 =
            PDETools::Assembler_Utilities::global_solution_to_local_solution<2>(2,
                                                                                dofs_data,
                                                                                local_count_dofs_2.num_total_dofs,
                                                                                local_count_dofs_2.offsets_DOFs,
                                                                                count_dofs.offsets_DOFs,
                                                                                count_dofs.offsets_Strongs,
                                                                                global_solution,
                                                                                global_neumann_solution);

        const Eigen::VectorXd Dmatrix_2 =
            dofs_values_2.segment(0, local_count_dofs_2.num_total_dofs - num_local_dofs_pressure);
        const Eigen::VectorXd Dmatrix_pressure_2 =
            dofs_values_2.segment(local_count_dofs_2.num_total_dofs - num_local_dofs_pressure, num_local_dofs_pressure);

        Eigen::MatrixXd global_matrix_2 = Eigen::MatrixXd::Zero(11, 11);
        Eigen::MatrixXd neumann_matrix_2 = Eigen::MatrixXd::Zero(11, 0);
        Eigen::VectorXd global_rhs_2 = Eigen::VectorXd::Zero(11);
        Eigen::VectorXd local_solution_2 = Eigen::VectorXd::Zero(11);
        Eigen::VectorXd local_neuman_solution_2 = Eigen::VectorXd::Zero(0);


        local_solution_2 << Dmatrix_2, Dmatrix_pressure_2;

        global_rhs_2.bottomRows(3) = local_rhs[2];
        global_rhs_2.segment(0, 2) = -local_dirichlet[2][0].middleRows(0, 2);
        global_rhs_2.segment(2, 2) = -local_dirichlet[2][1].middleRows(2, 2);
        global_rhs_2.segment(4, 2) = -local_dirichlet[2][2].middleRows(4, 2);



        global_matrix_2.block(0, 0, 8, 8) = local_A[2];
        global_matrix_2.block(0, 8, 8, 3) = - local_B[2].transpose();


        global_matrix_2.block(8, 0, 3, 8) = local_B[2];


        global_matrix_2.block(8, 8, 3, 3) = local_M[2];

        Eigen::VectorXd error_2 = global_matrix_2 * local_solution_2 - (global_rhs_2 - neumann_matrix_2 * local_neuman_solution_2);
        const double max_error_2 = error_2.cwiseAbs().maxCoeff();

        std::cout << std::scientific;
        std::cout.precision(16);
        std::cout << "c: " << 2 << " " << max_error_2 << std::endl;

        ASSERT_TRUE(max_error_2 < 1.0e-10);
    }

    {

        const unsigned int num_local_dofs_pressure = dofs_data[1].CellsGlobalDOFs[2].at(0).size();

        const auto local_count_dofs_0 = Polydim::PDETools::Assembler_Utilities::local_count_dofs<2>(0, dofs_data);


        const Eigen::VectorXd dofs_values_0=
            PDETools::Assembler_Utilities::global_solution_to_local_solution<2>(0,
                                                                                dofs_data,
                                                                                local_count_dofs_0.num_total_dofs,
                                                                                local_count_dofs_0.offsets_DOFs,
                                                                                count_dofs.offsets_DOFs,
                                                                                count_dofs.offsets_Strongs,
                                                                                global_solution,
                                                                                global_neumann_solution);

        const Eigen::VectorXd Dmatrix_0 =
            dofs_values_0.segment(0, local_count_dofs_0.num_total_dofs - num_local_dofs_pressure);
        const Eigen::VectorXd Dmatrix_pressure_0 =
            dofs_values_0.segment(local_count_dofs_0.num_total_dofs - num_local_dofs_pressure, num_local_dofs_pressure);

        Eigen::MatrixXd global_matrix_0 = Eigen::MatrixXd::Zero(11, 11);
        Eigen::MatrixXd neumann_matrix_0 = Eigen::MatrixXd::Zero(11, 0);
        Eigen::VectorXd global_rhs_0 = Eigen::VectorXd::Zero(11);
        Eigen::VectorXd local_solution_0 = Eigen::VectorXd::Zero(11);
        Eigen::VectorXd local_neuman_solution_0 = Eigen::VectorXd::Zero(0);


        local_solution_0 << Dmatrix_0, Dmatrix_pressure_0;

        global_rhs_0.bottomRows(3) = local_rhs[0];
        global_rhs_0.segment(0, 2) = -local_dirichlet[0][0].middleRows(0, 2);
        global_rhs_0.segment(2, 2) = -local_dirichlet[0][1].middleRows(2, 2);
        global_rhs_0.segment(4, 2) = -local_dirichlet[0][2].middleRows(4, 2);

        global_matrix_0.block(0, 0, 8, 8) = local_A[0];
        global_matrix_0.block(0, 8, 8, 3) = - local_B[0].transpose();

        global_matrix_0.block(8, 0, 3, 8) = local_B[0];

        global_matrix_0.block(8, 8, 3, 3) = local_M[0];

        Eigen::VectorXd error_0 = global_matrix_0 * local_solution_0 - (global_rhs_0 - neumann_matrix_0 * local_neuman_solution_0);
        const double max_error_0 = error_0.cwiseAbs().maxCoeff();

        std::cout << std::scientific;
        std::cout.precision(16);
        std::cout << "c: " << 0 << " " << max_error_0 << std::endl;

        ASSERT_TRUE(max_error_0 < 1.0e-10);
    }

    {

        const unsigned int num_local_dofs_pressure = dofs_data[1].CellsGlobalDOFs[2].at(1).size();

        const auto local_count_dofs_1 = Polydim::PDETools::Assembler_Utilities::local_count_dofs<2>(1, dofs_data);


        const Eigen::VectorXd dofs_values_1 =
            PDETools::Assembler_Utilities::global_solution_to_local_solution<2>(1,
                                                                                dofs_data,
                                                                                local_count_dofs_1.num_total_dofs,
                                                                                local_count_dofs_1.offsets_DOFs,
                                                                                count_dofs.offsets_DOFs,
                                                                                count_dofs.offsets_Strongs,
                                                                                global_solution,
                                                                                global_neumann_solution);

        const Eigen::VectorXd Dmatrix_1 =
            dofs_values_1.segment(0, local_count_dofs_1.num_total_dofs - num_local_dofs_pressure);
        const Eigen::VectorXd Dmatrix_pressure_1 =
            dofs_values_1.segment(local_count_dofs_1.num_total_dofs - num_local_dofs_pressure, num_local_dofs_pressure);

        Eigen::MatrixXd global_matrix_1 = Eigen::MatrixXd::Zero(5, 5);
        Eigen::MatrixXd neumann_matrix_1 = Eigen::MatrixXd::Zero(5, 6);
        Eigen::VectorXd global_rhs_1 = Eigen::VectorXd::Zero(5);
        Eigen::VectorXd local_solution_1 = Eigen::VectorXd::Zero(5);
        Eigen::VectorXd local_neuman_solution_1 = Eigen::VectorXd::Zero(6);

        local_neuman_solution_1 = Dmatrix_1.topRows(6);
        local_solution_1 << Dmatrix_1.bottomRows(2), Dmatrix_pressure_1;

        global_rhs_1.bottomRows(3) = local_rhs[1];


        neumann_matrix_1.block(0, 0, 2, 6) = local_A[1].block(6, 0, 2, 6);
        neumann_matrix_1.block(2, 0, 3, 6) = local_B[1].leftCols(6);

        global_matrix_1.block(0, 0, 2, 2) = local_A[1].block(6, 6, 2, 2);
        global_matrix_1.block(0, 2, 2, 3) = - local_B[1].rightCols(2).transpose();
        global_matrix_1.block(2, 0, 3, 2) = local_B[1].rightCols(2);
        global_matrix_1.block(2, 2, 3, 3) = local_M[1];

        const double max_error_1 = (global_matrix_1 * local_solution_1 - (global_rhs_1 - neumann_matrix_1 * local_neuman_solution_1)).cwiseAbs().maxCoeff();

        std::cout << std::scientific;
        std::cout.precision(16);
        std::cout << "c: " << 1 << " " << max_error_1 << std::endl;

        ASSERT_TRUE(max_error_1 < 1.0e-10);
    }

    {

        const unsigned int num_local_dofs_pressure = dofs_data[1].CellsGlobalDOFs[2].at(3).size();

        const auto local_count_dofs_3 = Polydim::PDETools::Assembler_Utilities::local_count_dofs<2>(3, dofs_data);


        const Eigen::VectorXd dofs_values_3 =
            PDETools::Assembler_Utilities::global_solution_to_local_solution<2>(3,
                                                                                dofs_data,
                                                                                local_count_dofs_3.num_total_dofs,
                                                                                local_count_dofs_3.offsets_DOFs,
                                                                                count_dofs.offsets_DOFs,
                                                                                count_dofs.offsets_Strongs,
                                                                                global_solution,
                                                                                global_neumann_solution);

        const Eigen::VectorXd Dmatrix_3 =
            dofs_values_3.segment(0, local_count_dofs_3.num_total_dofs - num_local_dofs_pressure);
        const Eigen::VectorXd Dmatrix_pressure_3 =
            dofs_values_3.segment(local_count_dofs_3.num_total_dofs - num_local_dofs_pressure, num_local_dofs_pressure);

        Eigen::MatrixXd global_matrix_3 = Eigen::MatrixXd::Zero(5, 5);
        Eigen::MatrixXd neumann_matrix_3 = Eigen::MatrixXd::Zero(5, 6);
        Eigen::VectorXd global_rhs_3 = Eigen::VectorXd::Zero(5);
        Eigen::VectorXd local_solution_3 = Eigen::VectorXd::Zero(5);
        Eigen::VectorXd local_neuman_solution_3 = Eigen::VectorXd::Zero(6);

        local_neuman_solution_3 = Dmatrix_3.topRows(6);
        local_solution_3 << Dmatrix_3.bottomRows(2), Dmatrix_pressure_3;

        global_rhs_3.bottomRows(3) = local_rhs[3];


        neumann_matrix_3.block(0, 0, 2, 6) = local_A[3].block(6, 0, 2, 6);
        neumann_matrix_3.block(2, 0, 3, 6) = local_B[3].leftCols(6);

        global_matrix_3.block(0, 0, 2, 2) = local_A[3].block(6, 6, 2, 2);
        global_matrix_3.block(0, 2, 2, 3) = - local_B[3].rightCols(2).transpose();
        global_matrix_3.block(2, 0, 3, 2) = local_B[3].rightCols(2);
        global_matrix_3.block(2, 2, 3, 3) = local_M[3];

        const double max_error_3 = (global_matrix_3 * local_solution_3 - (global_rhs_3 - neumann_matrix_3 * local_neuman_solution_3)).cwiseAbs().maxCoeff();

        std::cout << std::scientific;
        std::cout.precision(16);
        std::cout << "c: " << 3 << " " << max_error_3 << std::endl;

        ASSERT_TRUE(max_error_3 < 1.0e-10);
    }

}

} // namespace UnitTesting
} // namespace Polydim

#endif
