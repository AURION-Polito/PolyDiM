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

#include "FEM_MCC_2D_LocalSpace.hpp"
#include "FEM_Triangle_RT_MCC_2D_LocalSpace.hpp"
#include "FEM_Triangle_RT_MCC_2D_ReferenceElement.hpp"
#include "GeometryUtilities.hpp"
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
            reference_element.EvaluateVelociytBasisFunctionsDivergence(referenceQuadraturePoints, reference_element_data);

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

        const Eigen::MatrixXd ConsistencyErrorX =
            monomials_values -
            velocityBasisValues[0] * Dmatrix.leftCols(reference_element_data.rt_triangle_reference_element_data.Nk);
        const Eigen::MatrixXd ConsistencyErrorY =
            monomials_values -
            velocityBasisValues[1] * Dmatrix.rightCols(reference_element_data.rt_triangle_reference_element_data.Nk);

        const double max_error_x = ConsistencyErrorX.cwiseAbs().maxCoeff();
        const double max_error_y = ConsistencyErrorY.cwiseAbs().maxCoeff();

        ASSERT_TRUE(max_error_x < 1.0e-10);
        ASSERT_TRUE(max_error_y < 1.0e-10);
    }
}

} // namespace UnitTesting
} // namespace Polydim

#endif
