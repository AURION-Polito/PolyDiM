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
#include "Quadrature_Gauss2D_Triangle.hpp"
#include "VEM_Quadrature_2D.hpp"

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
    std::vector<std::array<unsigned int, 2>> DofTypes;

    std::map<std::pair<unsigned int, unsigned int>, std::pair<unsigned int, bool>> Edges_by_vertices;

    Gedim::Quadrature::QuadratureData ReferenceTriangleQuadrature;
    Gedim::Quadrature::QuadratureData ReferenceSquareQuadrature;

    Eigen::MatrixXd ReferenceBasisFunctionValues;
    std::vector<Eigen::MatrixXd> ReferenceBasisFunctionDerivativeValues;
    std::array<Eigen::MatrixXd, 4> ReferenceBasisFunctionSecondDerivativeValues;

    FEM_PCC_1D_ReferenceElement_Data BoundaryReferenceElement_Data;
};

struct FEM_Quadrilateral_PCC_2D_ReferenceElement final
{

    Eigen::MatrixXd Vertices;

    FEM_Quadrilateral_PCC_2D_ReferenceElement()
    {
        Vertices = Eigen::MatrixXd::Zero(3, 4);
        Vertices << 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0;
    }
    ~FEM_Quadrilateral_PCC_2D_ReferenceElement() {};

    FEM_Quadrilateral_PCC_2D_ReferenceElement_Data Create(const unsigned int order) const
    {

        if (order <= 0)
            throw std::runtime_error("not valid order");

        FEM_Quadrilateral_PCC_2D_ReferenceElement_Data result;

        result.Dimension = 2;
        result.Order = order;
        result.NumBasisFunctions = (order + 1) * (order + 1);
        result.NumDofs0D = 1;
        result.NumDofs1D = order - 1;
        result.NumDofs2D = (order - 1) * (order - 1);

        FEM_PCC_1D_ReferenceElement boundary_reference_element;
        result.BoundaryReferenceElement_Data = boundary_reference_element.Create(order, FEM_PCC_1D_Types::Equispaced);
        const Eigen::VectorXd reference_edge_dofs_poisitions = result.BoundaryReferenceElement_Data.DofPositions.row(0);

        // Edge directions
        result.Edges_by_vertices = {{{0, 1}, {0, true}},
                                    {{1, 0}, {0, false}},
                                    {{1, 2}, {1, true}},
                                    {{2, 1}, {1, false}},
                                    {{2, 3}, {2, false}},
                                    {{3, 2}, {2, true}},
                                    {{3, 0}, {3, false}},
                                    {{0, 3}, {3, true}}};

        // Reordering Dofs using convention [point, edge, cell]
        result.DofPositions.setZero(3, result.NumBasisFunctions);
        result.DofTypes.resize(result.NumBasisFunctions);

        for (unsigned int d1 = 0; d1 < 2; d1++)
        {
            if (d1 == 0)
            {
                for (unsigned int d2 = 0; d2 < 2; d2++)
                {
                    result.DofPositions.col(2 * d1 + d2) << reference_edge_dofs_poisitions(d2),
                        reference_edge_dofs_poisitions(d1), 0.0;
                    result.DofTypes[2 * d1 + d2] = {d2, d1};
                }
            }
            else
            {
                for (unsigned int d2 = 0; d2 < 2; d2++)
                {
                    const unsigned int index = (d2 == 0) ? 1 : 0;
                    result.DofPositions.col(2 * d1 + d2) << reference_edge_dofs_poisitions(index),
                        reference_edge_dofs_poisitions(d1), 0.0;
                    result.DofTypes[2 * d1 + d2] = {index, d1};
                }
            }
        }

        if (order > 1)
        {

            for (unsigned int d1 = 0; d1 < 2; d1++)
            {
                result.DofPositions.row(0).segment(4 * result.NumDofs0D + 2 * d1 * result.NumDofs1D, result.NumDofs1D) =
                    reference_edge_dofs_poisitions.segment(2, result.NumDofs1D);
                result.DofPositions.row(1).segment(4 * result.NumDofs0D + 2 * d1 * result.NumDofs1D, result.NumDofs1D) =
                    Eigen::VectorXd::Ones(result.NumDofs1D) * reference_edge_dofs_poisitions(d1);

                for (unsigned int d2 = 2; d2 < result.BoundaryReferenceElement_Data.NumBasisFunctions; d2++)
                    result.DofTypes[4 * result.NumDofs0D + 2 * d1 * result.NumDofs1D + d2 - 2] = {d2, d1};
            }

            for (int d1 = 0; d1 < 2; d1++)
            {
                const unsigned index = (d1 == 1) ? 0 : 1;
                result.DofPositions.row(0).segment(4 * result.NumDofs0D + (2 * d1 + 1) * result.NumDofs1D, result.NumDofs1D) =
                    Eigen::VectorXd::Ones(result.NumDofs1D) * reference_edge_dofs_poisitions(index);
                result.DofPositions.row(1).segment(4 * result.NumDofs0D + (2 * d1 + 1) * result.NumDofs1D, result.NumDofs1D) =
                    reference_edge_dofs_poisitions.segment(2, result.NumDofs1D);

                for (unsigned int d2 = 2; d2 < result.BoundaryReferenceElement_Data.NumBasisFunctions; d2++)
                    result.DofTypes[4 * result.NumDofs0D + (2 * d1 + 1) * result.NumDofs1D + d2 - 2] = {index, d2};
            }

            unsigned int dof = 4 + 4 * result.NumDofs1D;
            for (unsigned int d1 = 2; d1 < result.BoundaryReferenceElement_Data.NumBasisFunctions; d1++)
            {
                for (unsigned int d2 = 2; d2 < result.BoundaryReferenceElement_Data.NumBasisFunctions; d2++)
                {
                    result.DofPositions.col(dof) << reference_edge_dofs_poisitions(d2), reference_edge_dofs_poisitions(d1), 0.0;

                    result.DofTypes[dof] = {d2, d1};

                    dof++;
                }
            }
        }

        std::vector<Eigen::Matrix3d> polygonTriangulationVertices(2);
        polygonTriangulationVertices[0] = Vertices.leftCols(3);
        polygonTriangulationVertices[1] << Vertices.rightCols(2), Vertices.col(0);
        result.ReferenceTriangleQuadrature =
            Gedim::Quadrature::Quadrature_Gauss2D_Triangle::FillPointsAndWeights(2 * (order + 1));
        VEM::Quadrature::VEM_Quadrature_2D quadrature;
        result.ReferenceSquareQuadrature =
            quadrature.PolygonInternalQuadrature(result.ReferenceTriangleQuadrature, polygonTriangulationVertices);

        result.ReferenceBasisFunctionValues = EvaluateBasisFunctions(result.ReferenceSquareQuadrature.Points, result);
        result.ReferenceBasisFunctionDerivativeValues =
            EvaluateBasisFunctionDerivatives(result.ReferenceSquareQuadrature.Points, result);
        //        result.ReferenceBasisFunctionSecondDerivativeValues =
        //            EvaluateBasisFunctionSecondDerivatives(result.ReferenceSquareQuadrature.Points, result);

        return result;
    }

    // ***************************************************************************
    Eigen::MatrixXd EvaluateBasisFunctions(const Eigen::MatrixXd &points,
                                           const FEM_Quadrilateral_PCC_2D_ReferenceElement_Data &reference_element_data) const
    {

        const unsigned int num_points = points.cols();
        Eigen::MatrixXd x = Eigen::MatrixXd::Zero(3, num_points);
        x.row(0) = points.row(0);
        Eigen::MatrixXd y = Eigen::MatrixXd::Zero(3, num_points);
        y.row(0) = points.row(1);

        FEM_PCC_1D_ReferenceElement boundary_reference_element;
        const Eigen::MatrixXd values_x =
            boundary_reference_element.EvaluateBasisFunctions(x, reference_element_data.BoundaryReferenceElement_Data);
        const Eigen::MatrixXd values_y =
            boundary_reference_element.EvaluateBasisFunctions(y, reference_element_data.BoundaryReferenceElement_Data);

        Eigen::MatrixXd values = Eigen::MatrixXd::Ones(num_points, reference_element_data.NumBasisFunctions);

        for (unsigned int d = 0; d < reference_element_data.NumBasisFunctions; d++)
        {
            const auto &dofType = reference_element_data.DofTypes[d];
            values.col(d) = values_x.col(dofType[0]).array() * values_y.col(dofType[1]).array();
        }

        return values;
    }
    // ***************************************************************************
    std::vector<Eigen::MatrixXd> EvaluateBasisFunctionDerivatives(const Eigen::MatrixXd &points,
                                                                  const FEM_Quadrilateral_PCC_2D_ReferenceElement_Data &reference_element_data) const
    {
        const unsigned int num_points = points.cols();
        Eigen::MatrixXd x = Eigen::MatrixXd::Zero(3, num_points);
        x.row(0) = points.row(0);
        Eigen::MatrixXd y = Eigen::MatrixXd::Zero(3, num_points);
        y.row(0) = points.row(1);

        FEM_PCC_1D_ReferenceElement boundary_reference_element;
        const Eigen::MatrixXd values_x =
            boundary_reference_element.EvaluateBasisFunctions(x, reference_element_data.BoundaryReferenceElement_Data);
        const std::vector<Eigen::MatrixXd> values_x_dx =
            boundary_reference_element.EvaluateBasisFunctionDerivatives(x, reference_element_data.BoundaryReferenceElement_Data);
        const Eigen::MatrixXd values_y =
            boundary_reference_element.EvaluateBasisFunctions(y, reference_element_data.BoundaryReferenceElement_Data);
        const std::vector<Eigen::MatrixXd> values_y_dy =
            boundary_reference_element.EvaluateBasisFunctionDerivatives(y, reference_element_data.BoundaryReferenceElement_Data);

        std::vector<Eigen::MatrixXd> grad_values(reference_element_data.Dimension,
                                                 Eigen::MatrixXd::Ones(num_points, reference_element_data.NumBasisFunctions));

        for (unsigned int d = 0; d < reference_element_data.NumBasisFunctions; d++)
        {
            const auto &dofType = reference_element_data.DofTypes[d];
            grad_values[0].col(d) = values_x_dx[0].col(dofType[0]).array() * values_y.col(dofType[1]).array();
            grad_values[1].col(d) = values_x.col(dofType[0]).array() * values_y_dy[0].col(dofType[1]).array();
        }

        return grad_values;
    }
    // ***************************************************************************
    std::array<Eigen::MatrixXd, 4> EvaluateBasisFunctionSecondDerivatives(const Eigen::MatrixXd &,
                                                                          const FEM_Quadrilateral_PCC_2D_ReferenceElement_Data &) const
    {
        throw std::runtime_error("not implemented method");
    }
};
} // namespace PCC
} // namespace FEM
} // namespace Polydim

#endif
