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

#ifndef __FEM_Hexahedron_PCC_3D_ReferenceElement_HPP
#define __FEM_Hexahedron_PCC_3D_ReferenceElement_HPP

#include "Eigen/Eigen"
#include "FEM_Quadrilateral_PCC_2D_ReferenceElement.hpp"
#include "QuadratureData.hpp"
#include "Quadrature_Gauss3D_Hexahedron.hpp"

namespace Polydim
{
namespace FEM
{
namespace PCC
{
struct FEM_Hexahedron_PCC_3D_ReferenceElement_Data final
{
    unsigned int Dimension;
    unsigned int Order;
    unsigned int NumDofs0D;
    unsigned int NumDofs1D;
    unsigned int NumDofs2D;
    unsigned int NumDofs3D;

    std::map<std::pair<unsigned int, unsigned int>, std::pair<unsigned int, bool>> Edges_by_vertices;
    std::map<std::pair<unsigned int, unsigned int>, unsigned int> Faces_by_edges;

    unsigned int NumBasisFunctions;
    Eigen::MatrixXd DofPositions;
    std::vector<std::array<unsigned int, 3>> DofTypes;

    Gedim::Quadrature::QuadratureData ReferenceHexahedronQuadrature;

    Eigen::MatrixXd ReferenceBasisFunctionValues;
    std::vector<Eigen::MatrixXd> ReferenceBasisFunctionDerivativeValues;

    FEM_PCC_1D_ReferenceElement_Data EdgeReferenceElement_Data;
    FEM_Quadrilateral_PCC_2D_ReferenceElement_Data BoundaryReferenceElement_Data;
};

struct FEM_Hexahedron_PCC_3D_ReferenceElement final
{
  public:
    Eigen::MatrixXd Vertices;

    FEM_Hexahedron_PCC_3D_ReferenceElement()
    {
        Vertices = Eigen::MatrixXd::Zero(3, 8);
        Vertices << 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0,
            1.0, 1.0, 1.0, 1.0;
    }

    ~FEM_Hexahedron_PCC_3D_ReferenceElement(){};

    FEM_Hexahedron_PCC_3D_ReferenceElement_Data Create(const unsigned int order) const
    {
        if (order <= 0)
            throw std::runtime_error("not valid order");

        FEM_Hexahedron_PCC_3D_ReferenceElement_Data result;

        result.Order = order;
        result.Dimension = 3;
        result.NumBasisFunctions = (order + 1) * (order + 1) * (order + 1);
        result.NumDofs0D = 1;
        result.NumDofs1D = order - 1;
        result.NumDofs2D = (order - 1) * (order - 1);
        result.NumDofs3D = (order - 1) * (order - 1) * (order - 1);

        // Edge directions
        result.Edges_by_vertices = {
            {{0, 1}, {0, true}},  {{1, 0}, {0, false}},  {{1, 2}, {1, true}},  {{2, 1}, {1, false}},
            {{2, 3}, {2, false}}, {{3, 2}, {2, true}},   {{3, 0}, {3, false}}, {{0, 3}, {3, true}},
            {{4, 5}, {4, true}},  {{5, 4}, {4, false}},  {{5, 6}, {5, true}},  {{6, 5}, {5, false}},
            {{6, 7}, {6, false}}, {{7, 6}, {6, true}},   {{7, 4}, {7, false}}, {{4, 7}, {7, true}},
            {{0, 4}, {8, true}},  {{4, 0}, {8, false}},  {{1, 5}, {9, true}},  {{5, 1}, {9, false}},
            {{3, 7}, {10, true}}, {{7, 3}, {10, false}}, {{2, 6}, {11, true}}, {{6, 2}, {11, false}}};

        result.Faces_by_edges = {{{0, 2}, 0},  {{2, 0}, 0},  {{1, 3}, 0},   {{3, 1}, 0},   {{4, 6}, 1},  {{6, 4}, 1},
                                 {{5, 7}, 1},  {{7, 5}, 1},  {{0, 4}, 2},   {{4, 0}, 2},   {{8, 9}, 2},  {{9, 8}, 2},
                                 {{2, 6}, 3},  {{6, 2}, 3},  {{10, 11}, 3}, {{11, 10}, 3}, {{7, 3}, 4},  {{3, 7}, 4},
                                 {{8, 10}, 4}, {{10, 8}, 4}, {{1, 5}, 5},   {{9, 11}, 5},  {{11, 9}, 5}, {{5, 1}, 5}};

        FEM_PCC_1D_ReferenceElement edge_reference_element;
        result.EdgeReferenceElement_Data = edge_reference_element.Create(order, FEM_PCC_1D_Types::Equispaced);
        const Eigen::VectorXd reference_edge_dofs_poisitions = result.EdgeReferenceElement_Data.DofPositions.row(0);

        // Reordering Dofs using convention [point, edge, cell]
        result.DofPositions.setZero(3, result.NumBasisFunctions);
        result.DofTypes.resize(result.NumBasisFunctions);

        for (unsigned int d3 = 0; d3 < 2; d3++)
        {
            for (unsigned int d1 = 0; d1 < 2; d1++)
            {
                if (d1 == 0)
                {
                    for (unsigned int d2 = 0; d2 < 2; d2++)
                    {
                        result.DofPositions.col(4 * d3 + 2 * d1 + d2) << reference_edge_dofs_poisitions(d2),
                            reference_edge_dofs_poisitions(d1), reference_edge_dofs_poisitions(d3);
                        result.DofTypes[4 * d3 + 2 * d1 + d2] = {d2, d1, d3};
                    }
                }
                else
                {
                    for (unsigned int d2 = 0; d2 < 2; d2++)
                    {
                        const unsigned int index = (d2 == 0) ? 1 : 0;
                        result.DofPositions.col(4 * d3 + 2 * d1 + d2) << reference_edge_dofs_poisitions(index),
                            reference_edge_dofs_poisitions(d1), reference_edge_dofs_poisitions(d3);
                        result.DofTypes[4 * d3 + 2 * d1 + d2] = {index, d1, d3};
                    }
                }
            }
        }

        if (order > 1)
        {
            for (unsigned int d3 = 0; d3 < 2; d3++)
            {
                for (unsigned int d1 = 0; d1 < 2; d1++)
                {
                    result.DofPositions.row(0).segment(8 * result.NumDofs0D + (2 * d1 + 4 * d3) * result.NumDofs1D,
                                                       result.NumDofs1D) =
                        reference_edge_dofs_poisitions.segment(2, result.NumDofs1D);
                    result.DofPositions.row(1).segment(8 * result.NumDofs0D + (2 * d1 + 4 * d3) * result.NumDofs1D,
                                                       result.NumDofs1D) =
                        Eigen::VectorXd::Ones(result.NumDofs1D) * reference_edge_dofs_poisitions(d1);
                    result.DofPositions.row(2).segment(8 * result.NumDofs0D + (2 * d1 + 4 * d3) * result.NumDofs1D,
                                                       result.NumDofs1D) =
                        Eigen::VectorXd::Ones(result.NumDofs1D) * reference_edge_dofs_poisitions(d3);

                    for (unsigned int d2 = 2; d2 < result.EdgeReferenceElement_Data.NumBasisFunctions; d2++)
                        result.DofTypes[8 * result.NumDofs0D + (2 * d1 + 4 * d3) * result.NumDofs1D + d2 - 2] = {d2, d1, d3};
                }
            }

            for (unsigned int d3 = 0; d3 < 2; d3++)
            {
                for (unsigned int d1 = 0; d1 < 2; d1++)
                {
                    const unsigned index = (d1 == 1) ? 0 : 1;
                    result.DofPositions.row(0).segment(8 * result.NumDofs0D + (2 * d1 + 1 + 4 * d3) * result.NumDofs1D,
                                                       result.NumDofs1D) =
                        Eigen::VectorXd::Ones(result.NumDofs1D) * reference_edge_dofs_poisitions(index);
                    result.DofPositions.row(1).segment(8 * result.NumDofs0D + (2 * d1 + 1 + 4 * d3) * result.NumDofs1D,
                                                       result.NumDofs1D) =
                        reference_edge_dofs_poisitions.segment(2, result.NumDofs1D);
                    result.DofPositions.row(2).segment(8 * result.NumDofs0D + (2 * d1 + 1 + 4 * d3) * result.NumDofs1D,
                                                       result.NumDofs1D) =
                        Eigen::VectorXd::Ones(result.NumDofs1D) * reference_edge_dofs_poisitions(d3);

                    for (unsigned int d2 = 2; d2 < result.EdgeReferenceElement_Data.NumBasisFunctions; d2++)
                        result.DofTypes[8 * result.NumDofs0D + (2 * d1 + 1 + 4 * d3) * result.NumDofs1D + d2 - 2] = {index, d2, d3};
                }
            }

            // Verticel edges:
            for (unsigned int d3 = 0; d3 < 2; d3++) // y
            {
                for (unsigned int d1 = 0; d1 < 2; d1++) // x
                {
                    result.DofPositions.row(0).segment(8 * result.NumDofs0D + (8 + d1 + 2 * d3) * result.NumDofs1D,
                                                       result.NumDofs1D) =
                        Eigen::VectorXd::Ones(result.NumDofs1D) * reference_edge_dofs_poisitions(d1);
                    result.DofPositions.row(1).segment(8 * result.NumDofs0D + (8 + d1 + 2 * d3) * result.NumDofs1D,
                                                       result.NumDofs1D) =
                        Eigen::VectorXd::Ones(result.NumDofs1D) * reference_edge_dofs_poisitions(d3);
                    result.DofPositions.row(2).segment(8 * result.NumDofs0D + (8 + d1 + 2 * d3) * result.NumDofs1D,
                                                       result.NumDofs1D) =
                        reference_edge_dofs_poisitions.segment(2, result.NumDofs1D);

                    for (unsigned int d2 = 2; d2 < result.EdgeReferenceElement_Data.NumBasisFunctions; d2++) // z
                        result.DofTypes[8 * result.NumDofs0D + (8 + d1 + 2 * d3) * result.NumDofs1D + d2 - 2] = {d1, d3, d2};
                }
            }

            unsigned int dof = 8 * result.NumDofs0D + 12 * result.NumDofs1D;

            // Face dofs
            for (unsigned int d3 = 0; d3 < 2; d3++) // z
            {
                for (unsigned int d1 = 2; d1 < result.EdgeReferenceElement_Data.NumBasisFunctions; d1++) // x
                {
                    for (unsigned int d2 = 2; d2 < result.EdgeReferenceElement_Data.NumBasisFunctions; d2++) // y
                    {
                        result.DofPositions.col(dof) << reference_edge_dofs_poisitions(d1),
                            reference_edge_dofs_poisitions(d2), reference_edge_dofs_poisitions(d3);
                        result.DofTypes[dof] = {d1, d2, d3};
                        dof++;
                    }
                }
            }

            for (unsigned int d3 = 0; d3 < 2; d3++) // y
            {
                for (unsigned int d1 = 2; d1 < result.EdgeReferenceElement_Data.NumBasisFunctions; d1++) // x
                {
                    for (unsigned int d2 = 2; d2 < result.EdgeReferenceElement_Data.NumBasisFunctions; d2++) // z
                    {
                        result.DofPositions.col(dof) << reference_edge_dofs_poisitions(d1),
                            reference_edge_dofs_poisitions(d3), reference_edge_dofs_poisitions(d2);
                        result.DofTypes[dof] = {d1, d3, d2};
                        dof++;
                    }
                }
            }

            for (unsigned int d3 = 0; d3 < 2; d3++) // x
            {
                for (unsigned int d1 = 2; d1 < result.EdgeReferenceElement_Data.NumBasisFunctions; d1++) // y
                {
                    for (unsigned int d2 = 2; d2 < result.EdgeReferenceElement_Data.NumBasisFunctions; d2++) // z
                    {
                        result.DofPositions.col(dof) << reference_edge_dofs_poisitions(d3),
                            reference_edge_dofs_poisitions(d1), reference_edge_dofs_poisitions(d2);
                        result.DofTypes[dof] = {d3, d1, d2};
                        dof++;
                    }
                }
            }

            // Interni
            for (unsigned int d3 = 2; d3 < result.EdgeReferenceElement_Data.NumBasisFunctions; d3++) // x
            {
                for (unsigned int d1 = 2; d1 < result.EdgeReferenceElement_Data.NumBasisFunctions; d1++) // y
                {
                    for (unsigned int d2 = 2; d2 < result.EdgeReferenceElement_Data.NumBasisFunctions; d2++) // z
                    {
                        result.DofPositions.col(dof) << reference_edge_dofs_poisitions(d3),
                            reference_edge_dofs_poisitions(d1), reference_edge_dofs_poisitions(d2);
                        result.DofTypes[dof] = {d3, d1, d2};
                        dof++;
                    }
                }
            }
        }

        FEM_Quadrilateral_PCC_2D_ReferenceElement boundary_reference_element;
        result.BoundaryReferenceElement_Data = boundary_reference_element.Create(order);

        result.ReferenceHexahedronQuadrature =
            Gedim::Quadrature::Quadrature_Gauss3D_Hexahedron::FillPointsAndWeights(2 * (order + 1));

        result.ReferenceBasisFunctionValues = EvaluateBasisFunctions(result.ReferenceHexahedronQuadrature.Points, result);
        result.ReferenceBasisFunctionDerivativeValues =
            EvaluateBasisFunctionDerivatives(result.ReferenceHexahedronQuadrature.Points, result);

        return result;
    }
    // ***************************************************************************
    Eigen::MatrixXd EvaluateBasisFunctions(const Eigen::MatrixXd &points,
                                           const FEM_Hexahedron_PCC_3D_ReferenceElement_Data &reference_element_data) const
    {

        const unsigned int num_points = points.cols();
        Eigen::MatrixXd x = Eigen::MatrixXd::Zero(3, num_points);
        x.row(0) = points.row(0);
        Eigen::MatrixXd y = Eigen::MatrixXd::Zero(3, num_points);
        y.row(0) = points.row(1);
        Eigen::MatrixXd z = Eigen::MatrixXd::Zero(3, num_points);
        z.row(0) = points.row(2);

        FEM_PCC_1D_ReferenceElement boundary_reference_element;
        const Eigen::MatrixXd values_x =
            boundary_reference_element.EvaluateBasisFunctions(x, reference_element_data.EdgeReferenceElement_Data);
        const Eigen::MatrixXd values_y =
            boundary_reference_element.EvaluateBasisFunctions(y, reference_element_data.EdgeReferenceElement_Data);
        const Eigen::MatrixXd values_z =
            boundary_reference_element.EvaluateBasisFunctions(z, reference_element_data.EdgeReferenceElement_Data);

        Eigen::MatrixXd values = Eigen::MatrixXd::Ones(num_points, reference_element_data.NumBasisFunctions);

        for (unsigned int d = 0; d < reference_element_data.NumBasisFunctions; d++)
        {
            const auto &dofType = reference_element_data.DofTypes[d];
            values.col(d) =
                values_x.col(dofType[0]).array() * values_y.col(dofType[1]).array() * values_z.col(dofType[2]).array();
        }

        return values;
    }
    // ***************************************************************************
    std::vector<Eigen::MatrixXd> EvaluateBasisFunctionDerivatives(const Eigen::MatrixXd &points,
                                                                  const FEM_Hexahedron_PCC_3D_ReferenceElement_Data &reference_element_data) const
    {
        const unsigned int num_points = points.cols();
        Eigen::MatrixXd x = Eigen::MatrixXd::Zero(3, num_points);
        x.row(0) = points.row(0);
        Eigen::MatrixXd y = Eigen::MatrixXd::Zero(3, num_points);
        y.row(0) = points.row(1);
        Eigen::MatrixXd z = Eigen::MatrixXd::Zero(3, num_points);
        z.row(0) = points.row(2);

        FEM_PCC_1D_ReferenceElement boundary_reference_element;
        const Eigen::MatrixXd values_x =
            boundary_reference_element.EvaluateBasisFunctions(x, reference_element_data.EdgeReferenceElement_Data);
        const Eigen::MatrixXd values_y =
            boundary_reference_element.EvaluateBasisFunctions(y, reference_element_data.EdgeReferenceElement_Data);
        const Eigen::MatrixXd values_z =
            boundary_reference_element.EvaluateBasisFunctions(z, reference_element_data.EdgeReferenceElement_Data);

        const std::vector<Eigen::MatrixXd> values_x_dx =
            boundary_reference_element.EvaluateBasisFunctionDerivatives(x, reference_element_data.EdgeReferenceElement_Data);
        const std::vector<Eigen::MatrixXd> values_y_dy =
            boundary_reference_element.EvaluateBasisFunctionDerivatives(y, reference_element_data.EdgeReferenceElement_Data);
        const std::vector<Eigen::MatrixXd> values_z_dz =
            boundary_reference_element.EvaluateBasisFunctionDerivatives(z, reference_element_data.EdgeReferenceElement_Data);

        std::vector<Eigen::MatrixXd> grad_values(3, Eigen::MatrixXd::Ones(num_points, reference_element_data.NumBasisFunctions));

        for (unsigned int d = 0; d < reference_element_data.NumBasisFunctions; d++)
        {
            const auto &dofType = reference_element_data.DofTypes[d];
            grad_values[0].col(d) = values_x_dx[0].col(dofType[0]).array() * values_y.col(dofType[1]).array() *
                                    values_z.col(dofType[2]).array();
            grad_values[1].col(d) = values_x.col(dofType[0]).array() * values_y_dy[0].col(dofType[1]).array() *
                                    values_z.col(dofType[2]).array();
            grad_values[2].col(d) = values_x.col(dofType[0]).array() * values_y.col(dofType[1]).array() *
                                    values_z_dz[0].col(dofType[2]).array();
        }

        return grad_values;
    }
    // ***************************************************************************
    std::array<Eigen::MatrixXd, 9> EvaluateBasisFunctionSecondDerivatives(const Eigen::MatrixXd &,
                                                                          const FEM_Quadrilateral_PCC_2D_ReferenceElement_Data &) const
    {
        throw std::runtime_error("not implemented method");
    }
};
} // namespace PCC
} // namespace FEM
} // namespace Polydim

#endif
