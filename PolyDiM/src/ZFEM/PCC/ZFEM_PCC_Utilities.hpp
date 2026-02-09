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

#ifndef __ZFEM_PCC_Utilities_HPP
#define __ZFEM_PCC_Utilities_HPP

#include "Eigen/Eigen"
#include "FEM_Triangle_PCC_2D_LocalSpace.hpp"
#include "I_ZFEM_PCC_2D_ReferenceElement.hpp"
#include <algorithm>
#include <numeric>

namespace Polydim
{
namespace ZFEM
{
namespace PCC
{

struct ZFEM_PCC_Utilities final
{
    Eigen::VectorXd ComputeEdgeBasisCoefficients(const unsigned int &order, const Eigen::VectorXd &edgeInternalPoints) const
    {
        // Compute basis function coefficients on the generic edge.
        Eigen::VectorXd interpolation_points_x(order + 1);
        interpolation_points_x << 0.0, 1.0, edgeInternalPoints;
        return Interpolation::Lagrange::Lagrange_1D_coefficients(interpolation_points_x);
    }

    Eigen::MatrixXd ComputeValuesOnEdge(const Eigen::RowVectorXd &edgeInternalPoints,
                                        const unsigned int &order,
                                        const Eigen::VectorXd &edgeBasisCoefficients,
                                        const Eigen::VectorXd &pointsCurvilinearCoordinates) const
    {
        Eigen::VectorXd interpolation_points_x(order + 1);
        interpolation_points_x << 0.0, 1.0, edgeInternalPoints.transpose();
        return Interpolation::Lagrange::Lagrange_1D_values(interpolation_points_x, edgeBasisCoefficients, pointsCurvilinearCoordinates);
    }

    inline Eigen::MatrixXd ComputeBasisFunctionsValues(const unsigned int &NumBasisFunctions,
                                                       const unsigned int &NumVirtualBasisFunctions,
                                                       const Eigen::MatrixXd &fem_basis_functions_values,
                                                       const Eigen::MatrixXd &VirtualWeights) const
    {
        return fem_basis_functions_values.leftCols(NumBasisFunctions) +
               fem_basis_functions_values.rightCols(NumVirtualBasisFunctions) * VirtualWeights.transpose();
    }

    inline Eigen::MatrixXd ComputeFEMBasisFunctionsValues(const ZFEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                          const unsigned int &num_basis_functions,
                                                          const unsigned int &num_virtual_basis_functions,
                                                          const std::vector<Polydim::FEM::PCC::FEM_Triangle_PCC_2D_LocalSpace_Data> &fem_local_space_data,
                                                          const Eigen::MatrixXi &local_to_global) const
    {
        FEM::PCC::FEM_Triangle_PCC_2D_LocalSpace fem_local_space;
        const unsigned int num_quadrature =
            reference_element_data.fem_reference_element_data.ReferenceTriangleQuadrature.Points.cols();
        const unsigned int num_triangles = fem_local_space_data.size();

        Eigen::MatrixXd total_fem_basis_functions_values =
            Eigen::MatrixXd::Zero(num_quadrature * num_triangles, num_basis_functions + num_virtual_basis_functions);

        unsigned int offeset_quadrature_points = 0;
        for (unsigned int t = 0; t < num_triangles; t++)
        {
            const Eigen::MatrixXd fem_basis_function_values =
                fem_local_space.ComputeBasisFunctionsValues(reference_element_data.fem_reference_element_data,
                                                            fem_local_space_data[t]);

            for (unsigned int p = 0; p < local_to_global.cols(); p++)
                total_fem_basis_functions_values.block(offeset_quadrature_points, local_to_global(t, p), num_quadrature, 1) =
                    fem_basis_function_values.col(p);

            offeset_quadrature_points += num_quadrature;
        }

        return total_fem_basis_functions_values;
    }

    inline Eigen::MatrixXd ComputeFEMBasisFunctionsValues(const ZFEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                          const unsigned int &num_basis_functions,
                                                          const unsigned int &num_virtual_basis_functions,
                                                          const std::vector<Polydim::FEM::PCC::FEM_Triangle_PCC_2D_LocalSpace_Data> &fem_local_space_data,
                                                          const Eigen::MatrixXi &local_to_global,
                                                          const std::vector<Eigen::MatrixXd> &points) const
    {
        FEM::PCC::FEM_Triangle_PCC_2D_LocalSpace fem_local_space;
        unsigned int num_quadrature = 0;

        for (unsigned int t = 0; t < points.size(); t++)
            num_quadrature += points[t].cols();

        Eigen::MatrixXd total_fem_basis_functions_values =
            Eigen::MatrixXd::Zero(num_quadrature, num_basis_functions + num_virtual_basis_functions);

        unsigned int offeset_quadrature_points = 0;
        for (unsigned int t = 0; t < points.size(); t++)
        {
            if (points[t].size() > 0)
            {

                const unsigned int num_ref_quadrature = points[t].cols();

                const Eigen::MatrixXd fem_basis_function_values =
                    fem_local_space.ComputeBasisFunctionsValues(reference_element_data.fem_reference_element_data,
                                                                fem_local_space_data[t],
                                                                points[t]);

                for (unsigned int p = 0; p < local_to_global.cols(); p++)
                    total_fem_basis_functions_values.block(offeset_quadrature_points, local_to_global(t, p), num_ref_quadrature, 1) =
                        fem_basis_function_values.col(p);

                offeset_quadrature_points += num_ref_quadrature;
            }
        }

        return total_fem_basis_functions_values;
    }

    inline std::vector<Eigen::MatrixXd> ComputeFEMBasisFunctionsDerivativeValues(
        unsigned int dimension,
        const ZFEM_PCC_2D_ReferenceElement_Data &reference_element_data,
        const unsigned int &num_basis_functions,
        const unsigned int &num_virtual_basis_functions,
        const std::vector<Polydim::FEM::PCC::FEM_Triangle_PCC_2D_LocalSpace_Data> &fem_local_space_data,
        const Eigen::MatrixXi &local_to_global) const
    {
        FEM::PCC::FEM_Triangle_PCC_2D_LocalSpace fem_local_space;
        const unsigned int num_quadrature =
            reference_element_data.fem_reference_element_data.ReferenceTriangleQuadrature.Points.cols();
        const unsigned int num_triangles = fem_local_space_data.size();

        std::vector<Eigen::MatrixXd> total_fem_basis_functions_derivative_values(
            dimension,
            Eigen::MatrixXd::Zero(num_quadrature * num_triangles, num_basis_functions + num_virtual_basis_functions));

        unsigned int offeset_quadrature_points = 0;
        for (unsigned int t = 0; t < num_triangles; t++)
        {

            const std::vector<Eigen::MatrixXd> fem_basis_function_derivatives_values =
                fem_local_space.ComputeBasisFunctionsDerivativeValues(reference_element_data.fem_reference_element_data,
                                                                      fem_local_space_data[t]);

            for (unsigned int d = 0; d < dimension; d++)
            {
                for (unsigned int p = 0; p < local_to_global.cols(); p++)
                    total_fem_basis_functions_derivative_values[d].block(offeset_quadrature_points, local_to_global(t, p), num_quadrature, 1) =
                        fem_basis_function_derivatives_values[d].col(p);
            }

            offeset_quadrature_points += num_quadrature;
        }

        return total_fem_basis_functions_derivative_values;
    }

    inline std::vector<Eigen::MatrixXd> ComputeFEMBasisFunctionsDerivativeValues(
        unsigned int dimension,
        const ZFEM_PCC_2D_ReferenceElement_Data &reference_element_data,
        const unsigned int &num_basis_functions,
        const unsigned int &num_virtual_basis_functions,
        const std::vector<Polydim::FEM::PCC::FEM_Triangle_PCC_2D_LocalSpace_Data> &fem_local_space_data,
        const Eigen::MatrixXi &local_to_global,
        const std::vector<Eigen::MatrixXd> &points) const
    {
        FEM::PCC::FEM_Triangle_PCC_2D_LocalSpace fem_local_space;
        unsigned int num_quadrature = 0;

        for (unsigned int t = 0; t < points.size(); t++)
            num_quadrature += points[t].cols();

        std::vector<Eigen::MatrixXd> total_fem_basis_functions_derivative_values(
            dimension,
            Eigen::MatrixXd::Zero(num_quadrature, num_basis_functions + num_virtual_basis_functions));

        unsigned int offeset_quadrature_points = 0;
        for (unsigned int t = 0; t < points.size(); t++)
        {
            if (points[t].size() > 0)
            {

                const unsigned int num_ref_quadrature = points[t].cols();

                const std::vector<Eigen::MatrixXd> fem_basis_function_derivatives_values =
                    fem_local_space.ComputeBasisFunctionsDerivativeValues(reference_element_data.fem_reference_element_data,
                                                                          fem_local_space_data[t],
                                                                          points[t]);

                for (unsigned int d = 0; d < dimension; d++)
                {
                    for (unsigned int p = 0; p < local_to_global.cols(); p++)
                        total_fem_basis_functions_derivative_values[d].block(offeset_quadrature_points,
                                                                             local_to_global(t, p),
                                                                             num_ref_quadrature,
                                                                             1) =
                            fem_basis_function_derivatives_values[d].col(p);
                }

                offeset_quadrature_points += num_ref_quadrature;
            }
        }

        return total_fem_basis_functions_derivative_values;
    }

    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const unsigned int &Dimension,
                                                                              const unsigned int &NumBasisFunctions,
                                                                              const unsigned int &NumVirtualBasisFunctions,
                                                                              const std::vector<Eigen::MatrixXd> &fem_basis_functions_derivative_values,
                                                                              const Eigen::MatrixXd &VirtualWeights) const
    {
        std::vector<Eigen::MatrixXd> ZFEM_basis_functions_derivative_values(Dimension);

        for (unsigned int d = 0; d < Dimension; d++)
        {
            ZFEM_basis_functions_derivative_values[d] =
                fem_basis_functions_derivative_values[d].leftCols(NumBasisFunctions) +
                fem_basis_functions_derivative_values[d].rightCols(NumVirtualBasisFunctions) * VirtualWeights.transpose();
        }

        return ZFEM_basis_functions_derivative_values;
    }

    Eigen::MatrixXi CreateMaps(const unsigned int order,
                               const unsigned int num_vertices,
                               const unsigned int NumDOFs1D,
                               const unsigned int NumDOFs2D,
                               const unsigned int NumBasisFunctions) const
    {
        Eigen::MatrixXi local_to_total(num_vertices, 3 * (1 + NumDOFs1D) + NumDOFs2D);

        const Eigen::ArrayXi edge_reference_id_dofs = Eigen::VectorXi::LinSpaced(NumDOFs1D, 0, NumDOFs1D - 1);

        // Vertex and boundary DOFs
        for (unsigned int t = 0; t < num_vertices; t++)
        {
            const unsigned int next_t = (t + 1) % num_vertices;

            local_to_total(t, 1) = t;
            local_to_total(t, 2) = next_t;
            local_to_total.row(t).segment(3 + NumDOFs1D, NumDOFs1D) = edge_reference_id_dofs + t * NumDOFs1D + num_vertices;
        }

        if (order <= 2)
        {
            local_to_total.col(0) = Eigen::VectorXi::Constant(num_vertices, NumBasisFunctions);

            for (unsigned int t = 0; t < num_vertices; t++)
            {
                const unsigned int next_t = (t + 1) % num_vertices;

                local_to_total.row(t).segment(3, NumDOFs1D) = edge_reference_id_dofs + NumBasisFunctions + 1 + t * NumDOFs1D;
                local_to_total.row(t).segment(3 + 2 * NumDOFs1D, NumDOFs1D) =
                    edge_reference_id_dofs + NumBasisFunctions + 1 + next_t * NumDOFs1D;
            }
        }
        else if (order == 3)
        {
            local_to_total.col(0) = Eigen::VectorXi::Constant(num_vertices, NumBasisFunctions - 1);

            for (unsigned int t = 0; t < num_vertices; t++)
            {
                const unsigned int next_t = (t + 1) % num_vertices;

                local_to_total.row(t).segment(3, NumDOFs1D) = edge_reference_id_dofs + NumBasisFunctions + t * NumDOFs1D;
                local_to_total.row(t).segment(3 + 2 * NumDOFs1D, NumDOFs1D) =
                    edge_reference_id_dofs + NumBasisFunctions + next_t * NumDOFs1D;
                local_to_total(t, 3 + 3 * NumDOFs1D) = NumBasisFunctions + num_vertices * NumDOFs1D + t;
            }
        }
        else
        {
            local_to_total.col(0) = Eigen::VectorXi::Constant(num_vertices, NumBasisFunctions - NumDOFs2D);

            const unsigned int base = (NumDOFs2D - 1) / num_vertices;
            const unsigned int other = (NumDOFs2D - 1) % num_vertices;

            std::vector<unsigned int> num_pt_triangle(num_vertices, 0);

            const unsigned int num_selected_triangles = NumDOFs2D - 1;

            unsigned int count = 0;
            const unsigned int repeat_times = ceil(num_vertices / ((double)num_selected_triangles));
            for (unsigned int i = 0; i < repeat_times; i++)
            {
                for (unsigned int j = 0; j < num_selected_triangles; j++)
                {
                    if (i + j * repeat_times < num_vertices && count < num_vertices)
                    {
                        num_pt_triangle[i + j * repeat_times] = base + (count < other ? 1 : 0);
                        count++;
                    }
                }
            }

            unsigned int offset_internal = NumBasisFunctions + num_vertices * NumDOFs1D;
            unsigned int offset_internal_dof = NumBasisFunctions - NumDOFs2D + 1;

            std::vector<unsigned int> copy_ordered(NumDOFs2D);
            std::iota(copy_ordered.begin(), copy_ordered.end(), 0);

            std::vector<unsigned int> p_mod(NumDOFs2D);
            for (unsigned int t = 0; t < num_vertices; t++)
            {
                const unsigned int next_t = (t + 1) % num_vertices;

                local_to_total.row(t).segment(3, NumDOFs1D) = edge_reference_id_dofs + NumBasisFunctions + t * NumDOFs1D;
                local_to_total.row(t).segment(3 + 2 * NumDOFs1D, NumDOFs1D) =
                    edge_reference_id_dofs + NumBasisFunctions + next_t * NumDOFs1D;

                if (num_pt_triangle[t] > 0)
                {
                    unsigned int count = 0;
                    const unsigned int repeat_times = ceil(((double)NumDOFs2D) / num_pt_triangle[t]);
                    for (unsigned int i = 0; i < repeat_times; i++)
                    {
                        for (unsigned int j = 0; j < num_pt_triangle[t]; j++)
                        {
                            if (i + j * repeat_times < NumDOFs2D && count < NumDOFs2D)
                                p_mod[count++] = (i + j * repeat_times + t + 1) % NumDOFs2D;
                        }
                    }
                }
                else
                    p_mod = copy_ordered;

                for (unsigned int p = 0; p < num_pt_triangle[t]; p++)
                {
                    local_to_total(t, 3 * (NumDOFs1D + 1) + p_mod[p]) = offset_internal_dof;
                    offset_internal_dof += 1;
                }

                for (unsigned int p = num_pt_triangle[t]; p < NumDOFs2D; p++)
                {
                    local_to_total(t, 3 * (NumDOFs1D + 1) + p_mod[p]) = offset_internal;
                    offset_internal += 1;
                }
            }
        }

        return local_to_total;
    }

    void ComputeMinimizerSumOfSquaredWeightsMonomials(const Eigen::MatrixXd &virtual_vertices,
                                                      const unsigned int NumBasisFunctions,
                                                      const Eigen::MatrixXd &Dmatrix,
                                                      const Eigen::MatrixXd &VanderVirtuals,
                                                      Eigen::MatrixXd &weights) const
    {
        const unsigned int num_virtuals = virtual_vertices.cols();

        const Eigen::MatrixXd Mmatrix = Dmatrix.transpose() * Dmatrix;
        const Eigen::LLT<Eigen::MatrixXd> Mmatrix_LLT = Mmatrix.llt();

        weights.resize(NumBasisFunctions, num_virtuals);
        for (unsigned int j = 0; j < num_virtuals; j++)
        {
            const Eigen::VectorXd rhsP = VanderVirtuals.row(j);
            const Eigen::VectorXd y = Mmatrix_LLT.solve(rhsP);
            weights.col(j) = Dmatrix * y;
        }
    }
};
} // namespace PCC
} // namespace ZFEM
} // namespace Polydim

#endif
