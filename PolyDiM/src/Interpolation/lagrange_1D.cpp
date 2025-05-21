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

#include "lagrange_1D.hpp"

namespace Polydim
{
namespace Interpolation
{
namespace Lagrange
{
//****************************************************************************
Eigen::VectorXd Lagrange_1D_coefficients(const Eigen::VectorXd &interpolation_points_x)
{
    const unsigned int num_interpolation_points = interpolation_points_x.size();

    if (num_interpolation_points == 0)
        return Eigen::VectorXd(0);

    // vector below is defined as differences[i*(num_points-1)+j] = x_j - x_i ,
    // where x_k is the k-th interpolation point, for any j != i
    Eigen::MatrixXd differences(num_interpolation_points, num_interpolation_points - 1);
    for (unsigned int i = 0; i < num_interpolation_points; ++i)
    {
        unsigned int col = 0;
        for (unsigned int j = 0; j < num_interpolation_points; ++j)
        {
            if (j != i)
            {
                differences(i, col) = interpolation_points_x(i) - interpolation_points_x(j);
                col++;
            }
        }
    }
    // Compute coefficients[i]:
    // - compute prod_{j != i} (x_j - x_i)
    Eigen::VectorXd coefficients = differences.rowwise().prod();
    for (unsigned int i = 0; i < num_interpolation_points; ++i)
    {
        // - invert result.
        coefficients[i] = 1.0 / coefficients[i];
    }

    return coefficients;
}
//****************************************************************************
Eigen::MatrixXd Lagrange_1D_values(const Eigen::VectorXd &interpolation_points_x,
                                   const Eigen::VectorXd &lagrange_1D_coefficients,
                                   const Eigen::VectorXd &evaluation_points_x)
{
    const unsigned int num_interpolation_points = interpolation_points_x.size();
    const unsigned int num_evaluation_points = evaluation_points_x.size();

    if (num_interpolation_points == 1)
        return Eigen::VectorXd::Ones(num_evaluation_points);

    Eigen::MatrixXd differences(num_evaluation_points, num_interpolation_points);
    for (unsigned int i = 0; i < num_evaluation_points; ++i)
        differences.row(i) = evaluation_points_x(i) - interpolation_points_x.array();

    Eigen::MatrixXd values = Eigen::MatrixXd::Zero(num_evaluation_points, num_interpolation_points);

    for (unsigned int i = 0; i < num_interpolation_points; ++i)
    {
        values.col(i).setConstant(lagrange_1D_coefficients[i]);
        for (unsigned int j = 0; j < num_interpolation_points; ++j)
        {
            if (i != j)
                values.col(i) = values.col(i).cwiseProduct(differences.col(j));
        }
    }

    return values;
}
//****************************************************************************
Eigen::MatrixXd Lagrange_1D_derivative_values(const Eigen::VectorXd &interpolation_points_x,
                                              const Eigen::VectorXd &lagrange_1D_coefficients,
                                              const Eigen::VectorXd &evaluation_points_x)
{
    const unsigned int num_interpolation_points = interpolation_points_x.size();
    const unsigned int num_evaluation_points = evaluation_points_x.size();

    if (num_interpolation_points == 1)
        return Eigen::VectorXd::Zero(num_evaluation_points);

    Eigen::MatrixXd differences(num_evaluation_points, num_interpolation_points);
    for (unsigned int i = 0; i < num_evaluation_points; ++i)
        differences.row(i) = evaluation_points_x(i) - interpolation_points_x.array();

    Eigen::MatrixXd values = Eigen::MatrixXd::Zero(num_evaluation_points, num_interpolation_points);

    for (unsigned int i = 0; i < num_interpolation_points; ++i)
    {
        for (unsigned int j = 0; j < num_interpolation_points; ++j)
        {
            if (j == i)
                continue;

            Eigen::VectorXd col_value = Eigen::VectorXd::Constant(num_evaluation_points, lagrange_1D_coefficients[i]);

            for (unsigned int k = 0; k < num_interpolation_points; ++k)
            {
                if (k == j || k == i)
                    continue;

                col_value = col_value.cwiseProduct(differences.col(k));
            }

            values.col(i) += col_value;
        }
    }

    return values;
}
//****************************************************************************
} // namespace Lagrange
} // namespace Interpolation
} // namespace Polydim
