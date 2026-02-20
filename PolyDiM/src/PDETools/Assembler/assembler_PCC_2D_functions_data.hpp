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

#ifndef __PDETOOLS_ASSEMBLER_assembler_PCC_2D_functions_data_HPP
#define __PDETOOLS_ASSEMBLER_assembler_PCC_2D_functions_data_HPP

#include "Eigen/Eigen"

namespace Polydim
{
namespace PDETools
{
namespace Assembler_Utilities
{
namespace PCC_2D
{
// ***************************************************************************
struct Sparse_Matrix_Triplet final
{
    unsigned int i;
    unsigned int j;
    double value;
};
// ***************************************************************************
struct Sparse_Matrix_Data final
{
    std::array<unsigned int, 2> size;
    std::vector<unsigned int> rows;
    std::vector<unsigned int> cols;
    std::vector<double> values;
};
// ***************************************************************************
struct Evaluate_Function_On_DOFs_Data final
{
    Eigen::VectorXd function_dofs;
    Eigen::VectorXd function_strong;
};
// ***************************************************************************
struct Variational_Operator final
{
    Polydim::PDETools::Assembler_Utilities::PCC_2D::Sparse_Matrix_Data A;
    Polydim::PDETools::Assembler_Utilities::PCC_2D::Sparse_Matrix_Data A_Strong;
};
// ***************************************************************************
struct Evaluate_Solution_On_Quadrature_Points_Data final
{
    Eigen::MatrixXd quadrature_points;
    Eigen::VectorXd quadrature_weigths;
    Eigen::VectorXd numeric_solution;
    std::array<Eigen::VectorXd, 3> numeric_gradient_solution;
    Eigen::VectorXd exact_solution;
    std::array<Eigen::VectorXd, 3> exact_gradient_solution;
};
// ***************************************************************************
struct Post_Process_Data_Cell0Ds final
{
    Eigen::VectorXd numeric_solution;
    Eigen::VectorXd exact_solution;
    std::array<Eigen::VectorXd, 3> exact_gradient_solution;
};
// ***************************************************************************
struct Post_Process_Data_ErrorL2 final
{
    Eigen::VectorXd cell2Ds_exact_norm_L2;
    Eigen::VectorXd cell2Ds_numeric_norm_L2;
    Eigen::VectorXd cell2Ds_error_L2;
    double error_L2;
    double exact_norm_L2;
    double numeric_norm_L2;
};
// ***************************************************************************
struct Post_Process_Data_ErrorH1 final
{
    Eigen::VectorXd cell2Ds_exact_norm_H1;
    Eigen::VectorXd cell2Ds_numeric_norm_H1;
    Eigen::VectorXd cell2Ds_error_H1;
    double error_H1;
    double exact_norm_H1;
    double numeric_norm_H1;
};
// ***************************************************************************
} // namespace PCC_2D
} // namespace Assembler_Utilities
} // namespace PDETools
} // namespace Polydim

#endif
