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
struct Exact_Solution_Data final
{
    Eigen::VectorXd exact_solution;
    Eigen::VectorXd exact_solution_strong;
};
// ***************************************************************************
struct Variational_Operator final
{
    Sparse_Matrix_Data A;
    Sparse_Matrix_Data A_Strong;
};
// ***************************************************************************
struct Post_Process_Data final
{
    Eigen::VectorXd cell0Ds_numeric;
    Eigen::VectorXd cell0Ds_exact;
    Eigen::VectorXd cell2Ds_exact_norm_L2;
    Eigen::VectorXd cell2Ds_numeric_norm_L2;
    Eigen::VectorXd cell2Ds_error_L2;
    double mesh_size;
    double error_L2;
    double exact_norm_L2;
    double numeric_norm_L2;
};
// ***************************************************************************
} // namespace PCC_2D
} // namespace Assembler_Utilities
} // namespace PDETools
} // namespace Polydim

#endif
