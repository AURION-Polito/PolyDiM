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

#ifndef __PDETOOLS_ASSEMBLER_assembler_PCC_2D_functions_utilities_HPP
#define __PDETOOLS_ASSEMBLER_assembler_PCC_2D_functions_utilities_HPP

#include "Eigen_Array.hpp"
#include "Eigen_SparseArray.hpp"
#include "assembler_PCC_2D_functions_data.hpp"

namespace Polydim
{
namespace PDETools
{
namespace Assembler_Utilities
{
namespace PCC_2D
{
// ***************************************************************************
std::list<Eigen::Triplet<double>> to_triplets(const Eigen::SparseMatrix<double> &M);
// ***************************************************************************
inline Eigen::VectorXd to_VectorXd(const Gedim::Eigen_Array<> &v)
{
    return static_cast<const Eigen::VectorXd &>(v);
}
// ***************************************************************************
inline Gedim::Eigen_Array<> to_Eigen_Array(const Eigen::VectorXd &v)
{
    return Gedim::Eigen_Array<>(v);
}
// ***************************************************************************
Gedim::Eigen_SparseArray<> to_Eigen_SparseArray(const Sparse_Matrix_Data &A);
// ***************************************************************************
Gedim::Eigen_SparseArray<> to_Eigen_SparseArray(const Sparse_Matrix_Data &A,
                                                const std::array<unsigned int, 2> &new_size,
                                                const std::array<unsigned int, 2> &shifts = {0, 0},
                                                const bool transpose = false);
// ***************************************************************************
Sparse_Matrix_Data to_Sparse_Matrix_Data(const Gedim::Eigen_SparseArray<> &A);
// ***************************************************************************
inline Eigen::SparseMatrix<double> to_SparseMatrix(const Gedim::Eigen_SparseArray<> &A)
{
    return static_cast<const Eigen::SparseMatrix<double> &>(A);
}
// ***************************************************************************
inline Gedim::Eigen_SparseArray<> to_Eigen_SparseArray(const Eigen::SparseMatrix<double> &A)
{
    return Gedim::Eigen_SparseArray<>(A);
}
// ***************************************************************************
} // namespace PCC_2D
} // namespace Assembler_Utilities
} // namespace PDETools
} // namespace Polydim

#endif
