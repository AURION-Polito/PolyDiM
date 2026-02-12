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

#include "assembler_PCC_2D_functions_utilities.hpp"

namespace Polydim
{
namespace PDETools
{
namespace Assembler_Utilities
{
namespace PCC_2D
{
// ***************************************************************************
  std::list<Eigen::Triplet<double> > to_triplets(const Eigen::SparseMatrix<double>& M)
  {
    std::list<Eigen::Triplet<double>> v;

    for (int i = 0; i < M.outerSize(); i++)
    {
      for (typename Eigen::SparseMatrix<double>::InnerIterator it(M, i); it; ++it)
      {
        if (it.value() != 0.0)
          v.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
      }
    }

    return v;
  }
// ***************************************************************************
  Gedim::Eigen_SparseArray<> to_Eigen_SparseArray(const Sparse_Matrix_Data& A)
  {
    Gedim::Eigen_SparseArray<> eigen_A;
    eigen_A.SetSize(A.size.at(0), A.size.at(1), Gedim::ISparseArray::SparseArrayTypes::None);

    eigen_A.Triplets(A.rows, A.cols, A.values);

    eigen_A.Create();

    return eigen_A;
  }
// ***************************************************************************
  Sparse_Matrix_Data to_Sparse_Matrix_Data(const Gedim::Eigen_SparseArray<>& A)
  {
    Sparse_Matrix_Data result;

    const auto &eigen_A = static_cast<const Eigen::SparseMatrix<double> &>(A);
    const auto A_triplets = to_triplets(eigen_A);
    const auto num_triplets = A_triplets.size();

    result.size = {static_cast<unsigned int>(eigen_A.rows()), static_cast<unsigned int>(eigen_A.cols())};
    result.rows.resize(num_triplets);
    result.cols.resize(num_triplets);
    result.values.resize(num_triplets);

    unsigned int t = 0;
    for (const auto &triplet : A_triplets)
    {
      result.rows[t] = triplet.row();
      result.cols[t] = triplet.col();
      result.values[t] = triplet.value();
      t++;
    }

    return result;
  }
  // ***************************************************************************
} // namespace PCC_2D
} // namespace Assembler_Utilities
} // namespace PDETools
} // namespace Polydim
