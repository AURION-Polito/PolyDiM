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

#include "VEM_Monomials_Utilities.hpp"

#include "LAPACK_utilities.hpp"

using namespace std;
using namespace Eigen;

namespace Polydim
{
namespace VEM
{
namespace Utilities
{

//****************************************************************************
template <unsigned short dimension>
void VEM_Monomials_Utilities<dimension>::MGSOrthonormalize(const Eigen::VectorXd &weights,
                                                           const Eigen::MatrixXd &Vander,
                                                           Eigen::MatrixXd &Hmatrix,
                                                           Eigen::MatrixXd &QmatrixInv,
                                                           Eigen::MatrixXd &Qmatrix) const
{
    MatrixXd Q1;
    MatrixXd R1;
    LAPACK_utilities::MGS(Vander, Q1, R1);

    // L2(E)-re-orthogonalization process
    MatrixXd Q2;
    MatrixXd R2;
    LAPACK_utilities::MGS(weights.array().sqrt().matrix().asDiagonal() * Q1, Q2, R2);

    Hmatrix = Q2.transpose() * Q2;

    QmatrixInv = (R2 * R1).transpose();
    LAPACK_utilities::inverseTri(QmatrixInv, Qmatrix, 'L', 'N');
}
//****************************************************************************
} // namespace Utilities
} // namespace VEM
} // namespace Polydim
