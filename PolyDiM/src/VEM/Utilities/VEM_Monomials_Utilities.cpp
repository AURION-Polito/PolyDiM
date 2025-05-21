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
template struct VEM_Monomials_Utilities<1>;
template struct VEM_Monomials_Utilities<2>;
template struct VEM_Monomials_Utilities<3>;
//****************************************************************************
template <unsigned short dimension>
MatrixXi VEM_Monomials_Utilities<dimension>::Exponents(const VEM_Monomials_Data &data) const
{
    MatrixXi exponents(dimension, data.NumMonomials);

    for (unsigned int m = 0; m < data.NumMonomials; m++)
        exponents.col(m) << data.Exponents[m];

    return exponents;
}
//****************************************************************************
template <unsigned short dimension>
MatrixXd VEM_Monomials_Utilities<dimension>::Vander(const VEM_Monomials_Data &data,
                                                    const MatrixXd &points,
                                                    const Vector3d &centroid,
                                                    const double &diam) const
{
    MatrixXd vander;
    const unsigned int numPoints = points.cols();
    if (data.NumMonomials > 1)
    {
        // VanderPartial[i]'s rows contain (x-x_E)^i/h_E^i,
        // (y-y_E)^i/h_E^i and (possibly) (z-z_E)^i/h_E^i respectively.
        // Size is dimension x numPoints.
        vector<MatrixXd> VanderPartial(data.PolynomialDegree + 1, MatrixXd(dimension, numPoints));
        double inverseDiam = 1.0 / diam;
        VanderPartial[0].setOnes(dimension, numPoints);
        VanderPartial[1] = (points.colwise() - centroid) * inverseDiam;

        for (unsigned int i = 2; i <= data.PolynomialDegree; i++)
            VanderPartial[i] = VanderPartial[i - 1].cwiseProduct(VanderPartial[1]);

        vander.resize(numPoints, data.NumMonomials);
        vander.col(0).setOnes();
        for (unsigned int i = 1; i < data.NumMonomials; ++i)
        {
            const VectorXi expo = data.Exponents[i];

            vander.col(i) = (VanderPartial[expo[0]].row(0)).transpose();
            if (dimension > 1)
                vander.col(i) = vander.col(i).cwiseProduct(VanderPartial[expo[1]].row(1).transpose());
            if (dimension > 2)
                vander.col(i) = vander.col(i).cwiseProduct(VanderPartial[expo[2]].row(2).transpose());
        }
    }
    else
        vander.setOnes(numPoints, 1);

    return vander;
}
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
