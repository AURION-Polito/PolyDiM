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

#ifndef __Monomials_Utilities_HPP
#define __Monomials_Utilities_HPP

#include "LAPACK_utilities.hpp"
#include "Monomials_Data.hpp"

namespace Polydim
{
namespace Utilities
{
template <unsigned short dimension> struct Monomials_Utilities final
{
    Eigen::MatrixXi Exponents(const Monomials_Data &data) const
    {
        Eigen::MatrixXi exponents(dimension, data.NumMonomials);

        for (unsigned int m = 0; m < data.NumMonomials; m++)
            exponents.col(m) << data.Exponents[m];

        return exponents;
    }

    Eigen::MatrixXd Vander(const Monomials_Data &data, const Eigen::MatrixXd &points, const Eigen::Vector3d &centroid, const double &diam) const
    {
        Eigen::MatrixXd vander;
        const unsigned int numPoints = points.cols();
        if (data.NumMonomials > 1)
        {
            // VanderPartial[i]'s rows contain (x-x_E)^i/h_E^i,
            // (y-y_E)^i/h_E^i and (possibly) (z-z_E)^i/h_E^i respectively.
            // Size is dimension x numPoints.
            std::vector<Eigen::MatrixXd> VanderPartial(data.PolynomialDegree + 1, Eigen::MatrixXd(dimension, numPoints));
            double inverseDiam = 1.0 / diam;
            VanderPartial[0].setOnes(dimension, numPoints);
            VanderPartial[1] = (points.colwise() - centroid) * inverseDiam;

            for (unsigned int i = 2; i <= data.PolynomialDegree; i++)
                VanderPartial[i] = VanderPartial[i - 1].cwiseProduct(VanderPartial[1]);

            vander.resize(numPoints, data.NumMonomials);
            vander.col(0).setOnes();
            for (unsigned int i = 1; i < data.NumMonomials; ++i)
            {
                const Eigen::VectorXi expo = data.Exponents[i];

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

    template <typename MonomialType>
    std::vector<Eigen::MatrixXd> VanderDerivatives(const Monomials_Data &data,
                                                   const MonomialType &monomials,
                                                   const Eigen::MatrixXd &Vander,
                                                   const double &diam) const
    {
        std::vector<Eigen::MatrixXd> vanderDerivatives;
        vanderDerivatives.resize(dimension);
        for (unsigned int i = 0; i < dimension; i++)
        {
            vanderDerivatives[i].resizeLike(Vander);
            vanderDerivatives[i].col(0).setZero();
        }
        if (data.NumMonomials > 1)
        {
            double inverseDiam = 1.0 / diam;
            for (unsigned int k = 1; k < data.NumMonomials; k++)
            {
                std::vector<int> derIndices = monomials.DerivativeIndices(data, k);
                for (unsigned int i = 0; i < dimension; i++)
                {
                    if (derIndices[i] >= 0)
                        vanderDerivatives[i].col(k) =
                            inverseDiam * monomials.DerivativeMatrix(data, i)(k, derIndices[i]) * Vander.col(derIndices[i]);
                    else
                        vanderDerivatives[i].col(k).setZero();
                }
            }
        }

        return vanderDerivatives;
    }

    template <typename MonomialType>
    Eigen::MatrixXd VanderLaplacian(const Monomials_Data &data, const MonomialType &monomials, const Eigen::MatrixXd &Vander, const double &diam) const
    {
        Eigen::MatrixXd vanderLaplacian;

        vanderLaplacian.resizeLike(Vander);
        vanderLaplacian.block(0, 0, Vander.rows(), 3).setZero();
        Eigen::MatrixXd laplacian = data.Laplacian;

        if (data.NumMonomials > 3)
        {
            const double inverseDiamSqrd = 1.0 / (diam * diam);
            for (unsigned int k = 3; k < data.NumMonomials; k++)
            {
                std::vector<int> secondDerIndices = monomials.SecondDerivativeIndices(data, k);
                if (secondDerIndices[0] >= 0)
                    vanderLaplacian.col(k) =
                        inverseDiamSqrd * laplacian(k, secondDerIndices[0]) * Vander.col(secondDerIndices[0]);
                else
                    vanderLaplacian.col(k).setZero();
                for (unsigned int i = 1; i < dimension; i++)
                {
                    if (secondDerIndices[i] >= 0)
                        vanderLaplacian.col(k) +=
                            inverseDiamSqrd * laplacian(k, secondDerIndices[i]) * Vander.col(secondDerIndices[i]);
                }
            }
        }

        return vanderLaplacian;
    }

    void MGSOrthonormalize(const Eigen::VectorXd &weights,
                           const Eigen::MatrixXd &Vander,
                           Eigen::MatrixXd &Hmatrix,
                           Eigen::MatrixXd &QmatrixInv,
                           Eigen::MatrixXd &Qmatrix) const
    {
        Eigen::MatrixXd Q1;
        Eigen::MatrixXd R1;
        LAPACK_utilities::MGS(Vander, Q1, R1);

        // L2(E)-re-orthogonalization process
        Eigen::MatrixXd Q2;
        Eigen::MatrixXd R2;
        LAPACK_utilities::MGS(weights.array().sqrt().matrix().asDiagonal() * Q1, Q2, R2);

        Hmatrix = Q2.transpose() * Q2;

        QmatrixInv = (R2 * R1).transpose();
        LAPACK_utilities::inverseTri(QmatrixInv, Qmatrix, 'L', 'N');
    }
};
} // namespace Utilities
} // namespace Polydim

#endif
