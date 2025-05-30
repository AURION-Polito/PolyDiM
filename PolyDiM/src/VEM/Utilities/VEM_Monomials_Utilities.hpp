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

#ifndef __VEM_Monomials_Utilities_HPP
#define __VEM_Monomials_Utilities_HPP

#include "VEM_Monomials_Data.hpp"

namespace Polydim
{
namespace VEM
{
namespace Utilities
{
template <unsigned short dimension> struct VEM_Monomials_Utilities final
{
    Eigen::MatrixXi Exponents(const VEM_Monomials_Data &data) const;

    Eigen::MatrixXd Vander(const VEM_Monomials_Data &data,
                           const Eigen::MatrixXd &points,
                           const Eigen::Vector3d &centroid,
                           const double &diam) const;

    template <typename VEM_MonomialType>
    std::vector<Eigen::MatrixXd> VanderDerivatives(const VEM_Monomials_Data &data,
                                                   const VEM_MonomialType &monomials,
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

    template <typename VEM_MonomialType>
    Eigen::MatrixXd VanderLaplacian(const VEM_Monomials_Data &data,
                                    const VEM_MonomialType &monomials,
                                    const Eigen::MatrixXd &Vander,
                                    const double &diam) const
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
                           Eigen::MatrixXd &Qmatrix) const;
};
} // namespace Utilities
} // namespace VEM
} // namespace Polydim

#endif
