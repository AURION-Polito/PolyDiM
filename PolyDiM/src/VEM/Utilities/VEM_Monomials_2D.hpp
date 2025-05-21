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

#ifndef __VEM_Monomials_2D_HPP
#define __VEM_Monomials_2D_HPP

#include "VEM_Monomials_Utilities.hpp"

namespace Polydim
{
namespace VEM
{
namespace Utilities
{
class VEM_Monomials_2D final
{
  private:
    VEM_Monomials_Utilities<2> utilities;

  public:
    VEM_Monomials_Data Compute(const unsigned int polynomial_degree) const;

    inline Eigen::MatrixXi Exponents(const VEM_Monomials_Data &data) const
    {
        return utilities.Exponents(data);
    }

    inline Eigen::MatrixXd DerivativeMatrix(const VEM_Monomials_Data &data, const unsigned int &i) const
    {
        return data.DerivativeMatrices[i];
    }

    inline Eigen::MatrixXd D_x(const VEM_Monomials_Data &data) const
    {
        return data.DerivativeMatrices[0];
    }

    inline Eigen::MatrixXd D_y(const VEM_Monomials_Data &data) const
    {
        return data.DerivativeMatrices[1];
    }

    int Index(const Eigen::VectorXi &exponents) const;

    std::vector<int> DerivativeIndices(const VEM_Monomials_Data &data, const unsigned int &index) const;

    std::vector<int> SecondDerivativeIndices(const VEM_Monomials_Data &data, const unsigned int &index) const;

    inline Eigen::MatrixXd Vander(const VEM_Monomials_Data &data,
                                  const Eigen::MatrixXd &points,
                                  const Eigen::Vector3d &centroid,
                                  const double &diam) const
    {
        return utilities.Vander(data, points, centroid, diam);
    }

    inline std::vector<Eigen::MatrixXd> VanderDerivatives(const VEM_Monomials_Data &data,
                                                          const Eigen::MatrixXd &vander,
                                                          const double &diam) const
    {
        return utilities.VanderDerivatives(data, (*this), vander, diam);
    }

    inline Eigen::MatrixXd VanderLaplacian(const VEM_Monomials_Data &data, const Eigen::MatrixXd &vander, const double &diam) const
    {
        return utilities.VanderLaplacian(data, (*this), vander, diam);
    }

    inline void MGSOrthonormalize(const Eigen::VectorXd &weights,
                                  const Eigen::MatrixXd &Vander,
                                  Eigen::MatrixXd &Hmatrix,
                                  Eigen::MatrixXd &QmatrixInv,
                                  Eigen::MatrixXd &Qmatrix) const
    {
        return utilities.MGSOrthonormalize(weights, Vander, Hmatrix, QmatrixInv, Qmatrix);
    };
};
} // namespace Utilities
} // namespace VEM
} // namespace Polydim

#endif
