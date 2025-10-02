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

#ifndef __Monomials_3D_HPP
#define __Monomials_3D_HPP

#include "Monomials_Utilities.hpp"

namespace Polydim
{
namespace Utilities
{
class Monomials_3D final
{
  private:
    Polydim::Utilities::Monomials_Utilities<3> utilities;

  public:
    Polydim::Utilities::Monomials_Data Compute(const unsigned int polynomial_degree) const;

    inline Eigen::MatrixXi Exponents(const Polydim::Utilities::Monomials_Data &data) const
    {
        return utilities.Exponents(data);
    }

    inline Eigen::MatrixXd DerivativeMatrix(const Polydim::Utilities::Monomials_Data &data, const unsigned int &i) const
    {
        return data.DerivativeMatrices[i];
    }

    inline Eigen::MatrixXd D_x(const Polydim::Utilities::Monomials_Data &data) const
    {
        return data.DerivativeMatrices[0];
    }

    inline Eigen::MatrixXd D_y(const Polydim::Utilities::Monomials_Data &data) const
    {
        return data.DerivativeMatrices[1];
    }

    inline Eigen::MatrixXd D_z(const Polydim::Utilities::Monomials_Data &data) const
    {
        return data.DerivativeMatrices[2];
    }

    int Index(const Eigen::VectorXi &exponents) const;
    std::vector<int> DerivativeIndices(const Polydim::Utilities::Monomials_Data &data, const unsigned int &index) const;
    std::vector<int> SecondDerivativeIndices(const Polydim::Utilities::Monomials_Data &data, const unsigned int &index) const;

    inline Eigen::MatrixXd Vander(const Polydim::Utilities::Monomials_Data &data,
                                  const Eigen::MatrixXd &points,
                                  const Eigen::Vector3d &centroid,
                                  const double &diam) const
    {
        return utilities.Vander(data, points, centroid, diam);
    }
    inline std::vector<Eigen::MatrixXd> VanderDerivatives(const Polydim::Utilities::Monomials_Data &data,
                                                          const Eigen::MatrixXd &vander,
                                                          const double &diam) const
    {
        return utilities.VanderDerivatives(data, (*this), vander, diam);
    }

    inline Eigen::MatrixXd VanderLaplacian(const Polydim::Utilities::Monomials_Data &data,
                                           const Eigen::MatrixXd &vander,
                                           const double &diam) const
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
} // namespace Polydim

#endif
