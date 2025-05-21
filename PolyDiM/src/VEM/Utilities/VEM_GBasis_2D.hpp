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

#ifndef __VEM_GBasis_VEM_HPP
#define __VEM_GBasis_VEM_HPP

#include "VEM_GBasis_Data.hpp"
#include "VEM_Monomials_2D.hpp"

namespace Polydim
{
namespace VEM
{
namespace Utilities
{
class VEM_GBasis_2D final
{
  private:
    VEM_Monomials_2D monomials;

    std::vector<Eigen::Vector2i> VectorDecompositionIndices(const VEM_GBasis_Data &data, const Eigen::VectorXi &expo) const;

  public:
    VEM_GBasis_Data Compute(const unsigned int polynomial_degree);

    std::vector<Eigen::MatrixXd> VanderGBigOPlus(const VEM_GBasis_Data &data, const Eigen::MatrixXd &vander) const
    {
        std::vector<Eigen::MatrixXd> vanderGBigOPlus(2, Eigen::MatrixXd::Zero(vander.rows(), data.Nkm1));
        vanderGBigOPlus[0] = vander.leftCols(data.Nkm1).array().colwise() * vander.col(2).array();
        vanderGBigOPlus[1] = vander.leftCols(data.Nkm1).array().colwise() * (-1.0 * vander.col(1).array());
        return vanderGBigOPlus;
    };
};
} // namespace Utilities
} // namespace VEM
} // namespace Polydim

#endif
