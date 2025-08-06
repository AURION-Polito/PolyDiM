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

#include "VEM_DF_PCC_Utilities.hpp"


using namespace Eigen;
using namespace std;

namespace Polydim
{
namespace VEM
{
namespace DF_PCC
{
template struct VEM_DF_PCC_Utilities<2>;
template struct VEM_DF_PCC_Utilities<3>;
//****************************************************************************
template <unsigned short dimension>
Eigen::VectorXd VEM_DF_PCC_Utilities<dimension>::ComputeEdgeBasisCoefficients(const unsigned int &order,
                                                                              const Eigen::VectorXd &edgeInternalPoints) const

//****************************************************************************
template <unsigned short dimension>
MatrixXd VEM_DF_PCC_Utilities<dimension>::ComputeValuesOnEdge(const Eigen::RowVectorXd &edgeInternalPoints,
                                                              const unsigned int &order,
                                                              const Eigen::VectorXd &edgeBasisCoefficients,
                                                              const Eigen::VectorXd &pointsCurvilinearCoordinates) const

//****************************************************************************
template <unsigned short dimension>
MatrixXd VEM_DF_PCC_Utilities<dimension>::ComputeDofiDofiStabilizationMatrix(const std::vector<MatrixXd> &projector,
                                                                             const double &coefficient,
                                                                             const std::vector<Eigen::MatrixXd> &dmatrix) const

//****************************************************************************
} // namespace DF_PCC
} // namespace VEM
} // namespace Polydim
