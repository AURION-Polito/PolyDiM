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

#ifndef __VEM_MCC_PerformanceAnalysis_HPP
#define __VEM_MCC_PerformanceAnalysis_HPP

#include "Eigen/Eigen"

#include "LAPACK_utilities.hpp"
#include "VEM_MCC_Utilities.hpp"

namespace Polydim
{
namespace VEM
{
namespace MCC
{
struct VEM_MCC_PerformanceAnalysis_Data
{
    double VmatrixConditioning = -1.0;
    double HmatrixConditioning = -1.0;
    double Pi0kConditioning = -1.0;
    double GmatrixConditioning = -1.0;
    double ErrorPi0k = -1.0;
    double StabNorm = -1.0;
    double ErrorStabilization = -1.0;
    double ErrorGBD = -1.0;
};

struct VEM_MCC_PerformanceAnalysis final
{
    template <typename VEM_Monomials_Type, typename VEM_Monomials_Data_Type, typename VEM_LocalSpace_Type, typename VEM_LocalSpaceData_Type>
    VEM_MCC_PerformanceAnalysis_Data Compute(const VEM_Monomials_Type &vem_monomials,
                                             const VEM_Monomials_Data_Type &vem_monomials_data,
                                             const VEM_LocalSpace_Type &vem_local_space,
                                             const VEM_LocalSpaceData_Type &vem_local_space_data) const
    {
        VEM_MCC_PerformanceAnalysis_Data result;

        const Eigen::MatrixXd &Vmatrix = vem_local_space_data.Vmatrix;
        const Eigen::MatrixXd &Hmatrix = vem_local_space_data.Hmatrix;
        const Eigen::MatrixXd &Gmatrix = vem_local_space_data.Gmatrix;
        const Eigen::MatrixXd &pi0k = vem_local_space_data.Pi0k;

        result.VmatrixConditioning = LAPACK_utilities::cond(LAPACK_utilities::svd(Vmatrix));
        result.HmatrixConditioning = LAPACK_utilities::cond(LAPACK_utilities::svd(Hmatrix));
        result.Pi0kConditioning = LAPACK_utilities::cond(LAPACK_utilities::svd(pi0k));
        result.GmatrixConditioning = LAPACK_utilities::cond(LAPACK_utilities::svd(Gmatrix));

        const unsigned int Nk = vem_local_space_data.Nk;
        const unsigned int dimension = vem_local_space_data.Dimension;

        const Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(dimension * Nk, dimension * Nk);

        const Eigen::MatrixXd &polynomialBasisDofs = vem_local_space_data.Dmatrix;
        result.ErrorPi0k = (pi0k * polynomialBasisDofs - identity).norm() / identity.norm();

        const Eigen::MatrixXd stabilizationMatrix =
            vem_local_space.ComputeDofiDofiStabilizationMatrix(vem_local_space_data, ProjectionTypes::Pi0k);
        result.StabNorm = stabilizationMatrix.norm();
        result.ErrorStabilization = (stabilizationMatrix * polynomialBasisDofs).norm();

        if (vem_local_space_data.Gmatrix.size() > 0 && vem_local_space_data.Bmatrix.size() > 0)
            result.ErrorGBD = (vem_local_space_data.Gmatrix - vem_local_space_data.Bmatrix * polynomialBasisDofs).norm() /
                              vem_local_space_data.Gmatrix.norm();

        return result;
    }
};

} // namespace MCC
} // namespace VEM
} // namespace Polydim

#endif
