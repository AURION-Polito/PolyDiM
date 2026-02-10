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

#ifndef __ZFEM_PCC_PerformanceAnalysis_HPP
#define __ZFEM_PCC_PerformanceAnalysis_HPP

#include "Eigen/Eigen"
#include <vector>

namespace Polydim
{
namespace ZFEM
{
namespace PCC
{
struct ZFEM_PCC_PerformanceAnalysis_Data
{
    Eigen::VectorXd ConsistencyError;
    std::vector<Eigen::VectorXd> DerivativesConsistencyError;
};

struct ZFEM_PCC_PerformanceAnalysis final
{
    template <typename ZFEM_Monomials_Type, typename ZFEM_Monomials_Data_Type, typename ZFEM_LocalSpace_Type, typename ZFEM_LocalSpaceData_Type>
    ZFEM_PCC_PerformanceAnalysis_Data Compute2D(const ZFEM_Monomials_Type &monomials,
                                                const ZFEM_Monomials_Data_Type &monomials_data,
                                                const ZFEM_LocalSpace_Type &,
                                                const ZFEM_LocalSpaceData_Type &ZFEM_local_space_data) const
    {
        const Eigen::MatrixXd VanderInternal = monomials.Vander(monomials_data,
                                                                ZFEM_local_space_data.InternalQuadrature.Points,
                                                                ZFEM_local_space_data.KernelIncenter,
                                                                ZFEM_local_space_data.Diameter);

        const auto VanderInternalDerivatives =
            monomials.VanderDerivatives(monomials_data, VanderInternal, ZFEM_local_space_data.Diameter);

        ZFEM_PCC_PerformanceAnalysis_Data result;
        result.ConsistencyError = (ZFEM_local_space_data.ZFEM_basis_functions_values * ZFEM_local_space_data.Dmatrix - VanderInternal)
                                      .array()
                                      .square()
                                      .matrix()
                                      .transpose() *
                                  ZFEM_local_space_data.InternalQuadrature.Weights;

        result.DerivativesConsistencyError.resize(ZFEM_local_space_data.Dimension);
        for (unsigned int d = 0; d < ZFEM_local_space_data.Dimension; d++)
            result.DerivativesConsistencyError[d] =
                (ZFEM_local_space_data.ZFEM_basis_functions_derivative_values[d] * ZFEM_local_space_data.Dmatrix -
                 VanderInternalDerivatives[d])
                    .array()
                    .square()
                    .matrix()
                    .transpose() *
                ZFEM_local_space_data.InternalQuadrature.Weights;

        return result;
    }
};

} // namespace PCC
} // namespace ZFEM
} // namespace Polydim

#endif
