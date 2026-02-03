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

#include "LAPACK_utilities.hpp"
#include "VEM_PCC_Utilities.hpp"

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

    // 3D
    double PiNablaConditioning = -1.0;   ///< conditioning of piNabla
    double Pi0km1Conditioning = -1.0;    ///< conditioning of pi0km1
    double Pi0kConditioning = -1.0;      ///< conditioning of pi0k
    double ErrorPiNabla = -1.0;          ///< |piNabla * Dofs - I|
    double ErrorPi0km1 = -1.0;           ///< |pi0km1 * Dofs.leftCols(Nkm1) - I.topLeftCorner(Nkm1, Nkm1)|
    double ErrorPi0k = -1.0;             ///< |pi0k * Dofs - I|
    std::vector<double> ErrorPi0km1Grad; ///< Error of Pi0km1Grad, size geometric dimension
    double StabNorm = -1.0;              ///< Norm of S
    double ErrorStabilization = -1.0;    ///< |S * Dofs|
    double ErrorHCD = -1.0;              ///< |H - CD|
    double ErrorGBD = -1.0;              ///< |G - BD|
    std::vector<double> ErrorHED;        ///< |H - ED|
};

struct ZFEM_PCC_PerformanceAnalysis final
{
    template <typename ZFEM_Monomials_Type, typename ZFEM_Monomials_Data_Type, typename ZFEM_LocalSpace_Type, typename ZFEM_LocalSpaceData_Type>
    ZFEM_PCC_PerformanceAnalysis_Data Compute2D(const ZFEM_Monomials_Type &,
                                                const ZFEM_Monomials_Data_Type &,
                                                const ZFEM_LocalSpace_Type &,
                                                const ZFEM_LocalSpaceData_Type &ZFEM_local_space_data) const
    {
        ZFEM_PCC_PerformanceAnalysis_Data result;
        result.ConsistencyError =
            (ZFEM_local_space_data.ZFEM_basis_functions_values * ZFEM_local_space_data.Dmatrix - ZFEM_local_space_data.VanderInternal)
                .array()
                .square()
                .matrix()
                .transpose() *
            ZFEM_local_space_data.InternalQuadrature.Weights;

        result.DerivativesConsistencyError.resize(ZFEM_local_space_data.Dimension);
        for (unsigned int d = 0; d < ZFEM_local_space_data.Dimension; d++)
            result.DerivativesConsistencyError[d] =
                (ZFEM_local_space_data.ZFEM_basis_functions_derivative_values[d] * ZFEM_local_space_data.Dmatrix -
                 ZFEM_local_space_data.VanderInternalDerivatives[d])
                    .array()
                    .square()
                    .matrix()
                    .transpose() *
                ZFEM_local_space_data.InternalQuadrature.Weights;

        return result;
    }

    template <typename ZFEM_Monomials_Type, typename ZFEM_Monomials_Data_Type, typename ZFEM_LocalSpace_Type, typename ZFEM_LocalSpaceData_Type>
    ZFEM_PCC_PerformanceAnalysis_Data Compute3D(const ZFEM_Monomials_Type &ZFEM_monomials,
                                                const ZFEM_Monomials_Data_Type &ZFEM_monomials_data,
                                                const ZFEM_LocalSpace_Type &ZFEM_local_space,
                                                const ZFEM_LocalSpaceData_Type &ZFEM_local_space_data) const
    {
        ZFEM_PCC_PerformanceAnalysis_Data result;

        const double invDiameter = 1.0 / ZFEM_local_space_data.Diameter;

        const Eigen::MatrixXd &piNabla = ZFEM_local_space_data.PiNabla;
        const Eigen::MatrixXd &pi0km1 = ZFEM_local_space_data.Pi0km1;
        const Eigen::MatrixXd &pi0k = ZFEM_local_space_data.Pi0k;

        result.PiNablaConditioning = LAPACK_utilities::cond(LAPACK_utilities::svd(piNabla));
        result.Pi0kConditioning = LAPACK_utilities::cond(LAPACK_utilities::svd(pi0k));
        result.Pi0km1Conditioning = LAPACK_utilities::cond(LAPACK_utilities::svd(pi0km1));

        const Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(ZFEM_local_space_data.NumProjectorBasisFunctions,
                                                                   ZFEM_local_space_data.NumProjectorBasisFunctions);

        const unsigned int Nkm1 = pi0km1.rows();
        const unsigned int Nk = piNabla.rows();

        result.ErrorPiNabla = (piNabla * ZFEM_local_space_data.Dmatrix - identity).norm() / identity.norm();
        result.ErrorPi0km1 =
            (pi0km1 * ZFEM_local_space_data.Dmatrix.leftCols(Nkm1) - identity.topLeftCorner(Nkm1, Nkm1)).norm() /
            identity.topLeftCorner(Nkm1, Nkm1).norm();

        result.ErrorPi0k = (pi0k * ZFEM_local_space_data.Dmatrix - identity).norm() / identity.norm();

        result.ErrorPi0km1Grad.resize(ZFEM_local_space_data.Pi0km1Der.size());
        for (unsigned int d = 0; d < ZFEM_local_space_data.Pi0km1Der.size(); ++d)
        {
            const Eigen::MatrixXd &piDerkm1 = ZFEM_local_space_data.Pi0km1Der[d];
            const Eigen::MatrixXd derMatrix =
                invDiameter * (ZFEM_monomials.DerivativeMatrix(ZFEM_monomials_data, d).topLeftCorner(Nk, Nkm1)).transpose();
            double relErrDenominator = (Nkm1 > 1) ? derMatrix.norm() : 1.0;
            result.ErrorPi0km1Grad[d] = (piDerkm1 * ZFEM_local_space_data.Dmatrix - derMatrix).norm() / relErrDenominator;
        }

        const Eigen::MatrixXd StabMatrix =
            ZFEM_local_space.ComputeDofiDofiStabilizationMatrix(ZFEM_local_space_data, VEM::PCC::ProjectionTypes::PiNabla);

        const Eigen::MatrixXd &stabilizationMatrix = StabMatrix;
        result.StabNorm = StabMatrix.norm();
        result.ErrorStabilization = (stabilizationMatrix * ZFEM_local_space_data.Dmatrix).norm();

        if (ZFEM_local_space_data.Hmatrix.size() > 0 && ZFEM_local_space_data.Cmatrix.size() > 0)
            result.ErrorHCD =
                (ZFEM_local_space_data.Hmatrix - ZFEM_local_space_data.Cmatrix * ZFEM_local_space_data.Dmatrix).norm() /
                ZFEM_local_space_data.Hmatrix.norm();
        if (ZFEM_local_space_data.Gmatrix.size() > 0 && ZFEM_local_space_data.Bmatrix.size() > 0)
            result.ErrorGBD =
                (ZFEM_local_space_data.Gmatrix - ZFEM_local_space_data.Bmatrix * ZFEM_local_space_data.Dmatrix).norm() /
                ZFEM_local_space_data.Gmatrix.norm();

        result.ErrorHED.resize(ZFEM_local_space_data.Pi0km1Der.size(), -1.0);
        for (unsigned int d = 0; d < ZFEM_local_space_data.Pi0km1Der.size(); d++)
        {
            const Eigen::MatrixXd derMatrix =
                invDiameter *
                ZFEM_monomials.DerivativeMatrix(ZFEM_monomials_data, d).topLeftCorner(Nk, ZFEM_local_space_data.Nkm1) *
                ZFEM_local_space_data.Hmatrix.topLeftCorner(ZFEM_local_space_data.Nkm1, ZFEM_local_space_data.Nkm1);

            result.ErrorHED[d] =
                (derMatrix.transpose() - ZFEM_local_space_data.Ematrix[d] * ZFEM_local_space_data.Dmatrix).norm() /
                derMatrix.norm();
        }

        return result;
    }
};

} // namespace PCC
} // namespace ZFEM
} // namespace Polydim

#endif
