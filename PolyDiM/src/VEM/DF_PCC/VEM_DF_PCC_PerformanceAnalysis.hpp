#ifndef __VEM_DF_PCC_PerformanceAnalysis_H
#define __VEM_DF_PCC_PerformanceAnalysis_H

#include "Eigen/Eigen"
#include <vector>

#include "LAPACK_utilities.hpp"
#include "VEM_DF_PCC_Utilities.hpp"

namespace Polydim
{
namespace VEM
{
namespace DF_PCC
{
struct VEM_DF_PCC_PerformanceAnalysis_Data
{
    std::vector<double> PiNablaConditioning; ///< conditioning of piNabla
    std::vector<double> Pi0km2Conditioning;  ///< conditioning of piNabla
    std::vector<double> Pi0kConditioning;    ///< conditioning of piNabla
    std::vector<double> ErrorPiNabla;        ///< |piNabla * Dofs - I|
    std::vector<double> ErrorPi0km2;         ///< |pi0km2 * Dofs.leftCols(Nkm2) - I.topLeftCorner(Nkm2, Nkm2)|
    std::vector<double> ErrorPi0k;           ///< |pi0k * Dofs - I|
    std::vector<double> ErrorPi0km1Grad;     ///< Error of Pi0km1Grad, size geometric dimension
    double StabNorm = -1.0;
    double ErrorStabilization = -1.0; ///< |S * Dofs|
    std::vector<double> ErrorHCD;     ///< |H - CD|
    std::vector<double> ErrorGBD;     ///< |G - BD|
    std::vector<double> ErrorHED;     ///< |H - ED|
};

struct VEM_DF_PCC_PerformanceAnalysis final
{
    template <typename VEM_Monomials_Type, typename VEM_Monomials_Data_Type, typename VEM_vem_local_space_data_Type, typename VEM_vem_local_space_dataData_Type>
    VEM_DF_PCC_PerformanceAnalysis_Data Compute(const VEM_Monomials_Type &vem_monomials,
                                                const VEM_Monomials_Data_Type &vem_monomials_data,
                                                const VEM_vem_local_space_data_Type &vem_local_space,
                                                const VEM_vem_local_space_dataData_Type &vem_local_space_data) const
    {
        VEM_DF_PCC_PerformanceAnalysis_Data result;

        const double invDiameter = 1.0 / vem_local_space_data.Diameter;

        const std::vector<Eigen::MatrixXd> &piNabla = vem_local_space_data.PiNabla;
        const std::vector<Eigen::MatrixXd> &pi0km2 = vem_local_space_data.Pi0km2;
        const std::vector<Eigen::MatrixXd> &pi0k = vem_local_space_data.Pi0k;
        const std::vector<Eigen::MatrixXd> &Dmatrix = vem_local_space_data.Dmatrix;

        result.PiNablaConditioning.resize(vem_local_space_data.Dimension, -1.0);
        result.Pi0kConditioning.resize(vem_local_space_data.Dimension, -1.0);
        result.Pi0km2Conditioning.resize(vem_local_space_data.Dimension, -1.0);
        for (unsigned int d = 0; d < vem_local_space_data.Dimension; d++)
        {
            result.PiNablaConditioning[d] = LAPACK_utilities::cond(LAPACK_utilities::svd(piNabla[0]));
            result.Pi0kConditioning[d] = LAPACK_utilities::cond(LAPACK_utilities::svd(pi0k[0]));
            result.Pi0km2Conditioning[d] = LAPACK_utilities::cond(LAPACK_utilities::svd(pi0km2[0]));
        }

        const Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(vem_local_space_data.Nk, vem_local_space_data.Nk);

        result.ErrorPiNabla.resize(vem_local_space_data.Dimension, 1.0);
        result.ErrorPi0km2.resize(vem_local_space_data.Dimension, 1.0);
        result.ErrorPi0k.resize(vem_local_space_data.Dimension, 1.0);
        result.ErrorPi0km1Grad.resize(vem_local_space_data.Dimension * vem_local_space_data.Dimension, 1.0);
        for (unsigned int d1 = 0; d1 < vem_local_space_data.Dimension; d1++)
        {
            result.ErrorPiNabla[d1] = (piNabla[d1] * Dmatrix[d1] - identity).norm() / identity.norm();
            result.ErrorPi0km2[d1] = (pi0km2[d1] * Dmatrix[d1].leftCols(vem_local_space_data.Nkm2) -
                                      identity.topLeftCorner(vem_local_space_data.Nkm2, vem_local_space_data.Nkm2))
                                         .norm() /
                                     identity.topLeftCorner(vem_local_space_data.Nkm2, vem_local_space_data.Nkm2).norm();

            result.ErrorPi0k[d1] = (pi0k[d1] * Dmatrix[d1] - identity).norm() / identity.norm();

            for (unsigned int d2 = 0; d2 < vem_local_space_data.Dimension; d2++)
            {
                const Eigen::MatrixXd &piDerkm1 = vem_local_space_data.Pi0km1Der[vem_local_space_data.Dimension * d1 + d2];
                const Eigen::MatrixXd derMatrix =
                    invDiameter * vem_monomials.DerivativeMatrix(vem_monomials_data, d2)
                                      .topLeftCorner(vem_local_space_data.Nk, vem_local_space_data.Nkm1);

                result.ErrorPi0km1Grad[vem_local_space_data.Dimension * d1 + d2] =
                    (vem_local_space_data.Pi0km1Der[vem_local_space_data.Dimension * d1 + d2] * vem_local_space_data.Dmatrix[d1] -
                     derMatrix.transpose())
                        .norm() /
                    derMatrix.norm();
            }
        }

        const Eigen::MatrixXd StabMatrix =
            vem_local_space.ComputeDofiDofiStabilizationMatrix(vem_local_space_data, ProjectionTypes::PiNabla);

        const Eigen::MatrixXd &stabilizationMatrix = StabMatrix;
        result.StabNorm = StabMatrix.norm();

        Eigen::MatrixXd polynomialBasisDofs =
            Eigen::MatrixXd::Zero(vem_local_space_data.NumBasisFunctions,
                                  vem_local_space_data.Dimension * vem_local_space_data.Nk);
        for (unsigned int d = 0; d < vem_local_space_data.Dimension; d++)
            polynomialBasisDofs.middleCols(d * vem_local_space_data.Nk, vem_local_space_data.Nk) =
                vem_local_space_data.Dmatrix[d];

        result.ErrorStabilization = (stabilizationMatrix * polynomialBasisDofs).norm();

        result.ErrorHCD.resize(vem_local_space_data.Dimension, 1.0);
        result.ErrorGBD.resize(vem_local_space_data.Dimension, 1.0);
        result.ErrorHED.resize(vem_local_space_data.Dimension * vem_local_space_data.Dimension, 1.0);

        if (vem_local_space_data.Hmatrix.size() > 0 && vem_local_space_data.Cmatrix.size() > 0)
        {
            const double HMatrixInvNorm = 1.0 / vem_local_space_data.Hmatrix.norm();
            for (unsigned int d1 = 0; d1 < vem_local_space_data.Dimension; d1++)
                result.ErrorHCD[d1] =
                    (vem_local_space_data.Hmatrix - vem_local_space_data.Cmatrix[d1] * vem_local_space_data.Dmatrix[d1]).norm() * HMatrixInvNorm;
        }

        if (vem_local_space_data.Gmatrix.size() > 0 && vem_local_space_data.Bmatrix.size() > 0)
        {
            const double GMatrixInvNorm = 1.0 / vem_local_space_data.Gmatrix.norm();
            for (unsigned int d1 = 0; d1 < vem_local_space_data.Dimension; d1++)
                result.ErrorGBD[d1] =
                    (vem_local_space_data.Gmatrix - vem_local_space_data.Bmatrix[d1] * vem_local_space_data.Dmatrix[d1]).norm() * GMatrixInvNorm;
        }

        if (vem_local_space_data.Hmatrix.size() > 0 && vem_local_space_data.Ematrix.size() > 0)
        {
            for (unsigned int d1 = 0; d1 < vem_local_space_data.Dimension; d1++)
            {
                for (unsigned int d2 = 0; d2 < vem_local_space_data.Dimension; d2++)
                {
                    const Eigen::MatrixXd &piDerkm1 = vem_local_space_data.Pi0km1Der[vem_local_space_data.Dimension * d1 + d2];
                    const Eigen::MatrixXd derMatrix =
                        invDiameter *
                        vem_monomials.DerivativeMatrix(vem_monomials_data, d2)
                            .topLeftCorner(vem_local_space_data.Nk, vem_local_space_data.Nkm1) *
                        vem_local_space_data.Hmatrix.topLeftCorner(vem_local_space_data.Nkm1, vem_local_space_data.Nkm1);

                    result.ErrorHED[vem_local_space_data.Dimension * d1 + d2] =
                        (vem_local_space_data.Ematrix[vem_local_space_data.Dimension * d1 + d2] * vem_local_space_data.Dmatrix[d1] -
                         derMatrix.transpose())
                            .norm() /
                        derMatrix.norm();
                }
            }
        }

        return result;
    }
};

} // namespace DF_PCC
} // namespace VEM
} // namespace Polydim

#endif
