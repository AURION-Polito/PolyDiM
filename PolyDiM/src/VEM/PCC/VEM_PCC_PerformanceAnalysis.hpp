#ifndef __VEM_PCC_PerformanceAnalysis_H
#define __VEM_PCC_PerformanceAnalysis_H

#include "Eigen/Eigen"
#include <vector>

#include "LAPACK_utilities.hpp"

namespace Polydim
{
namespace VEM
{
namespace PCC
{
struct VEM_PCC_PerformanceAnalysis_Data
{
    double PiNablaConditioning = -1.0; ///< conditioning of piNabla
    double Pi0km1Conditioning = -1.0; ///< conditioning of piNabla
    double Pi0kConditioning = -1.0; ///< conditioning of piNabla
    double ErrorPiNabla = -1.0; ///< |piNabla * Dofs - I|
    double ErrorPi0km1 = -1.0; ///< |pi0km1 * Dofs.leftCols(Nkm1) - I.topLeftCorner(Nkm1, Nkm1)|
    double ErrorPi0k = -1.0; ///< |pi0k * Dofs - I|
    std::vector<double> ErrorPi0km1Grad; ///< Error of Pi0km1Grad, size geometric dimension
    double StabNorm = -1.0;
    double ErrorStabilization = -1.0; ///< |S * Dofs|
    double ErrorHCD = -1.0; ///< |H - CD|
    double ErrorGBD = -1.0; ///< |G - BD|
    std::vector<double> ErrorHED; ///< |H - ED|
};

struct VEM_PCC_PerformanceAnalysis final
{
    template<typename VEM_Monomials_Type,
             typename VEM_Monomials_Data_Type,
             typename VEM_LocalSpace_Type,
             typename VEM_LocalSpaceData_Type>
    VEM_PCC_PerformanceAnalysis_Data Compute(const double& polytopeMeasure,
                                             const double& polytopeDiameter,
                                             const VEM_Monomials_Type& vem_monomials,
                                             const VEM_Monomials_Data_Type& vem_monomials_data,
                                             const VEM_LocalSpace_Type& vem_local_space,
                                             const VEM_LocalSpaceData_Type& vem_local_space_data) const
    {
        VEM_PCC_PerformanceAnalysis_Data result;

        double invDiameter = 1.0 / polytopeDiameter;

        const Eigen::MatrixXd& piNabla = vem_local_space_data.PiNabla;
        const Eigen::MatrixXd& pi0km1 = vem_local_space_data.Pi0km1;
        const Eigen::MatrixXd& pi0k = vem_local_space_data.Pi0k;

        result.PiNablaConditioning = LAPACK_utilities::cond(LAPACK_utilities::svd(piNabla));
        result.Pi0kConditioning = LAPACK_utilities::cond(LAPACK_utilities::svd(pi0k));
        result.Pi0km1Conditioning = LAPACK_utilities::cond(LAPACK_utilities::svd(pi0km1));

        const Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(vem_local_space_data.NumProjectorBasisFunctions,
                                                                   vem_local_space_data.NumProjectorBasisFunctions);

        const unsigned int Nkm1 = pi0km1.rows();
        const unsigned int Nk = piNabla.rows();

        result.ErrorPiNabla = (piNabla * vem_local_space_data.Dmatrix - identity).norm() / identity.norm();
        result.ErrorPi0km1 = (pi0km1 * vem_local_space_data.Dmatrix.leftCols(Nkm1) -
                              identity.topLeftCorner(Nkm1, Nkm1)).norm() / identity.topLeftCorner(Nkm1, Nkm1).norm();

        result.ErrorPi0k = (pi0k * vem_local_space_data.Dmatrix - identity).norm() / identity.norm();

        result.ErrorPi0km1Grad.resize(vem_local_space_data.Pi0km1Der.size());
        for(unsigned int d = 0; d < vem_local_space_data.Pi0km1Der.size(); ++d)
        {
            const Eigen::MatrixXd& piDerkm1 = vem_local_space_data.Pi0km1Der[d];
            const Eigen::MatrixXd derMatrix = invDiameter * (vem_local_space_data.Qmatrix *
                                                             vem_monomials.DerivativeMatrix(vem_monomials_data,
                                                                                            d).
                                                             topLeftCorner(Nk, Nkm1) * vem_local_space_data.QmatrixInv.topLeftCorner(Nkm1, Nkm1)).transpose();
            double relErrDenominator = (Nkm1 > 1) ? derMatrix.norm() : 1.0;
            result.ErrorPi0km1Grad[d] = (piDerkm1 * vem_local_space_data.Dmatrix - derMatrix).norm() /
                                        relErrDenominator;
        }

        if (vem_local_space_data.StabMatrix.size() > 0)
        {
            const Eigen::MatrixXd& stabilizationMatrix = vem_local_space_data.StabMatrix;
            result.StabNorm = vem_local_space_data.StabMatrix.norm();
            result.ErrorStabilization = (stabilizationMatrix * vem_local_space_data.Dmatrix).norm();
        }

        if (vem_local_space_data.Hmatrix.size() > 0 && vem_local_space_data.Cmatrix.size() > 0)
            result.ErrorHCD = (vem_local_space_data.Hmatrix - vem_local_space_data.Cmatrix * vem_local_space_data.Dmatrix).norm() / vem_local_space_data.Hmatrix.norm();
        if (vem_local_space_data.Gmatrix.size() > 0 && vem_local_space_data.Bmatrix.size() > 0)
            result.ErrorGBD = (vem_local_space_data.Gmatrix - vem_local_space_data.Bmatrix * vem_local_space_data.Dmatrix).norm() / vem_local_space_data.Gmatrix.norm();

        result.ErrorHED.resize(vem_local_space_data.Pi0km1Der.size(), -1.0);
        for(unsigned int d = 0; d < vem_local_space_data.Pi0km1Der.size(); d++)
        {
            const Eigen::MatrixXd derMatrix = invDiameter *
                                              vem_monomials.DerivativeMatrix(vem_monomials_data,
                                                                             d).
                                              topLeftCorner(Nk,
                                                            vem_local_space_data.Nkm1) *
                                              vem_local_space_data.Hmatrix.topLeftCorner(vem_local_space_data.Nkm1, vem_local_space_data.Nkm1);

            result.ErrorHED[d] = (derMatrix.transpose() - vem_local_space_data.Ematrix[d] * vem_local_space_data.Dmatrix).norm()
                                 / derMatrix.norm();
        }

        return result;
    }
};

}
}
}

#endif
