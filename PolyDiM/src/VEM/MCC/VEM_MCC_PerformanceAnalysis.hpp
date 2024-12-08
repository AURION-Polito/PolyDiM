#ifndef __VEM_MCC_PerformanceAnalysis_H
#define __VEM_MCC_PerformanceAnalysis_H

#include "Eigen/Eigen"
#include <vector>

#include "LAPACK_utilities.hpp"

namespace Polydim
{
namespace VEM
{
namespace MCC
{
struct VEM_MCC_PerformanceAnalysis_Data
{
    double VmatrixConditioning = -1.0; ///< conditioning of piNabla
    double HmatrixConditioning = -1.0; ///< conditioning of piNabla
    double Pi0kConditioning = -1.0; ///< conditioning of piNabla
    double GmatrixConditioning = -1.0; ///< conditioning of piNabla
    double ErrorPi0k = -1.0; ///< |pi0k * Dofs - I|
    double StabNorm = -1.0;
    double ErrorStabilization = -1.0; ///< |S * Dofs|
    double ErrorGBD = -1.0; ///< |G - BD|
};

struct VEM_MCC_PerformanceAnalysis final
{
    template<typename VEM_Monomials_Type,
             typename VEM_Monomials_Data_Type,
             typename VEM_LocalSpace_Type,
             typename VEM_LocalSpaceData_Type>
    VEM_MCC_PerformanceAnalysis_Data Compute(const double& polytopeMeasure,
                                             const double& polytopeDiameter,
                                             const VEM_Monomials_Type& vem_monomials,
                                             const VEM_Monomials_Data_Type& vem_monomials_data,
                                             const VEM_LocalSpace_Type& vem_local_space,
                                             const VEM_LocalSpaceData_Type& vem_local_space_data) const
    {
        VEM_MCC_PerformanceAnalysis_Data result;

        double invDiameter = 1.0 / polytopeDiameter;

        const Eigen::MatrixXd& Vmatrix = vem_local_space_data.Vmatrix;
        const Eigen::MatrixXd& Hmatrix = vem_local_space_data.Hmatrix;
        const Eigen::MatrixXd& Gmatrix = vem_local_space_data.Gmatrix;
        const Eigen::MatrixXd& pi0k = vem_local_space_data.Pi0k;

        result.VmatrixConditioning = LAPACK_utilities::cond(LAPACK_utilities::svd(Vmatrix));
        result.HmatrixConditioning = LAPACK_utilities::cond(LAPACK_utilities::svd(Hmatrix));
        result.Pi0kConditioning = LAPACK_utilities::cond(LAPACK_utilities::svd(pi0k));
        result.GmatrixConditioning = LAPACK_utilities::cond(LAPACK_utilities::svd(Gmatrix));

        const unsigned int Nk = vem_local_space_data.Nk;
        const unsigned int dimension = vem_local_space_data.Dimension;

        const Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(dimension * Nk,
                                                                   dimension * Nk);

        const Eigen::MatrixXd &polynomialBasisDofs = vem_local_space_data.Dmatrix;
        result.ErrorPi0k = (pi0k * polynomialBasisDofs - identity).norm() / identity.norm();

        if (vem_local_space_data.StabMatrix.size() > 0)
        {
            const Eigen::MatrixXd& stabilizationMatrix = vem_local_space_data.StabMatrix;
            result.StabNorm = vem_local_space_data.StabMatrix.norm();
            result.ErrorStabilization = (stabilizationMatrix * polynomialBasisDofs).norm();
        }

        if (vem_local_space_data.Gmatrix.size() > 0 && vem_local_space_data.Bmatrix.size() > 0)
            result.ErrorGBD = (vem_local_space_data.Gmatrix - vem_local_space_data.Bmatrix * polynomialBasisDofs).norm() / vem_local_space_data.Gmatrix.norm();


        return result;
    }
};

}
}
}

#endif
