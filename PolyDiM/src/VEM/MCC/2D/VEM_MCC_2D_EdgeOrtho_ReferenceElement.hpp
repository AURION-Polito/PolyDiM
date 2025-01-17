#ifndef __VEM_MCC_2D_EdgeOrtho_ReferenceElement_H
#define __VEM_MCC_2D_EdgeOrtho_ReferenceElement_H

#include "I_VEM_MCC_2D_ReferenceElement.hpp"
#include "LAPACK_utilities.hpp"
#include "VEM_Monomials_1D.hpp"
#include "VEM_Monomials_2D.hpp"
#include "VEM_Quadrature_2D.hpp"

namespace Polydim
{
namespace VEM
{
namespace MCC
{
/// \brief Base class for Primal Conforming Virtual Element Method of Constant Degree.
class VEM_MCC_2D_EdgeOrtho_Pressure_ReferenceElement final : public I_VEM_MCC_2D_Pressure_ReferenceElement
{
  public:
    VEM_MCC_2D_Pressure_ReferenceElement_Data Create(const unsigned int order) const
    {
        Monomials::VEM_Monomials_2D monomials;
        Quadrature::VEM_Quadrature_2D quadrature;

        VEM_MCC_2D_Pressure_ReferenceElement_Data result;

        result.Monomials = monomials.Compute(order);
        result.Quadrature = quadrature.Compute_MCC_EdgeOrtho_2D(order);

        result.Dimension = 2;
        result.Order = order;
        result.NumDofs0D = 0;
        result.NumDofs1D = 0;
        result.NumDofs2D = (order + 1) * (order + 2) / 2;

        return result;
    }
};

class VEM_MCC_2D_EdgeOrtho_Velocity_ReferenceElement final : public I_VEM_MCC_2D_Velocity_ReferenceElement
{
  public:
    VEM_MCC_2D_Velocity_ReferenceElement_Data Create(const unsigned int order) const
    {
        Monomials::VEM_Monomials_1D monomials_1D;
        Monomials::VEM_Monomials_2D monomials_2D;
        Quadrature::VEM_Quadrature_2D quadrature;

        VEM_MCC_2D_Velocity_ReferenceElement_Data result;

        result.Monomials_1D = monomials_1D.Compute(order + 1);
        result.MonomialsKp1 = monomials_2D.Compute(order + 1);
        result.Quadrature = quadrature.Compute_MCC_EdgeOrtho_2D(order);

        result.Dimension = 2;
        result.Order = order;
        result.NumDofs0D = 0;
        result.NumDofs1D = order + 1;
        result.NumDofs2D = order * (order + 2);

        {
            Gedim::Quadrature::QuadratureData referenceQuadrature1D = result.Quadrature.ReferenceSegmentQuadrature;

            const Eigen::VectorXd sqrtInternalQuadratureWeights1D = referenceQuadrature1D.Weights.array().sqrt();
            const Eigen::MatrixXd VanderBoundary1Dkp1 =
                monomials_1D.Vander(result.Monomials_1D, referenceQuadrature1D.Points, Eigen::Vector3d::Constant(0.0), 1.0);

            Eigen::MatrixXd Q1_1D;
            Eigen::MatrixXd R1_1D;
            LAPACK_utilities::MGS(VanderBoundary1Dkp1, Q1_1D, R1_1D);

            // L2(E)-re-orthogonalization process
            Eigen::MatrixXd Q2_1D;
            Eigen::MatrixXd R2_1D;
            Eigen::MatrixXd temp_1D = sqrtInternalQuadratureWeights1D.asDiagonal() * Q1_1D;
            LAPACK_utilities::MGS(temp_1D, Q2_1D, R2_1D);

            result.edge_ortho.Hmatrix1D = Q2_1D.transpose() * Q2_1D;

            result.edge_ortho.QmatrixInvKp1_1D = (R2_1D * R1_1D).transpose();
            LAPACK_utilities::inverseTri(result.edge_ortho.QmatrixInvKp1_1D, result.edge_ortho.QmatrixKp1_1D, 'L', 'N');
        }

        return result;
    }
};
} // namespace MCC
} // namespace VEM
} // namespace Polydim

#endif
