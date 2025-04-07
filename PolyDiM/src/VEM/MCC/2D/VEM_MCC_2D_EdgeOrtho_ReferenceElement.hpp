#ifndef __VEM_MCC_2D_EdgeOrtho_ReferenceElement_HPP
#define __VEM_MCC_2D_EdgeOrtho_ReferenceElement_HPP

#include "I_VEM_MCC_2D_ReferenceElement.hpp"
#include "VEM_Monomials_1D.hpp"
#include "VEM_Monomials_2D.hpp"
#include "VEM_Quadrature_2D.hpp"

namespace Polydim
{
namespace VEM
{
namespace MCC
{
class VEM_MCC_2D_EdgeOrtho_Pressure_ReferenceElement final : public I_VEM_MCC_2D_Pressure_ReferenceElement
{
  public:
    VEM_MCC_2D_Pressure_ReferenceElement_Data Create(const unsigned int order) const
    {
        Utilities::VEM_Monomials_2D monomials;
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
        Utilities::VEM_Monomials_1D monomials_1D;
        Utilities::VEM_Monomials_2D monomials_2D;
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
            const Eigen::MatrixXd VanderBoundary1Dkp1 =
                monomials_1D.Vander(result.Monomials_1D, referenceQuadrature1D.Points, Eigen::Vector3d::Constant(0.0), 1.0);

            monomials_1D.MGSOrthonormalize(referenceQuadrature1D.Weights,
                                           VanderBoundary1Dkp1,
                                           result.edge_ortho.Hmatrix1D,
                                           result.edge_ortho.QmatrixInvKp1_1D,
                                           result.edge_ortho.QmatrixKp1_1D);
        }

        return result;
    }
};
} // namespace MCC
} // namespace VEM
} // namespace Polydim

#endif
