#ifndef __VEM_DF_PCC_2D_ReferenceElement_HPP
#define __VEM_DF_PCC_2D_ReferenceElement_HPP

#include "I_VEM_DF_PCC_2D_ReferenceElement.hpp"
#include "VEM_GBasis_2D.hpp"
#include "VEM_Quadrature_2D.hpp"

namespace Polydim
{
namespace VEM
{
namespace DF_PCC
{

class VEM_DF_PCC_2D_Pressure_ReferenceElement final : public I_VEM_DF_PCC_2D_Pressure_ReferenceElement
{
  public:
    VEM_DF_PCC_2D_Pressure_ReferenceElement_Data Create(const unsigned int order) const
    {
        Utilities::VEM_Monomials_2D monomials;
        Quadrature::VEM_Quadrature_2D quadrature;

        VEM_DF_PCC_2D_Pressure_ReferenceElement_Data result;

        result.Monomials = monomials.Compute(order - 1);
        result.Quadrature = quadrature.Compute_DF_PCC_2D(order);

        result.Dimension = 2;
        result.Order = order;
        result.NumDofs0D = 0;
        result.NumDofs1D = 0;
        result.NumDofs2D = order * (order + 1) / 2;

        return result;
    }
};

class VEM_DF_PCC_2D_Velocity_ReferenceElement final : public I_VEM_DF_PCC_2D_Velocity_ReferenceElement
{
  public:
    VEM_DF_PCC_2D_Velocity_ReferenceElement_Data Create(const unsigned int order) const
    {
        Utilities::VEM_GBasis_2D g_basis;
        Utilities::VEM_Monomials_2D monomials;
        Quadrature::VEM_Quadrature_2D quadrature;

        VEM_DF_PCC_2D_Velocity_ReferenceElement_Data result;

        result.Monomials = monomials.Compute(order);
        result.GBasis = g_basis.Compute(order);
        result.Quadrature = quadrature.Compute_DF_PCC_2D(order);

        result.Dimension = 2;
        result.Order = order;
        result.NumDofs0D = 1;
        result.NumDofs1D = order - 1;
        result.NumDofs2D_Divergence = order * (order + 1) / 2 - 1;
        result.NumDofs2D_BigOPlus = (order - 2) * (order - 1) / 2;

        return result;
    }
};
} // namespace DF_PCC
} // namespace VEM
} // namespace Polydim

#endif
