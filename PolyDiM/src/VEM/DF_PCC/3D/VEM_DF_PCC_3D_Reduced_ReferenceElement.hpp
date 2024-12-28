#ifndef __VEM_DF_PCC_3D_Reduced_ReferenceElement_H
#define __VEM_DF_PCC_3D_Reduced_ReferenceElement_H

#include "I_VEM_DF_PCC_3D_ReferenceElement.hpp"
#include "VEM_GBasis_3D.hpp"
#include "VEM_Quadrature_3D.hpp"

namespace Polydim
{
namespace VEM
{
namespace DF_PCC
{

class VEM_DF_PCC_3D_Reduced_Pressure_ReferenceElement final : public I_VEM_DF_PCC_3D_Pressure_ReferenceElement
{
  public:
    VEM_DF_PCC_3D_Pressure_ReferenceElement_Data Create(const unsigned int order) const
    {
        Monomials::VEM_Monomials_3D monomials;
        Quadrature::VEM_Quadrature_3D quadrature;

        VEM_DF_PCC_3D_Pressure_ReferenceElement_Data result;

        result.Monomials = monomials.Compute(order - 1);
        result.Quadrature = quadrature.Compute_DF_PCC_3D(order);

        result.Dimension = 3;
        result.Order = order;
        result.NumDofs0D = 0;
        result.NumDofs1D = 0;
        result.NumDofs2D = 0;
        result.NumDofs3D = 1;

        return result;
    }
};

class VEM_DF_PCC_3D_Reduced_Velocity_ReferenceElement final : public I_VEM_DF_PCC_3D_Velocity_ReferenceElement
{
  public:
    VEM_DF_PCC_3D_Velocity_ReferenceElement_Data Create(const unsigned int order) const
    {
        Monomials::VEM_GBasis_3D g_basis;
        Monomials::VEM_Monomials_3D monomials;
        Quadrature::VEM_Quadrature_3D quadrature;

        VEM_DF_PCC_3D_Velocity_ReferenceElement_Data result;

        result.Monomials = monomials.Compute(order);
        result.GBasis = g_basis.Compute(order);
        result.Quadrature = quadrature.Compute_DF_PCC_3D(order);

        result.Dimension = 3;
        result.Order = order;
        result.NumDofs0D = 1;
        result.NumDofs1D = order - 1;
        result.NumDofs2D = order * (order - 1) / 2;
        result.NumDofs3D_Divergence = 0;
        result.NumDofs3D_BigOPlus = order * (order + 1) * (order - 1) / 2 - result.NumDofs3D_Divergence;

        return result;
    }
};
} // namespace DF_PCC
} // namespace VEM
} // namespace Polydim

#endif
