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

#ifndef __VEM_DF_PCC_2D_Reduced_ReferenceElement_HPP
#define __VEM_DF_PCC_2D_Reduced_ReferenceElement_HPP

#include "GBasis_2D.hpp"
#include "I_VEM_DF_PCC_2D_ReferenceElement.hpp"
#include "VEM_Quadrature_2D.hpp"

namespace Polydim
{
namespace VEM
{
namespace DF_PCC
{

class VEM_DF_PCC_2D_Reduced_Pressure_ReferenceElement final : public Polydim::VEM::DF_PCC::I_VEM_DF_PCC_2D_Pressure_ReferenceElement
{
  public:
    VEM_DF_PCC_2D_Pressure_ReferenceElement_Data Create(const unsigned int order) const
    {
        Polydim::Utilities::Monomials_2D monomials;
        Polydim::VEM::Quadrature::VEM_Quadrature_2D quadrature;

        Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Pressure_ReferenceElement_Data result;

        result.Monomials = monomials.Compute(0);
        result.Quadrature = quadrature.Compute_DF_PCC_2D(order);

        result.Dimension = 2;
        result.Order = order;
        result.NumDofs0D = 0;
        result.NumDofs1D = 0;
        result.NumDofs2D = 1;

        return result;
    }
};

class VEM_DF_PCC_2D_Reduced_Velocity_ReferenceElement final : public Polydim::VEM::DF_PCC::I_VEM_DF_PCC_2D_Velocity_ReferenceElement
{
  public:
    VEM_DF_PCC_2D_Velocity_ReferenceElement_Data Create(const unsigned int order) const
    {
        Polydim::Utilities::GBasis_2D g_basis;
        Polydim::Utilities::Monomials_2D monomials;
        Polydim::VEM::Quadrature::VEM_Quadrature_2D quadrature;

        Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_ReferenceElement_Data result;

        result.Monomials = monomials.Compute(order);
        result.GBasis = g_basis.Compute(order);
        result.Quadrature = quadrature.Compute_DF_PCC_2D(order);

        result.Dimension = 2;
        result.Order = order;
        result.NumDofs0D = 1;
        result.NumDofs1D = order - 1;
        result.NumDofs2D_Divergence = 0;
        result.NumDofs2D_BigOPlus = (order - 2) * (order - 1) / 2;

        return result;
    }
};
} // namespace DF_PCC
} // namespace VEM
} // namespace Polydim

#endif
