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

#ifndef __VEM_MCC_2D_ReferenceElement_HPP
#define __VEM_MCC_2D_ReferenceElement_HPP

#include "I_VEM_MCC_2D_ReferenceElement.hpp"
#include "VEM_Monomials_2D.hpp"
#include "VEM_Quadrature_2D.hpp"

namespace Polydim
{
namespace VEM
{
namespace MCC
{

class VEM_MCC_2D_Pressure_ReferenceElement final : public I_VEM_MCC_2D_Pressure_ReferenceElement
{
  public:
    VEM_MCC_2D_Pressure_ReferenceElement_Data Create(const unsigned int order) const
    {
        Utilities::VEM_Monomials_2D monomials;
        Quadrature::VEM_Quadrature_2D quadrature;

        VEM_MCC_2D_Pressure_ReferenceElement_Data result;

        result.Monomials = monomials.Compute(order);
        result.Quadrature = quadrature.Compute_MCC_2D(order);

        result.Dimension = 2;
        result.Order = order;
        result.NumDofs0D = 0;
        result.NumDofs1D = 0;
        result.NumDofs2D = (order + 1) * (order + 2) / 2;

        return result;
    }
};

class VEM_MCC_2D_Velocity_ReferenceElement final : public I_VEM_MCC_2D_Velocity_ReferenceElement
{
  public:
    VEM_MCC_2D_Velocity_ReferenceElement_Data Create(const unsigned int order) const
    {
        Utilities::VEM_Monomials_2D monomials;
        Quadrature::VEM_Quadrature_2D quadrature;

        VEM_MCC_2D_Velocity_ReferenceElement_Data result;

        result.MonomialsKp1 = monomials.Compute(order + 1);
        result.Quadrature = quadrature.Compute_MCC_2D(order);

        result.Dimension = 2;
        result.Order = order;
        result.NumDofs0D = 0;
        result.NumDofs1D = order + 1;
        result.NumDofs2D = order * (order + 2);

        return result;
    }
};
} // namespace MCC
} // namespace VEM
} // namespace Polydim

#endif
