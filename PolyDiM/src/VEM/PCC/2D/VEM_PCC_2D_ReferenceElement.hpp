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

#ifndef __VEM_PCC_2D_ReferenceElement_HPP
#define __VEM_PCC_2D_ReferenceElement_HPP

#include "I_VEM_PCC_2D_ReferenceElement.hpp"
#include "VEM_Monomials_2D.hpp"

namespace Polydim
{
namespace VEM
{
namespace PCC
{

class VEM_PCC_2D_ReferenceElement final : public Polydim::VEM::PCC::I_VEM_PCC_2D_ReferenceElement
{
  public:
    VEM_PCC_2D_ReferenceElement_Data Create(const unsigned int order) const
    {
        Polydim::VEM::Utilities::VEM_Monomials_2D monomials;
        Polydim::VEM::Quadrature::VEM_Quadrature_2D quadrature;

        Polydim::VEM::PCC::VEM_PCC_2D_ReferenceElement_Data result;

        result.Monomials = monomials.Compute(order);
        result.Quadrature = quadrature.Compute_PCC_2D(order);

        result.Dimension = 2;
        result.Order = order;
        result.NumDofs0D = 1;
        result.NumDofs1D = order - 1;
        result.NumDofs2D = order * (order - 1) / 2;

        return result;
    }
};
} // namespace PCC
} // namespace VEM
} // namespace Polydim

#endif
