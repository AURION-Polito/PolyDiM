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

#ifndef __ZFEM_PCC_2D_ReferenceElement_HPP
#define __ZFEM_PCC_2D_ReferenceElement_HPP

#include "I_ZFEM_PCC_2D_ReferenceElement.hpp"
#include "Monomials_2D.hpp"

namespace Polydim
{
namespace ZFEM
{
namespace PCC
{

class ZFEM_PCC_2D_ReferenceElement final : public I_ZFEM_PCC_2D_ReferenceElement
{
  public:
    ZFEM_PCC_2D_ReferenceElement_Data Create(const unsigned int order) const
    {
        ZFEM_PCC_2D_ReferenceElement_Data result;

        result.Dimension = 2;
        result.Order = order;
        result.NumDofs0D = 1;
        result.NumDofs1D = order - 1;
        result.NumDofs2D = (order - 1) * (order - 2) / 2;

        Polydim::FEM::PCC::FEM_Triangle_PCC_2D_ReferenceElement fem_reference_element;
        result.fem_reference_element_data = fem_reference_element.Create(order);

        Utilities::Monomials_2D monomials;
        result.monomials_data = monomials.Compute(order);

        return result;
    }
};
} // namespace PCC
} // namespace ZFEM
} // namespace Polydim

#endif
