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

#ifndef __I_VEM_DF_PCC_2D_ReferenceElement_HPP
#define __I_VEM_DF_PCC_2D_ReferenceElement_HPP

#include "VEM_GBasis_Data.hpp"
#include "VEM_Quadrature_2D.hpp"

namespace Polydim
{
namespace VEM
{
namespace DF_PCC
{
struct VEM_DF_PCC_2D_Pressure_ReferenceElement_Data final
{
    unsigned int Dimension;
    unsigned int Order;
    unsigned int NumDofs0D;
    unsigned int NumDofs1D;
    unsigned int NumDofs2D;

    Utilities::VEM_Monomials_Data Monomials;
    Quadrature::VEM_QuadratureData_2D Quadrature;
};

class I_VEM_DF_PCC_2D_Pressure_ReferenceElement
{
  public:
    virtual VEM_DF_PCC_2D_Pressure_ReferenceElement_Data Create(const unsigned int order) const = 0;
};

struct VEM_DF_PCC_2D_Velocity_ReferenceElement_Data final
{
    unsigned int Dimension;
    unsigned int Order;
    unsigned int NumDofs0D;
    unsigned int NumDofs1D;
    unsigned int NumDofs2D_BigOPlus;
    unsigned int NumDofs2D_Divergence;

    Utilities::VEM_Monomials_Data Monomials;
    Utilities::VEM_GBasis_Data GBasis;
    Quadrature::VEM_QuadratureData_2D Quadrature;
};

class I_VEM_DF_PCC_2D_Velocity_ReferenceElement
{
  public:
    virtual VEM_DF_PCC_2D_Velocity_ReferenceElement_Data Create(const unsigned int order) const = 0;
};
} // namespace DF_PCC
} // namespace VEM
} // namespace Polydim

#endif
