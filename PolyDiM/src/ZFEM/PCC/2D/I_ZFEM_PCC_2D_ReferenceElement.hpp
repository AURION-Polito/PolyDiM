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

#ifndef __I_ZFEM_PCC_2D_ReferenceElement_HPP
#define __I_ZFEM_PCC_2D_ReferenceElement_HPP

#include "FEM_Triangle_PCC_2D_ReferenceElement.hpp"
#include "MeshUtilities.hpp"
#include "Monomials_Data.hpp"

namespace Polydim
{
namespace ZFEM
{
namespace PCC
{
struct ZFEM_PCC_2D_ReferenceElement_Data final
{
    unsigned int Dimension;
    unsigned int Order;
    unsigned int NumDofs0D;
    unsigned int NumDofs1D;
    unsigned int NumDofs2D;

    Polydim::FEM::PCC::FEM_Triangle_PCC_2D_ReferenceElement_Data fem_reference_element_data;
    Utilities::Monomials_Data monomials_data;

    Gedim::MeshUtilities::MeshGeometricData2DConfig mesh_geometric_data_config = Gedim::MeshUtilities::MeshGeometricData2DConfig(
        {false, false, true, false, true, true, false, true, true, true, true, true, true});
};

class I_ZFEM_PCC_2D_ReferenceElement
{
  public:
    virtual ZFEM_PCC_2D_ReferenceElement_Data Create(const unsigned int order) const = 0;
};

} // namespace PCC
} // namespace ZFEM
} // namespace Polydim

#endif
