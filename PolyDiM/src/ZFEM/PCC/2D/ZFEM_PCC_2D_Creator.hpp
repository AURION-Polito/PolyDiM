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

#ifndef __ZFEM_PCC_2D_Creator_HPP
#define __ZFEM_PCC_2D_Creator_HPP

#include "I_ZFEM_PCC_2D_ReferenceElement.hpp"
#include "ZFEM_PCC_2D_ReferenceElement.hpp"

#include "I_ZFEM_PCC_2D_LocalSpace.hpp"
#include "ZFEM_PCC_2D_LocalSpace.hpp"
#include <memory>

namespace Polydim
{
namespace ZFEM
{
namespace PCC
{
enum struct ZFEM_PCC_2D_LocalSpace_Types
{
    ZFEM_PCC_2D_LocalSpace = 1
};

inline std::unique_ptr<I_ZFEM_PCC_2D_ReferenceElement> create_ZFEM_PCC_2D_reference_element(const ZFEM_PCC_2D_LocalSpace_Types &type)
{
    switch (type)
    {
    case ZFEM_PCC_2D_LocalSpace_Types::ZFEM_PCC_2D_LocalSpace:
        return std::make_unique<ZFEM_PCC_2D_ReferenceElement>();
    default:
        throw std::runtime_error("ZFEM type " + std::to_string((unsigned int)type) + " not supported");
    }
}

inline std::unique_ptr<I_ZFEM_PCC_2D_LocalSpace> create_ZFEM_PCC_2D_local_space(const ZFEM_PCC_2D_LocalSpace_Types &type)
{
    switch (type)
    {
    case ZFEM_PCC_2D_LocalSpace_Types::ZFEM_PCC_2D_LocalSpace:
        return std::make_unique<ZFEM_PCC_2D_LocalSpace>();
    default:
        throw std::runtime_error("ZFEM type " + std::to_string((unsigned int)type) + " not supported");
    }
}

} // namespace PCC
} // namespace ZFEM
} // namespace Polydim

#endif
