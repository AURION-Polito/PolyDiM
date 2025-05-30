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

#ifndef __VEM_PCC_2D_Creator_HPP
#define __VEM_PCC_2D_Creator_HPP

#include "I_VEM_PCC_2D_ReferenceElement.hpp"
#include "VEM_PCC_2D_ReferenceElement.hpp"

#include "I_VEM_PCC_2D_LocalSpace.hpp"
#include "VEM_PCC_2D_Inertia_LocalSpace.hpp"
#include "VEM_PCC_2D_LocalSpace.hpp"
#include "VEM_PCC_2D_Ortho_LocalSpace.hpp"
#include <memory>

namespace Polydim
{
namespace VEM
{
namespace PCC
{
enum struct VEM_PCC_2D_LocalSpace_Types
{
    VEM_PCC_2D_LocalSpace = 1,
    VEM_PCC_2D_Inertia_LocalSpace = 2,
    VEM_PCC_2D_Ortho_LocalSpace = 3
};

inline std::unique_ptr<I_VEM_PCC_2D_ReferenceElement> create_VEM_PCC_2D_reference_element(const VEM_PCC_2D_LocalSpace_Types &type)
{
    switch (type)
    {
    case VEM_PCC_2D_LocalSpace_Types::VEM_PCC_2D_LocalSpace:
    case VEM_PCC_2D_LocalSpace_Types::VEM_PCC_2D_Inertia_LocalSpace:
    case VEM_PCC_2D_LocalSpace_Types::VEM_PCC_2D_Ortho_LocalSpace:
        return std::make_unique<VEM_PCC_2D_ReferenceElement>();
    default:
        throw std::runtime_error("VEM type " + std::to_string((unsigned int)type) + " not supported");
    }
}

inline std::unique_ptr<I_VEM_PCC_2D_LocalSpace> create_VEM_PCC_2D_local_space(const VEM_PCC_2D_LocalSpace_Types &type)
{
    switch (type)
    {
    case VEM_PCC_2D_LocalSpace_Types::VEM_PCC_2D_LocalSpace:
        return std::make_unique<VEM_PCC_2D_LocalSpace>();
    case VEM_PCC_2D_LocalSpace_Types::VEM_PCC_2D_Inertia_LocalSpace:
        return std::make_unique<VEM_PCC_2D_Inertia_LocalSpace>();
    case VEM_PCC_2D_LocalSpace_Types::VEM_PCC_2D_Ortho_LocalSpace:
        return std::make_unique<VEM_PCC_2D_Ortho_LocalSpace>();
    default:
        throw std::runtime_error("VEM type " + std::to_string((unsigned int)type) + " not supported");
    }
}

} // namespace PCC
} // namespace VEM
} // namespace Polydim

#endif
