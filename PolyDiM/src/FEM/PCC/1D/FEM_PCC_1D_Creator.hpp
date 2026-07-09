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

#ifndef __FEM_PCC_1D_Creator_HPP
#define __FEM_PCC_1D_Creator_HPP

#include "FEM_PCC_1D_ReferenceElement.hpp"

#include "FEM_PCC_1D_LocalSpace.hpp"
#include <memory>

namespace Polydim
{
namespace FEM
{
namespace PCC
{

/// @brief Enumeration of the available 1D PCC (primal/continuous) local space types.
///
/// Selects which concrete 1D PCC local space (and its matching reference element) is
/// instantiated by the factory functions below.
enum struct FEM_PCC_1D_LocalSpace_Types
{
    FEM_PCC_1D_LocalSpace = 1 ///< Standard 1D PCC local space.
};

/// @brief Factory building the reference element for a 1D PCC local space type.
/// @param type The 1D PCC local space type to instantiate.
/// @return A unique pointer to the corresponding reference element.
/// @throws std::runtime_error if the requested type is not supported.
inline std::unique_ptr<FEM_PCC_1D_ReferenceElement> create_FEM_PCC_1D_reference_element(const FEM_PCC_1D_LocalSpace_Types &type)
{
    switch (type)
    {
    case FEM_PCC_1D_LocalSpace_Types::FEM_PCC_1D_LocalSpace:
        return std::make_unique<FEM_PCC_1D_ReferenceElement>();
    default:
        throw std::runtime_error("FEM type " + std::to_string((unsigned int)type) + " not supported");
    }
}

/// @brief Factory building the local space for a 1D PCC local space type.
/// @param type The 1D PCC local space type to instantiate.
/// @return A unique pointer to the corresponding local space.
/// @throws std::runtime_error if the requested type is not supported.
inline std::unique_ptr<FEM_PCC_1D_LocalSpace> create_FEM_PCC_1D_local_space(const FEM_PCC_1D_LocalSpace_Types &type)
{
    switch (type)
    {
    case FEM_PCC_1D_LocalSpace_Types::FEM_PCC_1D_LocalSpace:
        return std::make_unique<FEM_PCC_1D_LocalSpace>();
    default:
        throw std::runtime_error("FEM type " + std::to_string((unsigned int)type) + " not supported");
    }
}

} // namespace PCC
} // namespace FEM
} // namespace Polydim

#endif
