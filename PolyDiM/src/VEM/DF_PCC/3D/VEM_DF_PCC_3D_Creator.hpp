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

#ifndef __VEM_DF_PCC_3D_Creator_HPP
#define __VEM_DF_PCC_3D_Creator_HPP

#include "VEM_DF_PCC_3D_Pressure_LocalSpace.hpp"
#include "VEM_DF_PCC_3D_Reduced_Pressure_LocalSpace.hpp"
#include "VEM_DF_PCC_3D_Reduced_ReferenceElement.hpp"
#include "VEM_DF_PCC_3D_Reduced_Velocity_LocalSpace.hpp"
#include "VEM_DF_PCC_3D_ReferenceElement.hpp"
#include "VEM_DF_PCC_3D_Velocity_LocalSpace.hpp"
#include "VEM_PCC_2D_ReferenceElement.hpp"
#include <memory>

namespace Polydim
{
namespace VEM
{
namespace DF_PCC
{
enum struct VEM_DF_PCC_3D_LocalSpace_Types
{
    VEM_DF_PCC_3D_LocalSpace = 1,
    VEM_DF_PCC_3D_Reduced_LocalSpace = 2
};

inline std::unique_ptr<PCC::I_VEM_PCC_2D_ReferenceElement> create_VEM_DF_PCC_3D_velocity_reference_element_2D(const VEM_DF_PCC_3D_LocalSpace_Types &type)
{
    switch (type)
    {
    case VEM_DF_PCC_3D_LocalSpace_Types::VEM_DF_PCC_3D_LocalSpace:
    case VEM_DF_PCC_3D_LocalSpace_Types::VEM_DF_PCC_3D_Reduced_LocalSpace:
        return std::make_unique<PCC::VEM_PCC_2D_ReferenceElement>();
    default:
        throw std::runtime_error("VEM type " + std::to_string((unsigned int)type) + " not supported");
    }
}

inline std::unique_ptr<I_VEM_DF_PCC_3D_Pressure_ReferenceElement> create_VEM_DF_PCC_3D_full_pressure_reference_element_3D(
    const VEM_DF_PCC_3D_LocalSpace_Types &type)
{
    switch (type)
    {
    case VEM_DF_PCC_3D_LocalSpace_Types::VEM_DF_PCC_3D_Reduced_LocalSpace:
        return std::make_unique<VEM_DF_PCC_3D_Pressure_ReferenceElement>();
    default:
        throw std::runtime_error("VEM type " + std::to_string((unsigned int)type) + " not supported");
    }
}

inline std::unique_ptr<I_VEM_DF_PCC_3D_Velocity_ReferenceElement> create_VEM_DF_PCC_3D_full_velocity_reference_element_3D(
    const VEM_DF_PCC_3D_LocalSpace_Types &type)
{
    switch (type)
    {
    case VEM_DF_PCC_3D_LocalSpace_Types::VEM_DF_PCC_3D_Reduced_LocalSpace:
        return std::make_unique<VEM_DF_PCC_3D_Velocity_ReferenceElement>();
    default:
        throw std::runtime_error("VEM type " + std::to_string((unsigned int)type) + " not supported");
    }
}

inline std::unique_ptr<I_VEM_DF_PCC_3D_Pressure_LocalSpace> create_VEM_DF_PCC_3D_full_pressure_local_space_3D(const VEM_DF_PCC_3D_LocalSpace_Types &type)
{
    switch (type)
    {
    case VEM_DF_PCC_3D_LocalSpace_Types::VEM_DF_PCC_3D_Reduced_LocalSpace:
        return std::make_unique<VEM_DF_PCC_3D_Pressure_LocalSpace>();
    default:
        throw std::runtime_error("VEM type " + std::to_string((unsigned int)type) + " not supported");
    }
}

inline std::unique_ptr<I_VEM_DF_PCC_3D_Velocity_LocalSpace> create_VEM_DF_PCC_3D_full_velocity_local_space_3D(const VEM_DF_PCC_3D_LocalSpace_Types &type)
{
    switch (type)
    {
    case VEM_DF_PCC_3D_LocalSpace_Types::VEM_DF_PCC_3D_Reduced_LocalSpace:
        return std::make_unique<VEM_DF_PCC_3D_Velocity_LocalSpace>();
    default:
        throw std::runtime_error("VEM type " + std::to_string((unsigned int)type) + " not supported");
    }
}

inline std::unique_ptr<I_VEM_DF_PCC_3D_Pressure_ReferenceElement> create_VEM_DF_PCC_3D_pressure_reference_element_3D(const VEM_DF_PCC_3D_LocalSpace_Types &type)
{
    switch (type)
    {
    case VEM_DF_PCC_3D_LocalSpace_Types::VEM_DF_PCC_3D_LocalSpace:
        return std::make_unique<VEM_DF_PCC_3D_Pressure_ReferenceElement>();
    case VEM_DF_PCC_3D_LocalSpace_Types::VEM_DF_PCC_3D_Reduced_LocalSpace:
        return std::make_unique<VEM_DF_PCC_3D_Reduced_Pressure_ReferenceElement>();
    default:
        throw std::runtime_error("VEM type " + std::to_string((unsigned int)type) + " not supported");
    }
}

inline std::unique_ptr<I_VEM_DF_PCC_3D_Velocity_ReferenceElement> create_VEM_DF_PCC_3D_velocity_reference_element_3D(const VEM_DF_PCC_3D_LocalSpace_Types &type)
{
    switch (type)
    {
    case VEM_DF_PCC_3D_LocalSpace_Types::VEM_DF_PCC_3D_LocalSpace:
        return std::make_unique<VEM_DF_PCC_3D_Velocity_ReferenceElement>();
    case VEM_DF_PCC_3D_LocalSpace_Types::VEM_DF_PCC_3D_Reduced_LocalSpace:
        return std::make_unique<VEM_DF_PCC_3D_Reduced_Velocity_ReferenceElement>();
    default:
        throw std::runtime_error("VEM type " + std::to_string((unsigned int)type) + " not supported");
    }
}

inline std::unique_ptr<I_VEM_DF_PCC_3D_Pressure_LocalSpace> create_VEM_DF_PCC_3D_pressure_local_space_3D(const VEM_DF_PCC_3D_LocalSpace_Types &type)
{
    switch (type)
    {
    case VEM_DF_PCC_3D_LocalSpace_Types::VEM_DF_PCC_3D_LocalSpace:
        return std::make_unique<VEM_DF_PCC_3D_Pressure_LocalSpace>();
    case VEM_DF_PCC_3D_LocalSpace_Types::VEM_DF_PCC_3D_Reduced_LocalSpace:
        return std::make_unique<VEM_DF_PCC_3D_Reduced_Pressure_LocalSpace>();
    default:
        throw std::runtime_error("VEM type " + std::to_string((unsigned int)type) + " not supported");
    }
}

inline std::unique_ptr<I_VEM_DF_PCC_3D_Velocity_LocalSpace> create_VEM_DF_PCC_3D_velocity_local_space_3D(const VEM_DF_PCC_3D_LocalSpace_Types &type)
{
    switch (type)
    {
    case VEM_DF_PCC_3D_LocalSpace_Types::VEM_DF_PCC_3D_LocalSpace:
        return std::make_unique<VEM_DF_PCC_3D_Velocity_LocalSpace>();
    case VEM_DF_PCC_3D_LocalSpace_Types::VEM_DF_PCC_3D_Reduced_LocalSpace:
        return std::make_unique<VEM_DF_PCC_3D_Reduced_Velocity_LocalSpace>();
    default:
        throw std::runtime_error("VEM type " + std::to_string((unsigned int)type) + " not supported");
    }
}

} // namespace DF_PCC
} // namespace VEM
} // namespace Polydim

#endif
