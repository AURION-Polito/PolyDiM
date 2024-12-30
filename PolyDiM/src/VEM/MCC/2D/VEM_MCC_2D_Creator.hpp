#ifndef __VEM_MCC_2D_Creator_HPP
#define __VEM_MCC_2D_Creator_HPP

#include "I_VEM_MCC_2D_ReferenceElement.hpp"
#include "I_VEM_MCC_2D_Velocity_LocalSpace.hpp"
#include "VEM_MCC_2D_EdgeOrtho_ReferenceElement.hpp"
#include "VEM_MCC_2D_EdgeOrtho_Velocity_LocalSpace.hpp"
#include "VEM_MCC_2D_Ortho_EdgeOrtho_Velocity_LocalSpace.hpp"
#include "VEM_MCC_2D_Ortho_Velocity_LocalSpace.hpp"
#include "VEM_MCC_2D_Partial_Velocity_LocalSpace.hpp"
#include "VEM_MCC_2D_Pressure_LocalSpace.hpp"
#include "VEM_MCC_2D_ReferenceElement.hpp"
#include "VEM_MCC_2D_Velocity_LocalSpace.hpp"
#include <memory>

namespace Polydim
{
namespace VEM
{
namespace MCC
{

enum struct VEM_MCC_2D_LocalSpace_Types
{
    VEM_MCC_2D_LocalSpace = 1,
    VEM_MCC_2D_Partial_LocalSpace = 2,
    VEM_MCC_2D_Ortho_LocalSpace = 3,
    VEM_MCC_2D_EdgeOrtho_LocalSpace = 4,
    VEM_MCC_2D_Ortho_EdgeOrtho_LocalSpace = 5
};

inline std::unique_ptr<I_VEM_MCC_2D_Pressure_ReferenceElement> create_VEM_MCC_2D_pressure_reference_element(const VEM_MCC_2D_LocalSpace_Types &type)
{
    switch (type)
    {
    case VEM_MCC_2D_LocalSpace_Types::VEM_MCC_2D_LocalSpace:
    case VEM_MCC_2D_LocalSpace_Types::VEM_MCC_2D_Partial_LocalSpace:
    case VEM_MCC_2D_LocalSpace_Types::VEM_MCC_2D_Ortho_LocalSpace:
        return std::make_unique<VEM_MCC_2D_Pressure_ReferenceElement>();
    case VEM_MCC_2D_LocalSpace_Types::VEM_MCC_2D_EdgeOrtho_LocalSpace:
    case VEM_MCC_2D_LocalSpace_Types::VEM_MCC_2D_Ortho_EdgeOrtho_LocalSpace:
        return std::make_unique<VEM_MCC_2D_EdgeOrtho_Pressure_ReferenceElement>();
    default:
        throw std::runtime_error("VEM type " + std::to_string((unsigned int)type) + " not supported");
    }
}

inline std::unique_ptr<I_VEM_MCC_2D_Velocity_ReferenceElement> create_VEM_MCC_2D_velocity_reference_element(const VEM_MCC_2D_LocalSpace_Types &type)
{
    switch (type)
    {
    case VEM_MCC_2D_LocalSpace_Types::VEM_MCC_2D_LocalSpace:
    case VEM_MCC_2D_LocalSpace_Types::VEM_MCC_2D_Partial_LocalSpace:
    case VEM_MCC_2D_LocalSpace_Types::VEM_MCC_2D_Ortho_LocalSpace:
        return std::make_unique<VEM_MCC_2D_Velocity_ReferenceElement>();
    case VEM_MCC_2D_LocalSpace_Types::VEM_MCC_2D_EdgeOrtho_LocalSpace:
    case VEM_MCC_2D_LocalSpace_Types::VEM_MCC_2D_Ortho_EdgeOrtho_LocalSpace:
        return std::make_unique<VEM_MCC_2D_EdgeOrtho_Velocity_ReferenceElement>();
    default:
        throw std::runtime_error("VEM type " + std::to_string((unsigned int)type) + " not supported");
    }
}

inline std::unique_ptr<I_VEM_MCC_2D_Velocity_LocalSpace> create_VEM_MCC_2D_pressure_local_space(const VEM_MCC_2D_LocalSpace_Types &type)
{
    switch (type)
    {
    case VEM_MCC_2D_LocalSpace_Types::VEM_MCC_2D_LocalSpace:
    case VEM_MCC_2D_LocalSpace_Types::VEM_MCC_2D_Partial_LocalSpace:
    case VEM_MCC_2D_LocalSpace_Types::VEM_MCC_2D_Ortho_LocalSpace:
    case VEM_MCC_2D_LocalSpace_Types::VEM_MCC_2D_EdgeOrtho_LocalSpace:
    case VEM_MCC_2D_LocalSpace_Types::VEM_MCC_2D_Ortho_EdgeOrtho_LocalSpace:
        return std::make_unique<VEM_MCC_2D_Pressure_LocalSpace>();
    default:
        throw std::runtime_error("VEM type " + std::to_string((unsigned int)type) + " not supported");
    }
}

inline std::unique_ptr<I_VEM_MCC_2D_Velocity_LocalSpace> create_VEM_MCC_2D_velocity_local_space(const VEM_MCC_2D_LocalSpace_Types &type)
{
    switch (type)
    {
    case VEM_MCC_2D_LocalSpace_Types::VEM_MCC_2D_LocalSpace:
        return std::make_unique<VEM_MCC_2D_Velocity_LocalSpace>();
    case VEM_MCC_2D_LocalSpace_Types::VEM_MCC_2D_Partial_LocalSpace:
        return std::make_unique<VEM_MCC_2D_Partial_Velocity_LocalSpace>();
    case VEM_MCC_2D_LocalSpace_Types::VEM_MCC_2D_Ortho_LocalSpace:
        return std::make_unique<VEM_MCC_2D_Ortho_Velocity_LocalSpace>();
    case VEM_MCC_2D_LocalSpace_Types::VEM_MCC_2D_EdgeOrtho_LocalSpace:
        return std::make_unique<VEM_MCC_2D_EdgeOrtho_Velocity_LocalSpace>();
    case VEM_MCC_2D_LocalSpace_Types::VEM_MCC_2D_Ortho_EdgeOrtho_LocalSpace:
        return std::make_unique<VEM_MCC_2D_Ortho_EdgeOrtho_Velocity_LocalSpace>();
    default:
        throw std::runtime_error("VEM type " + std::to_string((unsigned int)type) + " not supported");
    }
}

} // namespace MCC
} // namespace VEM
} // namespace Polydim

#endif
