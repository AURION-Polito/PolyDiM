#ifndef __VEM_PCC_2D_LocalSpace_creator_HPP
#define __VEM_PCC_2D_LocalSpace_creator_HPP

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

std::unique_ptr<I_VEM_PCC_2D_LocalSpace> create_VEM_PCC_2D_local_space(const VEM_PCC_2D_LocalSpace_Types &type)
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
