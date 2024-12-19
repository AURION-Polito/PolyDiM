#ifndef __VEM_MCC_3D_Velocity_LocalSpace_Creator_HPP
#define __VEM_MCC_3D_Velocity_LocalSpace_Creator_HPP

#include "I_VEM_MCC_3D_Velocity_LocalSpace.hpp"
#include "VEM_MCC_3D_Velocity_LocalSpace.hpp"
#include <memory>


namespace Polydim
{
namespace VEM
{
namespace MCC
{
enum struct VEM_MCC_3D_LocalSpace_Types
{
    VEM_MCC_3D_LocalSpace = 1
};

inline std::unique_ptr<I_VEM_MCC_3D_Velocity_LocalSpace> create_VEM_MCC_3D_local_space(const VEM_MCC_3D_LocalSpace_Types &type)
{
    switch (type)
    {
    case VEM_MCC_3D_LocalSpace_Types::VEM_MCC_3D_LocalSpace:
        return std::make_unique<VEM_MCC_3D_Velocity_LocalSpace>();
    default:
        throw std::runtime_error("VEM type " + std::to_string((unsigned int)type) + " not supported");
    }
}

} // namespace MCC
} // namespace VEM
} // namespace Polydim

#endif
