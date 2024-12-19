#ifndef __VEM_DF_PCC_2D_Velocity_LocalSpace_Creator_HPP
#define __VEM_DF_PCC_2D_Velocity_LocalSpace_Creator_HPP

#include "Eigen/Eigen"
#include "I_VEM_DF_PCC_2D_Velocity_LocalSpace.hpp"
#include "VEM_DF_PCC_2D_Velocity_LocalSpace.hpp"
#include <memory>

namespace Polydim
{
namespace VEM
{
namespace DF_PCC
{

enum struct VEM_DF_PCC_2D_LocalSpace_Types
{
    VEM_DF_PCC_2D_LocalSpace = 1,
};

inline std::unique_ptr<I_VEM_DF_PCC_2D_Velocity_LocalSpace> create_VEM_DF_PCC_2D_local_space(const VEM_DF_PCC_2D_LocalSpace_Types &type)
{
    switch (type)
    {
    case VEM_DF_PCC_2D_LocalSpace_Types::VEM_DF_PCC_2D_LocalSpace:
        return std::make_unique<VEM_DF_PCC_2D_Velocity_LocalSpace>();
    default:
        throw std::runtime_error("VEM type " + std::to_string((unsigned int)type) + " not supported");
    }
}

} // namespace DF_PCC
} // namespace VEM
} // namespace Polydim

#endif
