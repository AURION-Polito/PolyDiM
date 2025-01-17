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
enum struct FEM_PCC_1D_LocalSpace_Types
{
    FEM_PCC_1D_LocalSpace = 1
};

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
