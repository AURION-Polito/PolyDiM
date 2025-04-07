#ifndef __I_VEM_PCC_3D_ReferenceElement_HPP
#define __I_VEM_PCC_3D_ReferenceElement_HPP

#include "VEM_Monomials_Data.hpp"
#include "VEM_Quadrature_3D.hpp"

namespace Polydim
{
namespace VEM
{
namespace PCC
{
struct VEM_PCC_3D_ReferenceElement_Data final
{
    unsigned int Dimension;
    unsigned int Order;
    unsigned int NumDofs0D;
    unsigned int NumDofs1D;
    unsigned int NumDofs2D;
    unsigned int NumDofs3D;

    Utilities::VEM_Monomials_Data Monomials;
    Quadrature::VEM_QuadratureData_3D Quadrature;
};

class I_VEM_PCC_3D_ReferenceElement
{
  public:
    virtual VEM_PCC_3D_ReferenceElement_Data Create(const unsigned int order) const = 0;
};

} // namespace PCC
} // namespace VEM
} // namespace Polydim

#endif
