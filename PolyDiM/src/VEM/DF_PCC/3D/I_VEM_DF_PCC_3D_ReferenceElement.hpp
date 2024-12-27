#ifndef __I_VEM_DF_PCC_3D_ReferenceElement_H
#define __I_VEM_DF_PCC_3D_ReferenceElement_H

#include "VEM_GBasis_3D.hpp"
#include "VEM_Quadrature_3D.hpp"

namespace Polydim
{
namespace VEM
{
namespace DF_PCC
{
struct VEM_DF_PCC_3D_Pressure_ReferenceElement_Data final
{
    unsigned int Dimension;
    unsigned int Order;
    unsigned int NumDofs0D; ///< Number of dofs for each vertex.
    unsigned int NumDofs1D; ///< Number of dofs internal to each edge.
    unsigned int NumDofs2D; ///< Number of dofs internal to each polygon.
    unsigned int NumDofs3D; ///< Number of dofs internal to each polygon.

    Monomials::VEM_Monomials_Data Monomials;
    Quadrature::VEM_QuadratureData_3D Quadrature;
};

class I_VEM_DF_PCC_3D_Pressure_ReferenceElement
{
  public:
    virtual VEM_DF_PCC_3D_Pressure_ReferenceElement_Data Create(const unsigned int order) const = 0;
};

struct VEM_DF_PCC_3D_Velocity_ReferenceElement_Data final
{
    unsigned int Dimension;
    unsigned int Order;
    unsigned int NumDofs0D;
    unsigned int NumDofs1D;
    unsigned int NumDofs2D;
    unsigned int NumDofs3D_BigOPlus;
    unsigned int NumDofs3D_Divergence;

    Monomials::VEM_Monomials_Data Monomials;
    Monomials::VEM_GBasis_Data GBasis;
    Quadrature::VEM_QuadratureData_3D Quadrature;
};

class I_VEM_DF_PCC_3D_Velocity_ReferenceElement
{
  public:
    virtual VEM_DF_PCC_3D_Velocity_ReferenceElement_Data Create(const unsigned int order) const = 0;
};
} // namespace DF_PCC
} // namespace VEM
} // namespace Polydim

#endif
