#ifndef __I_VEM_PCC_2D_ReferenceElement_H
#define __I_VEM_PCC_2D_ReferenceElement_H

#include "VEM_Monomials_Data.hpp"
#include "VEM_Quadrature_2D.hpp"

namespace Polydim
{
namespace VEM
{
namespace PCC
{
/// \brief Base class for storing information related to \ref VEM::PCC::I_VEM_PCC_2D_ReferenceElement
struct VEM_PCC_2D_ReferenceElement_Data final
{
    unsigned int Dimension; ///< Geometric dimension
    unsigned int Order;     ///< Order of the method
    unsigned int NumDofs0D; ///< Number of dofs for each vertex.
    unsigned int NumDofs1D; ///< Number of dofs internal to each edge.
    unsigned int NumDofs2D; ///< Number of dofs internal to each polygon.

    Monomials::VEM_Monomials_Data Monomials;      ///< Monomials used as support for building vem local matrices
    Quadrature::VEM_QuadratureData_2D Quadrature; ///< Quadrature used as support for building vem local matrices
};

class I_VEM_PCC_2D_ReferenceElement
{
  public:
    virtual VEM_PCC_2D_ReferenceElement_Data Create(const unsigned int order) const = 0;
};

} // namespace PCC
} // namespace VEM
} // namespace Polydim

#endif
