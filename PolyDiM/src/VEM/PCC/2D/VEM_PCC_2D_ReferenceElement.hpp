#ifndef __VEM_PCC_2D_ReferenceElement_H
#define __VEM_PCC_2D_ReferenceElement_H

#include "I_VEM_PCC_2D_ReferenceElement.hpp"
#include "VEM_Monomials_2D.hpp"

namespace Polydim
{
namespace VEM
{
namespace PCC
{

/// \brief Base class for Primal Conforming Virtual Element Method of Constant Degree.
class VEM_PCC_2D_ReferenceElement final : public I_VEM_PCC_2D_ReferenceElement
{
  public:
    VEM_PCC_2D_ReferenceElement_Data Create(const unsigned int order) const
    {
        Monomials::VEM_Monomials_2D monomials;
        Quadrature::VEM_Quadrature_2D quadrature;

        VEM_PCC_2D_ReferenceElement_Data result;

        result.Monomials = monomials.Compute(order);
        result.Quadrature = quadrature.Compute_PCC_2D(order);

        result.Dimension = 2;
        result.Order = order;
        result.NumDofs0D = 1;
        result.NumDofs1D = order - 1;
        result.NumDofs2D = order * (order - 1) / 2;

        return result;
    }
};
} // namespace PCC
} // namespace VEM
} // namespace Polydim

#endif
