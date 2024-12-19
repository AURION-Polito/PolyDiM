#ifndef __VEM_PCC_3D_ReferenceElement_H
#define __VEM_PCC_3D_ReferenceElement_H

#include "VEM_Monomials_3D.hpp"
#include "I_VEM_PCC_3D_ReferenceElement.hpp"

namespace Polydim
{
namespace VEM
{
namespace PCC
{

/// \brief Base class for Primal Conforming Virtual Element Method of Constant Degree.
class VEM_PCC_3D_ReferenceElement final : public I_VEM_PCC_3D_ReferenceElement
{
  public:
    VEM_PCC_3D_ReferenceElement_Data Create(const unsigned int order) const
    {
        Monomials::VEM_Monomials_3D monomials;
        Quadrature::VEM_Quadrature_3D quadrature;

        VEM_PCC_3D_ReferenceElement_Data result;

        result.Monomials = monomials.Compute(order);
        result.Quadrature = quadrature.Compute_PCC_3D(order);

        result.Dimension = 3;
        result.Order = order;
        result.NumDofs0D = 1;
        result.NumDofs1D = order - 1;
        result.NumDofs2D = order * (order - 1) / 2;
        result.NumDofs3D = (order + 1) * order * (order - 1) / 6;

        return result;
    }
};
} // namespace PCC
} // namespace VEM
} // namespace Polydim

#endif
