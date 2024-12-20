#ifndef __VEM_MCC_3D_ReferenceElement_H
#define __VEM_MCC_3D_ReferenceElement_H

#include "I_VEM_MCC_3D_ReferenceElement.hpp"
#include "VEM_Monomials_2D.hpp"
#include "VEM_Monomials_3D.hpp"
#include "VEM_Quadrature_3D.hpp"

namespace Polydim
{
namespace VEM
{
namespace MCC
{


/// \brief Base class for Primal Conforming Virtual Element Method of Constant Degree.
class VEM_MCC_3D_Pressure_ReferenceElement final : public I_VEM_MCC_3D_Pressure_ReferenceElement
{
  public:
    VEM_MCC_3D_Pressure_ReferenceElement_Data Create(const unsigned int order) const
    {
        Monomials::VEM_Monomials_3D monomials;
        Quadrature::VEM_Quadrature_3D quadrature;

        VEM_MCC_3D_Pressure_ReferenceElement_Data result;

        result.Monomials = monomials.Compute(order);
        result.Quadrature = quadrature.Compute_MCC_3D(order);

        result.Dimension = 3;
        result.Order = order;
        result.NumDofs0D = 0;
        result.NumDofs1D = 0;
        result.NumDofs2D = 0;
        result.NumDofs3D = (order + 1) * (order + 2) * (order + 3) / 6;

        return result;
    }
};



class VEM_MCC_3D_Velocity_ReferenceElement final : public I_VEM_MCC_3D_Velocity_ReferenceElement
{
  public:
    VEM_MCC_3D_Velocity_ReferenceElement_Data Create(const unsigned int order) const
    {
        Monomials::VEM_Monomials_2D monomials2D;
        Monomials::VEM_Monomials_3D monomials3D;
        Quadrature::VEM_Quadrature_3D quadrature;

        VEM_MCC_3D_Velocity_ReferenceElement_Data result;

        result.Monomials2D = monomials2D.Compute(order);
        result.MonomialsKp1 = monomials3D.Compute(order + 1);
        result.Quadrature = quadrature.Compute_MCC_3D(order);

        result.Dimension = 3;
        result.Order = order;
        result.NumDofs0D = 0;
        result.NumDofs1D = 0;
        result.NumDofs2D = (order + 1) * (order + 2) / 2;
        result.NumDofs3D = order * (order + 2) * (order + 3) / 2;

        return result;
    }
};
} // namespace MCC
} // namespace VEM
} // namespace Polydim

#endif
