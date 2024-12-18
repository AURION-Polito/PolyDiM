#ifndef __VEM_MCC_2D_ReferenceElement_H
#define __VEM_MCC_2D_ReferenceElement_H

#include "VEM_Monomials_2D.hpp"
#include "VEM_Quadrature_2D.hpp"

namespace Polydim
{
namespace VEM
{
namespace MCC
{
struct VEM_MCC_2D_Pressure_ReferenceElement_Data final
{
    unsigned int Dimension;
    unsigned int Order;
    unsigned int NumDofs0D; ///< Number of dofs for each vertex.
    unsigned int NumDofs1D; ///< Number of dofs internal to each edge.
    unsigned int NumDofs2D; ///< Number of dofs internal to each polygon.

    Monomials::VEM_Monomials_Data Monomials;
    Quadrature::VEM_QuadratureData_2D Quadrature;
};

/// \brief Base class for Primal Conforming Virtual Element Method of Constant Degree.
class VEM_MCC_2D_Pressure_ReferenceElement final
{
  public:
    VEM_MCC_2D_Pressure_ReferenceElement_Data Create(const unsigned int order) const
    {
        Monomials::VEM_Monomials_2D monomials;
        Quadrature::VEM_Quadrature_2D quadrature;

        VEM_MCC_2D_Pressure_ReferenceElement_Data result;

        result.Monomials = monomials.Compute(order);
        result.Quadrature = quadrature.Compute_MCC_2D(order);

        result.Dimension = 2;
        result.Order = order;
        result.NumDofs0D = 0;
        result.NumDofs1D = 0;
        result.NumDofs2D = (order + 1) * (order + 2) / 2;

        return result;
    }
};

struct VEM_MCC_2D_Velocity_ReferenceElement_Data final
{
    unsigned int Dimension;
    unsigned int Order;
    unsigned int NumDofs0D; ///< Number of dofs for each vertex.
    unsigned int NumDofs1D; ///< Number of dofs internal to each edge.
    unsigned int NumDofs2D; ///< Number of dofs internal to each polygon.

    Monomials::VEM_Monomials_Data MonomialsKp1;
    Quadrature::VEM_QuadratureData_2D Quadrature;
};

class VEM_MCC_2D_Velocity_ReferenceElement final
{
  public:
    VEM_MCC_2D_Velocity_ReferenceElement_Data Create(const unsigned int order) const
    {
        Monomials::VEM_Monomials_2D monomials;
        Quadrature::VEM_Quadrature_2D quadrature;

        VEM_MCC_2D_Velocity_ReferenceElement_Data result;

        result.MonomialsKp1 = monomials.Compute(order + 1);
        result.Quadrature = quadrature.Compute_MCC_2D(order);

        result.Dimension = 2;
        result.Order = order;
        result.NumDofs0D = 0;
        result.NumDofs1D = order + 1;
        result.NumDofs2D = order * (order + 2);

        return result;
    }
};
} // namespace MCC
} // namespace VEM
} // namespace Polydim

#endif
