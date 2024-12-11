#ifndef __VEM_DF_PCC_2D_ReferenceElement_H
#define __VEM_DF_PCC_2D_ReferenceElement_H

#include "VEM_GBasis_2D.hpp"
#include "VEM_Quadrature_2D.hpp"

namespace Polydim
{
namespace VEM
{
namespace DF_PCC
{
struct VEM_DF_PCC_2D_Pressure_ReferenceElement_Data final
{
    unsigned int Dimension;
    unsigned int Order;
    unsigned int NumDofs0D; ///< Number of dofs for each vertex.
    unsigned int NumDofs1D; ///< Number of dofs internal to each edge.
    unsigned int NumDofs2D; ///< Number of dofs internal to each polygon.

    Monomials::VEM_Monomials_Data Monomials;
    Quadrature::VEM_QuadratureData_2D Quadrature;
};

class VEM_DF_PCC_2D_Pressure_ReferenceElement final
{
public:
    VEM_DF_PCC_2D_Pressure_ReferenceElement_Data Create(const unsigned int order) const
    {
        Monomials::VEM_Monomials_2D monomials;
        Quadrature::VEM_Quadrature_2D quadrature;

        VEM_DF_PCC_2D_Pressure_ReferenceElement_Data result;

        result.Monomials = monomials.Compute(order - 1);
        result.Quadrature = quadrature.Compute_DF_PCC_2D(order);

        result.Dimension = 2;
        result.Order = order;
        result.NumDofs0D = 0;
        result.NumDofs1D = 0;
        result.NumDofs2D = order * (order + 1) / 2;

        return result;
    }
};

struct VEM_DF_PCC_2D_Velocity_ReferenceElement_Data final
{
    unsigned int Dimension;
    unsigned int Order;
    unsigned int NumDofs0D;
    unsigned int NumDofs1D;
    unsigned int NumDofs2D_BigOPlus;
    unsigned int NumDofs2D_Divergence;

    Monomials::VEM_Monomials_Data Monomials;
    Monomials::VEM_GBasis_Data GBasis;
    Quadrature::VEM_QuadratureData_2D Quadrature;
};

class VEM_DF_PCC_2D_Velocity_ReferenceElement final
{
public:
    VEM_DF_PCC_2D_Velocity_ReferenceElement_Data Create(const unsigned int order) const
    {
        Monomials::VEM_GBasis_2D g_basis;
        Monomials::VEM_Monomials_2D monomials;
        Quadrature::VEM_Quadrature_2D quadrature;

        VEM_DF_PCC_2D_Velocity_ReferenceElement_Data result;

        result.Monomials = monomials.Compute(order);
        result.GBasis = g_basis.Compute(order);
        result.Quadrature = quadrature.Compute_DF_PCC_2D(order);

        result.Dimension = 2;
        result.Order = order;
        result.NumDofs0D = 1;
        result.NumDofs1D = order - 1;
        result.NumDofs2D_Divergence = order * (order + 1) / 2 - 1;
        result.NumDofs2D_BigOPlus = (order - 2) * (order - 1) / 2;

        return result;
    }
};
}
}
}


#endif
