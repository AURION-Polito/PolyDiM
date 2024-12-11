#ifndef __VEM_DF_PCC_3D_ReferenceElement_H
#define __VEM_DF_PCC_3D_ReferenceElement_H

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

class VEM_DF_PCC_3D_Pressure_ReferenceElement final
{
public:
    VEM_DF_PCC_3D_Pressure_ReferenceElement_Data Create(const unsigned int order) const
    {
        Monomials::VEM_Monomials_3D monomials;
        Quadrature::VEM_Quadrature_3D quadrature;

        VEM_DF_PCC_3D_Pressure_ReferenceElement_Data result;

        result.Monomials = monomials.Compute(order - 1);
        result.Quadrature = quadrature.Compute_DF_PCC_3D(order);

        result.Dimension = 3;
        result.Order = order;
        result.NumDofs0D = 0;
        result.NumDofs1D = 0;
        result.NumDofs2D = 0;
        result.NumDofs2D = order * (order + 1) * (order + 2) / 6;

        return result;
    }
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

class VEM_DF_PCC_3D_Velocity_ReferenceElement final
{
public:
    VEM_DF_PCC_3D_Velocity_ReferenceElement_Data Create(const unsigned int order) const
    {
        Monomials::VEM_GBasis_3D g_basis;
        Monomials::VEM_Monomials_3D monomials;
        Quadrature::VEM_Quadrature_3D quadrature;

        VEM_DF_PCC_3D_Velocity_ReferenceElement_Data result;

        result.Monomials = monomials.Compute(order);
        result.GBasis = g_basis.Compute(order);
        result.Quadrature = quadrature.Compute_DF_PCC_3D(order);

        result.Dimension = 3;
        result.Order = order;
        result.NumDofs0D = 1;
        result.NumDofs1D = order - 1;
        result.NumDofs2D = order * (order - 1) / 2;
        result.NumDofs3D_Divergence = order * (order + 1) * (order + 2) / 6 - 1;
        result.NumDofs3D_BigOPlus = order * (order + 1) * (order - 1) / 2 - result.NumDofs3D_Divergence;

        return result;
    }
};
}
}
}


#endif
