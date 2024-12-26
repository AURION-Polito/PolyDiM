#ifndef __I_VEM_MCC_2D_ReferenceElement_H
#define __I_VEM_MCC_2D_ReferenceElement_H

#include "VEM_Monomials_Data.hpp"
#include "VEM_Quadrature_2D.hpp"

namespace Polydim
{
namespace VEM
{
namespace MCC
{

struct EdgeOrtho final
{
    Eigen::MatrixXd Hmatrix1D;
    Eigen::MatrixXd QmatrixInvKp1_1D;
    Eigen::MatrixXd QmatrixKp1_1D;
};

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

struct VEM_MCC_2D_Velocity_ReferenceElement_Data final
{
    unsigned int Dimension;
    unsigned int Order;
    unsigned int NumDofs0D; ///< Number of dofs for each vertex.
    unsigned int NumDofs1D; ///< Number of dofs internal to each edge.
    unsigned int NumDofs2D; ///< Number of dofs internal to each polygon.

    EdgeOrtho edge_ortho;
    Monomials::VEM_Monomials_Data Monomials_1D;
    Monomials::VEM_Monomials_Data MonomialsKp1;
    Quadrature::VEM_QuadratureData_2D Quadrature;
};

class I_VEM_MCC_2D_Velocity_ReferenceElement
{
  public:
    virtual VEM_MCC_2D_Velocity_ReferenceElement_Data Create(const unsigned int order) const = 0;
};

class I_VEM_MCC_2D_Pressure_ReferenceElement
{
  public:
    virtual VEM_MCC_2D_Pressure_ReferenceElement_Data Create(const unsigned int order) const = 0;
};
} // namespace MCC
} // namespace VEM
} // namespace Polydim

#endif
