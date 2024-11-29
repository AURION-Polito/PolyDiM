#ifndef __VEM_PCC_2D_ReferenceElement_H
#define __VEM_PCC_2D_ReferenceElement_H

#include "VEM_Monomials_2D.hpp"
#include "VEM_Quadrature_2D.hpp"

namespace Polydim
{
  namespace VEM
  {
    namespace PCC
    {
      struct VEM_PCC_2D_ReferenceElement_Data final
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
      class VEM_PCC_2D_ReferenceElement final
      {
        public:
          VEM_PCC_2D_ReferenceElement_Data Create(const unsigned int order) const
          {
            Monomials::VEM_Monomials_2D monomials;
            Quadrature::VEM_Quadrature_2D quadrature;

            VEM_PCC_2D_ReferenceElement_Data result;

            result.Monomials = monomials.Compute(order);
            result.Quadrature = quadrature.Compute(order);

            result.Dimension = 2;
            result.Order = order;
            result.NumDofs0D = 1;
            result.NumDofs1D = order - 1;
            result.NumDofs2D = order * (order - 1) / 2;

            return result;
          }
      };
    }
  }
}


#endif
