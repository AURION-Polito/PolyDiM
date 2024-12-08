#ifndef __VEM_MCC_2D_Ortho_LocalSpace_HPP
#define __VEM_MCC_2D_Ortho_LocalSpace_HPP

#include "Eigen/Eigen"
#include "VEM_Monomials_2D.hpp"
#include "VEM_MCC_VelocityLocalSpace_Data.hpp"
#include "VEM_MCC_2D_ReferenceElement.hpp"
#include "VEM_MCC_Utilities.hpp"
#include <vector>

namespace Polydim
{
namespace VEM
{
namespace MCC
{
class VEM_MCC_2D_Ortho_LocalSpace final
{
private:
    VEM_MCC_Utilities<2> utilities;
    Monomials::VEM_Monomials_2D monomials;

    inline void ComputeStabilizationMatrix(const double& polygonMeasure,
                                           VEM_MCC_VelocityLocalSpace_Data& localSpace) const
    {
        localSpace.StabMatrix = utilities.ComputeStabilizationMatrix(localSpace.Pi0k,
                                                                     polygonMeasure,
                                                                     localSpace.Dmatrix);
    }

public:

    /// \brief Compute matrix D: D_{ij} = dof_i(m_j).
    void ComputePolynomialBasisDofs(const double& polytopeMeasure,
                                    VEM_MCC_VelocityLocalSpace_Data& localSpace) const
    {
        localSpace.Dmatrix = utilities.ComputePolynomialBasisDofs(polytopeMeasure,
                                                                  localSpace.Order,
                                                                  localSpace.Nk,
                                                                  localSpace.NumBoundaryBasisFunctions,
                                                                  localSpace.NumNablaInternalBasisFunctions,
                                                                  localSpace.NumBigOPlusInternalBasisFunctions,
                                                                  localSpace.NumBasisFunctions,
                                                                  localSpace.GkVanderBoundaryTimesNormal,
                                                                  localSpace.Gmatrix);
    };

    VEM_MCC_VelocityLocalSpace_Data CreateLocalSpace(const VEM_MCC_2D_ReferenceElement_Data &reference_element_data,
                                                     const VEM_MCC_2D_Polygon_Geometry &polygon) const;
};
}
}
}

#endif
