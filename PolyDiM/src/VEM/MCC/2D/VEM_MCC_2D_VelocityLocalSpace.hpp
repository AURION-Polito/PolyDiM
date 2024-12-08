#ifndef __VEM_MCC_2D_VelocityLocalSpace_HPP
#define __VEM_MCC_2D_VelocityLocalSpace_HPP

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
/// \brief Class used for computing values of basis functions of 2D
/// Mixed Conforming Constant degree Virtual Element Methods.
class VEM_MCC_2D_VelocityLocalSpace final
{
private:
    MCC::VEM_MCC_Utilities<2> utilities;
    Monomials::VEM_Monomials_2D monomials;

    void InitializeProjectorsComputation(const VEM_MCC_2D_ReferenceElement_Data& reference_element_data,
                                         const Eigen::MatrixXd& polygonVertices,
                                         const Eigen::Vector3d& polygonCentroid,
                                         const double& polygonDiameter,
                                         const Eigen::MatrixXd& internalQuadraturePoints,
                                         const Eigen::VectorXd& internalQuadratureWeights,
                                         const Eigen::MatrixXd& boundaryQuadraturePoints,
                                         VEM_MCC_VelocityLocalSpace_Data& localSpace) const;



    inline void ComputeStabilizationMatrix(const double& polygonMeasure,
                                           VEM_MCC_VelocityLocalSpace_Data& localSpace) const
    {
        localSpace.StabMatrix = utilities.ComputeStabilizationMatrix(localSpace.Pi0k,
                                                                     polygonMeasure,
                                                                     localSpace.Dmatrix);
    }

    void ComputeL2Projectors(const double &polygonMeasure,
                             const Eigen::VectorXd &internalQuadratureWeights,
                             const Eigen::MatrixXd &B2Nabla,
                             VEM_MCC_VelocityLocalSpace_Data &localSpace) const;

    void ComputeDivergenceCoefficients(const double &polytopeMeasure,
                                       const Eigen::MatrixXd &W2,
                                       VEM_MCC_VelocityLocalSpace_Data &localSpace) const;

    void ComputeValuesOnBoundary(const Eigen::MatrixXd &polytopeVertices,
                                 const Eigen::MatrixXd &edgeNormals,
                                 const std::vector<bool> &edgeDirections,
                                 const Eigen::VectorXd &boundaryQuadratureWeights,
                                 Eigen::MatrixXd &W2,
                                 Eigen::MatrixXd &B2Nabla,
                                 VEM_MCC_VelocityLocalSpace_Data &localSpace) const;
public:
    VEM_MCC_VelocityLocalSpace_Data CreateLocalSpace(const VEM_MCC_2D_ReferenceElement_Data& reference_element_data,
                                                     const VEM_MCC_2D_Polygon_Geometry& polygon) const;

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



};
}
}
}

#endif
