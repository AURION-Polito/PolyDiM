#ifndef __VEM_DF_PCC_2D_LocalSpace_HPP
#define __VEM_DF_PCC_2D_LocalSpace_HPP

#include "Eigen/Eigen"
#include "VEM_Monomials_2D.hpp"
#include "VEM_DF_PCC_2D_LocalSpace_Data.hpp"
#include "VEM_DF_PCC_2D_ReferenceElement.hpp"
#include "VEM_DF_PCC_Utilities.hpp"
#include <vector>

namespace Polydim {
namespace VEM {
namespace DF_PCC {

/// \brief Class used for computing values of basis functions of 2D
/// Divergence Free Primal Conforming Constant degree Virtual Element Methods.
///
/// Please cite the following article:
///     - <a href="https://doi.org/10.1016/j.matcom.2023.10.003">"Improving high-order VEM stability on badly-shaped elements. Stefano Berrone, Gioana Teora and Fabio Vicini. (2024)"</a>

class VEM_DF_PCC_2D_LocalSpace final
{
private:
    VEM_DF_PCC_Utilities<2> utilities;
    Monomials::VEM_Monomials_2D monomials;
    Monomials::VEM_GBasis_2D g_basis;

    void InitializeProjectorsComputation(const VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                         const Eigen::MatrixXd &polygonVertices,
                                         const Eigen::Vector3d &polygonCentroid,
                                         const double &polygonDiameter,
                                         const std::vector<bool> &edgeDirections,
                                         const Eigen::MatrixXd &internalQuadraturePoints,
                                         const Eigen::VectorXd &internalQuadratureWeights,
                                         const Eigen::MatrixXd &boundaryQuadraturePoints,
                                         const Eigen::MatrixXd &boundaryDofQuadraturePoints,
                                         const Eigen::MatrixXd &referenceEdgeInternalPoints,
                                         const Eigen::MatrixXd &referenceEdgeDofInternalPoints,
                                         VEM_DF_PCC_2D_LocalSpace_Data &localSpace) const;

    void ComputeDivergenceCoefficients(const std::vector<Eigen::VectorXd> &boundaryDofQuadratureWeightsTimesNormal,
                                       VEM_DF_PCC_2D_LocalSpace_Data &localSpace) const;

    void ComputeStabilizationMatrix(VEM_DF_PCC_2D_LocalSpace_Data &localSpace) const;

    void ComputePolynomialBasisDofs(const Eigen::VectorXd &internalQuadratureWeights,
                                    VEM_DF_PCC_2D_LocalSpace_Data &localSpace) const;

    void ComputeCMatrixkm2(const double &polygonDiameter,
                           const std::vector<Eigen::VectorXd> &boundaryQuadratureWeightsTimesNormal,
                           VEM_DF_PCC_2D_LocalSpace_Data &localSpace) const;

    void ComputePiNabla(const VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                        const double &polygonMeasure,
                        const double &polygonDiameter,
                        const Eigen::VectorXd &internalQuadratureWeights,
                        const std::vector<Eigen::VectorXd> &boundaryDofQuadratureWeightsTimesNormal,
                        VEM_DF_PCC_2D_LocalSpace_Data &localSpace) const;

    void ComputeL2Projectors(const Eigen::VectorXd &internalQuadratureWeights,
                             VEM_DF_PCC_2D_LocalSpace_Data &localSpace) const;

    void ComputeL2ProjectorsOfDerivatives(const VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                          const double &polygonDiameter,
                                          const std::vector<Eigen::VectorXd> &boundaryDofQuadratureWeightsTimesNormal,
                                          VEM_DF_PCC_2D_LocalSpace_Data &localSpace) const;

public:

    VEM_DF_PCC_2D_LocalSpace_Data CreateLocalSpace(const VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                                   const VEM_DF_PCC_2D_Polygon_Geometry &polygon) const;

    inline Eigen::MatrixXd ComputeValuesOnEdge(const VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                               const Eigen::VectorXd &pointsCurvilinearCoordinates) const
    {
        const Eigen::VectorXd edgeInternalPoints = reference_element_data.Quadrature.ReferenceSegmentInternalPoints;
        const Eigen::VectorXd edgeBasisCoefficients = utilities.ComputeEdgeBasisCoefficients(reference_element_data.Order, edgeInternalPoints);

        return utilities.ComputeValuesOnEdge(edgeInternalPoints.transpose(),
                                             reference_element_data.Order,
                                             edgeBasisCoefficients,
                                             pointsCurvilinearCoordinates);
    }


};
} // namespace DF_PCC
} // namespace VEM
} // namespace Polydim

#endif
