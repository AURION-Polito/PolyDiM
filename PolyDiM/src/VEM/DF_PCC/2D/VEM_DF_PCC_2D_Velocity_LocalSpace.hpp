// _LICENSE_HEADER_
//
// Copyright (C) 2019 - 2025.
// Terms register on the GPL-3.0 license.
//
// This file can be redistributed and/or modified under the license terms.
//
// See top level LICENSE file for more details.
//
// This file can be used citing references in CITATION.cff file.

#ifndef __VEM_DF_PCC_2D_Velocity_LocalSpace_HPP
#define __VEM_DF_PCC_2D_Velocity_LocalSpace_HPP

#include "Eigen/Eigen"
#include "GBasis_2D.hpp"
#include "I_VEM_DF_PCC_2D_ReferenceElement.hpp"
#include "I_VEM_DF_PCC_2D_Velocity_LocalSpace.hpp"
#include "Monomials_2D.hpp"
#include "VEM_DF_PCC_2D_LocalSpace_Data.hpp"
#include "VEM_DF_PCC_Utilities.hpp"
#include <vector>

namespace Polydim
{
namespace VEM
{
namespace DF_PCC
{

class VEM_DF_PCC_2D_Velocity_LocalSpace final : public I_VEM_DF_PCC_2D_Velocity_LocalSpace
{
  private:
    VEM_DF_PCC_Utilities<2> utilities;
    Utilities::Monomials_2D monomials;
    Utilities::GBasis_2D g_basis;

    void InitializeProjectorsComputation(const VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                         const Eigen::MatrixXd &polygonVertices,
                                         const Eigen::Vector3d &polygonCentroid,
                                         const double &polygonMeasure,
                                         const double &polygonDiameter,
                                         const std::vector<bool> &edgeDirections,
                                         const Eigen::MatrixXd &internalQuadraturePoints,
                                         const Eigen::VectorXd &internalQuadratureWeights,
                                         const Eigen::MatrixXd &boundaryQuadraturePoints,
                                         const Eigen::MatrixXd &boundaryDofQuadraturePoints,
                                         const Eigen::MatrixXd &referenceEdgeInternalPoints,
                                         const Eigen::MatrixXd &referenceEdgeDofInternalPoints,
                                         VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace) const;

    void ComputeDivergenceCoefficients(const std::vector<Eigen::VectorXd> &boundaryDofQuadratureWeightsTimesNormal,
                                       VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace) const;

    void ComputePolynomialBasisDofs(const Eigen::VectorXd &internalQuadratureWeights,
                                    VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace) const;

    void ComputeCMatrixkm2(const double &polygonDiameter,
                           const std::vector<Eigen::VectorXd> &boundaryQuadratureWeightsTimesNormal,
                           VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace) const;

    void ComputePiNabla(const VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                        const double &polygonMeasure,
                        const double &polygonDiameter,
                        const Eigen::VectorXd &internalQuadratureWeights,
                        const std::vector<Eigen::VectorXd> &boundaryDofQuadratureWeightsTimesNormal,
                        VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace) const;

    void ComputeL2Projectors(const Eigen::VectorXd &internalQuadratureWeights, VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace) const;

    void ComputeL2ProjectorsOfDerivatives(const VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                          const double &polygonDiameter,
                                          const std::vector<Eigen::VectorXd> &boundaryDofQuadratureWeightsTimesNormal,
                                          VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace) const;

  public:
    virtual ~VEM_DF_PCC_2D_Velocity_LocalSpace()
    {
    }

    VEM_DF_PCC_2D_Velocity_LocalSpace_Data CreateLocalSpace(const VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                                            const VEM_DF_PCC_2D_Polygon_Geometry &polygon) const;

    inline Eigen::MatrixXd ComputeDofiDofiStabilizationMatrix(const VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace,
                                                              const ProjectionTypes &projectionType) const
    {
        switch (projectionType)
        {
        case ProjectionTypes::PiNabla:
            return utilities.ComputeDofiDofiStabilizationMatrix(localSpace.PiNabla, 1.0, localSpace.Dmatrix);
        case ProjectionTypes::Pi0k:
            return utilities.ComputeDofiDofiStabilizationMatrix(localSpace.Pi0k, localSpace.Measure, localSpace.Dmatrix);
        default:
            throw std::runtime_error("not valid projection type");
        }
    }

    inline Eigen::MatrixXd ComputeValuesOnEdge(const VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                               const Eigen::VectorXd &pointsCurvilinearCoordinates) const
    {
        Eigen::RowVectorXd edgeInternalPoints;
        if (reference_element_data.Quadrature.ReferenceEdgeDOFsInternalPoints.rows() > 0)
            edgeInternalPoints = reference_element_data.Quadrature.ReferenceEdgeDOFsInternalPoints.row(0);
        const Eigen::VectorXd edgeBasisCoefficients =
            utilities.ComputeEdgeBasisCoefficients(reference_element_data.Order, edgeInternalPoints);

        return utilities.ComputeValuesOnEdge(edgeInternalPoints, reference_element_data.Order, edgeBasisCoefficients, pointsCurvilinearCoordinates);
    }

    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsValues(const VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace,
                                                                    const ProjectionTypes &projectionType) const
    {
        return utilities.ComputeBasisFunctionsValues(projectionType,
                                                     localSpace.Nkm2,
                                                     localSpace.Pi0km2,
                                                     localSpace.Pi0k,
                                                     localSpace.VanderInternal);
    }

    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace,
                                                                              const ProjectionTypes &projectionType) const
    {
        return utilities.ComputeBasisFunctionsDerivativeValues(projectionType,
                                                               localSpace.Nkm1,
                                                               localSpace.VanderInternal,
                                                               localSpace.VanderInternalDerivatives,
                                                               localSpace.PiNabla,
                                                               localSpace.Pi0km1Der);
    }

    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsValues(const VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                                                    const VEM_DF_PCC_2D_Polygon_Geometry &polygon,
                                                                    const VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace,
                                                                    const ProjectionTypes &projectionType,
                                                                    const Eigen::MatrixXd &points) const
    {
        return utilities.ComputeBasisFunctionsValues(projectionType,
                                                     localSpace.Nkm2,
                                                     localSpace.Pi0km2,
                                                     localSpace.Pi0k,
                                                     ComputePolynomialsValues(reference_element_data, polygon, points));
    }

    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                                                              const VEM_DF_PCC_2D_Polygon_Geometry &polygon,
                                                                              const VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace,
                                                                              const ProjectionTypes &projectionType,
                                                                              const Eigen::MatrixXd &points) const
    {
        return utilities.ComputeBasisFunctionsDerivativeValues(
            projectionType,
            localSpace.Nkm1,
            ComputePolynomialsValues(reference_element_data, polygon, points),
            ComputePolynomialsDerivativeValues(reference_element_data, polygon, points),
            localSpace.PiNabla,
            localSpace.Pi0km1Der);
    }

    inline Eigen::MatrixXd ComputePolynomialsValues(const VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace) const
    {
        return utilities.ComputePolynomialsValues(localSpace.VanderInternal);
    }

    inline Eigen::MatrixXd ComputePolynomialsValues(const VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                                    const VEM_DF_PCC_2D_Polygon_Geometry &polygon,
                                                    const Eigen::MatrixXd &points) const
    {
        return utilities.ComputePolynomialsValues(reference_element_data.Monomials, monomials, polygon.Centroid, polygon.Diameter, points);
    }

    inline std::vector<Eigen::MatrixXd> ComputePolynomialsDerivativeValues(const VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace) const
    {
        return utilities.ComputePolynomialsDerivativeValues(localSpace.VanderInternalDerivatives);
    }

    inline std::vector<Eigen::MatrixXd> ComputePolynomialsDerivativeValues(const VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                                                           const VEM_DF_PCC_2D_Polygon_Geometry &polygon,
                                                                           const Eigen::MatrixXd &points) const
    {
        return utilities.ComputePolynomialsDerivativeValues(reference_element_data.Monomials,
                                                            monomials,
                                                            polygon.Diameter,
                                                            ComputePolynomialsValues(reference_element_data, polygon, points));
    }

    inline Eigen::MatrixXd ComputePolynomialsLaplacianValues(const VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                                             const VEM_DF_PCC_2D_Polygon_Geometry &polygon,
                                                             const Eigen::MatrixXd &points) const
    {
        return utilities.ComputePolynomialsLaplacianValues(reference_element_data.Monomials,
                                                           monomials,
                                                           polygon.Diameter,
                                                           ComputePolynomialsValues(reference_element_data, polygon, points));
    }

    inline Eigen::MatrixXd ComputeBasisFunctionsDivergenceValues(const VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace) const
    {
        return utilities.ComputeBasisFunctionsDivergenceValues(localSpace.Nkm1, localSpace.VanderInternal, localSpace.Vmatrix);
    }
};
} // namespace DF_PCC
} // namespace VEM
} // namespace Polydim

#endif
