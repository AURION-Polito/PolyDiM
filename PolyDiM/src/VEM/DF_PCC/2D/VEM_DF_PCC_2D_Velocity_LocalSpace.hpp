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

class VEM_DF_PCC_2D_Velocity_LocalSpace final : public Polydim::VEM::DF_PCC::I_VEM_DF_PCC_2D_Velocity_LocalSpace
{
  private:
    Polydim::VEM::DF_PCC::VEM_DF_PCC_Utilities utilities;
    Polydim::Utilities::Monomials_2D monomials;
    Polydim::Utilities::GBasis_2D g_basis;

    void InitializeProjectorsComputation(const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
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
                                         Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace) const;

    void ComputeDivergenceCoefficients(const std::vector<Eigen::VectorXd> &boundaryDofQuadratureWeightsTimesNormal,
                                       Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace) const;

    void ComputePolynomialBasisDofs(const Eigen::VectorXd &internalQuadratureWeights,
                                    Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace) const;

    void ComputeCMatrixkm2(const double &polygonDiameter,
                           const std::vector<Eigen::VectorXd> &boundaryQuadratureWeightsTimesNormal,
                           Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace) const;

    void ComputePiNabla(const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                        const double &polygonMeasure,
                        const double &polygonDiameter,
                        const Eigen::VectorXd &internalQuadratureWeights,
                        const std::vector<Eigen::VectorXd> &boundaryDofQuadratureWeightsTimesNormal,
                        Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace) const;

    void ComputeL2Projectors(const Eigen::VectorXd &internalQuadratureWeights,
                             Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace) const;

    void ComputeL2ProjectorsOfDerivatives(const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                          const double &polygonDiameter,
                                          const std::vector<Eigen::VectorXd> &boundaryDofQuadratureWeightsTimesNormal,
                                          Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace) const;

  public:
    virtual ~VEM_DF_PCC_2D_Velocity_LocalSpace()
    {
    }

    Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_LocalSpace_Data CreateLocalSpace(
        const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
        const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Polygon_Geometry &polygon) const;

    inline Eigen::MatrixXd ComputeDofiDofiStabilizationMatrix(const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace,
                                                              const ProjectionTypes &projectionType) const
    {
        switch (projectionType)
        {
        case ProjectionTypes::PiNabla:
            return utilities.ComputeDofiDofiStabilizationMatrix(localSpace.Dimension,
                                                                localSpace.PiNabla,
                                                                1.0,
                                                                localSpace.Dmatrix);
        case ProjectionTypes::Pi0k:
            return utilities.ComputeDofiDofiStabilizationMatrix(localSpace.Dimension,
                                                                localSpace.Pi0k,
                                                                localSpace.Measure,
                                                                localSpace.Dmatrix);
        default:
            throw std::runtime_error("not valid projection type");
        }
    }

    inline Eigen::MatrixXd ComputeValuesOnEdge(const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                               const Eigen::VectorXd &pointsCurvilinearCoordinates) const
    {
        Eigen::RowVectorXd edgeInternalPoints;
        if (reference_element_data.Quadrature.ReferenceEdgeDOFsInternalPoints.rows() > 0)
            edgeInternalPoints = reference_element_data.Quadrature.ReferenceEdgeDOFsInternalPoints.row(0);
        const Eigen::VectorXd edgeBasisCoefficients =
            utilities.ComputeEdgeBasisCoefficients(reference_element_data.Order, edgeInternalPoints);

        return utilities.ComputeValuesOnEdge(edgeInternalPoints, reference_element_data.Order, edgeBasisCoefficients, pointsCurvilinearCoordinates);
    }

    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsValues(const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace,
                                                                    const Polydim::VEM::DF_PCC::ProjectionTypes &projectionType) const
    {
        return utilities.ComputeBasisFunctionsValues(localSpace.Dimension,
                                                     projectionType,
                                                     localSpace.Nkm2,
                                                     localSpace.Pi0km2,
                                                     localSpace.Pi0k,
                                                     localSpace.VanderInternal);
    }

    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(
        const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace,
        const Polydim::VEM::DF_PCC::ProjectionTypes &projectionType) const
    {
        return utilities.ComputeBasisFunctionsDerivativeValues(localSpace.Dimension,
                                                               projectionType,
                                                               localSpace.Nkm1,
                                                               localSpace.VanderInternal,
                                                               localSpace.VanderInternalDerivatives,
                                                               localSpace.PiNabla,
                                                               localSpace.Pi0km1Der);
    }

    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsValues(
        const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
        const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Polygon_Geometry &polygon,
        const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace,
        const Polydim::VEM::DF_PCC::ProjectionTypes &projectionType,
        const Eigen::MatrixXd &points) const
    {
        return utilities.ComputeBasisFunctionsValues(localSpace.Dimension,
                                                     projectionType,
                                                     localSpace.Nkm2,
                                                     localSpace.Pi0km2,
                                                     localSpace.Pi0k,
                                                     ComputePolynomialsValues(reference_element_data, polygon, points));
    }

    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(
        const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
        const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Polygon_Geometry &polygon,
        const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace,
        const Polydim::VEM::DF_PCC::ProjectionTypes &projectionType,
        const Eigen::MatrixXd &points) const
    {
        return utilities.ComputeBasisFunctionsDerivativeValues(
            localSpace.Dimension,
            projectionType,
            localSpace.Nkm1,
            ComputePolynomialsValues(reference_element_data, polygon, points),
            ComputePolynomialsDerivativeValues(reference_element_data, polygon, points),
            localSpace.PiNabla,
            localSpace.Pi0km1Der);
    }

    inline Eigen::MatrixXd ComputePolynomialsValues(const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace) const
    {
        return utilities.ComputePolynomialsValues(localSpace.VanderInternal);
    }

    inline Eigen::MatrixXd ComputePolynomialsValues(const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                                    const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Polygon_Geometry &polygon,
                                                    const Eigen::MatrixXd &points) const
    {
        return utilities.ComputePolynomialsValues(reference_element_data.Monomials, monomials, polygon.Centroid, polygon.Diameter, points);
    }

    inline std::vector<Eigen::MatrixXd> ComputePolynomialsDerivativeValues(const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace) const
    {
        return utilities.ComputePolynomialsDerivativeValues(localSpace.VanderInternalDerivatives);
    }

    inline std::vector<Eigen::MatrixXd> ComputePolynomialsDerivativeValues(
        const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
        const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Polygon_Geometry &polygon,
        const Eigen::MatrixXd &points) const
    {
        return utilities.ComputePolynomialsDerivativeValues(reference_element_data.Monomials,
                                                            monomials,
                                                            polygon.Diameter,
                                                            ComputePolynomialsValues(reference_element_data, polygon, points));
    }

    inline Eigen::MatrixXd ComputePolynomialsLaplacianValues(const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                                             const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Polygon_Geometry &polygon,
                                                             const Eigen::MatrixXd &points) const
    {
        return utilities.ComputePolynomialsLaplacianValues(reference_element_data.Monomials,
                                                           monomials,
                                                           polygon.Diameter,
                                                           ComputePolynomialsValues(reference_element_data, polygon, points));
    }

    inline Eigen::MatrixXd ComputeBasisFunctionsDivergenceValues(const Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace) const
    {
        return utilities.ComputeBasisFunctionsDivergenceValues(localSpace.Nkm1, localSpace.VanderInternal, localSpace.Vmatrix);
    }
};
} // namespace DF_PCC
} // namespace VEM
} // namespace Polydim

#endif
