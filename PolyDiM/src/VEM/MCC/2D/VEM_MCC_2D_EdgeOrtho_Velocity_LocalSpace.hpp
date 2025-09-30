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

#ifndef __VEM_MCC_2D_EdgeOrtho_Velocity_LocalSpace_HPP
#define __VEM_MCC_2D_EdgeOrtho_Velocity_LocalSpace_HPP

#include "Eigen/Eigen"
#include "I_VEM_MCC_2D_ReferenceElement.hpp"
#include "I_VEM_MCC_2D_Velocity_LocalSpace.hpp"
#include "Monomials_2D.hpp"
#include "VEM_MCC_2D_LocalSpace_Data.hpp"
#include "VEM_MCC_Utilities.hpp"
#include <vector>

namespace Polydim
{
namespace VEM
{
namespace MCC
{
class VEM_MCC_2D_EdgeOrtho_Velocity_LocalSpace final : public I_VEM_MCC_2D_Velocity_LocalSpace
{
  private:
    MCC::VEM_MCC_Utilities<2> utilities;
    Utilities::Monomials_2D monomials;

    void InitializeProjectorsComputation(const VEM_MCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                         const unsigned int &numEdges,
                                         const Eigen::Vector3d &polygonCentroid,
                                         const double &polygonMeasure,
                                         const double &polygonDiameter,
                                         const Eigen::MatrixXd &internalQuadraturePoints,
                                         const Eigen::VectorXd &internalQuadratureWeights,
                                         const Eigen::MatrixXd &boundaryQuadraturePoints,
                                         VEM_MCC_2D_Velocity_LocalSpace_Data &localSpace) const;

    void ComputeL2Projectors(const double &polygonMeasure,
                             const Eigen::VectorXd &internalQuadratureWeights,
                             const Eigen::MatrixXd &B2Nabla,
                             VEM_MCC_2D_Velocity_LocalSpace_Data &localSpace) const;

    void ComputeDivergenceCoefficients(const double &polytopeMeasure,
                                       const Eigen::MatrixXd &W2,
                                       VEM_MCC_2D_Velocity_LocalSpace_Data &localSpace) const;

    void ComputeValuesOnBoundary(const VEM_MCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                 const Eigen::MatrixXd &polytopeVertices,
                                 const Eigen::VectorXd &edgeLengths,
                                 const Eigen::MatrixXd &edgeNormals,
                                 const std::vector<bool> &edgeDirections,
                                 const std::vector<Eigen::MatrixXd> &Cmatrixkp1,
                                 Eigen::MatrixXd &W2,
                                 Eigen::MatrixXd &B2Nabla,
                                 VEM_MCC_2D_Velocity_LocalSpace_Data &localSpace) const;

    void ComputePolynomialBasisDofs(const double &polytopeMeasure, VEM_MCC_2D_Velocity_LocalSpace_Data &localSpace) const
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

  public:
    VEM_MCC_2D_Velocity_LocalSpace_Data CreateLocalSpace(const VEM_MCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                                         const VEM_MCC_2D_Polygon_Geometry &polygon) const;

    inline Eigen::MatrixXd ComputeDofiDofiStabilizationMatrix(const VEM_MCC_2D_Velocity_LocalSpace_Data &localSpace,
                                                              const ProjectionTypes &projectionType) const
    {
        switch (projectionType)
        {
        case ProjectionTypes::Pi0k:
            return utilities.ComputeDofiDofiStabilizationMatrix(localSpace.Pi0k, localSpace.Measure, localSpace.Dmatrix);
        default:
            throw std::runtime_error("not valid projection type");
        }
    }

    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsValues(const VEM_MCC_2D_Velocity_LocalSpace_Data &localSpace,
                                                                    const ProjectionTypes &projectionType) const
    {
        switch (projectionType)
        {
        case ProjectionTypes::Pi0k:
            return utilities.ComputeBasisFunctionsValues(localSpace.Pi0k, localSpace.GkVanderInternal);
        default:
            throw std::runtime_error("not valid projectors");
        }
    }

    inline Eigen::MatrixXd ComputeBasisFunctionsDivergenceValues(const VEM_MCC_2D_Velocity_LocalSpace_Data &localSpace) const
    {
        return localSpace.VanderInternal * localSpace.Vmatrix;
    }

    inline Eigen::MatrixXd ComputePolynomialsValues(const VEM_MCC_2D_Velocity_LocalSpace_Data &localSpace) const
    {
        return localSpace.VanderInternal;
    }

    inline Eigen::MatrixXd ComputePolynomialsValues(const VEM_MCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                                    const VEM_MCC_2D_Velocity_LocalSpace_Data &localSpace,
                                                    const Eigen::MatrixXd &points) const
    {
        return monomials.Vander(reference_element_data.MonomialsKp1, points, localSpace.Centroid, localSpace.Diameter)
            .leftCols(localSpace.Nk);
    }
};
} // namespace MCC
} // namespace VEM
} // namespace Polydim

#endif
