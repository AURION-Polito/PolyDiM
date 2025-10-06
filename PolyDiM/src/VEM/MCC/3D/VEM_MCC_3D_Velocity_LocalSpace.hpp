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

#ifndef __VEM_MCC_3D_Velocity_LocalSpace_HPP
#define __VEM_MCC_3D_Velocity_LocalSpace_HPP

#include "Eigen/Eigen"
#include "I_VEM_MCC_3D_ReferenceElement.hpp"
#include "I_VEM_MCC_3D_Velocity_LocalSpace.hpp"
#include "Monomials_2D.hpp"
#include "Monomials_3D.hpp"
#include "VEM_MCC_3D_LocalSpace_Data.hpp"
#include "VEM_MCC_Utilities.hpp"
#include <vector>

namespace Polydim
{
namespace VEM
{
namespace MCC
{
class VEM_MCC_3D_Velocity_LocalSpace final : public I_VEM_MCC_3D_Velocity_LocalSpace
{
  private:
    Polydim::VEM::MCC::VEM_MCC_Utilities utilities;
    Polydim::Utilities::Monomials_3D monomials3D;
    Polydim::Utilities::Monomials_2D monomials2D;

    void InitializeProjectorsComputation(const Polydim::VEM::MCC::VEM_MCC_3D_Velocity_ReferenceElement_Data &reference_element_data,
                                         const unsigned int &numFaces,
                                         const Eigen::Vector3d &polyhedronCentroid,
                                         const double &polyhedronMeasure,
                                         const double &polyhedronDiameter,
                                         const Eigen::MatrixXd &internalQuadraturePoints,
                                         const Eigen::VectorXd &internalQuadratureWeights,
                                         const Eigen::MatrixXd &boundaryQuadraturePoints,
                                         Polydim::VEM::MCC::VEM_MCC_3D_Velocity_LocalSpace_Data &localSpace) const;

    void ComputeL2Projectors(const double &polyhedronMeasure,
                             const Eigen::VectorXd &internalQuadratureWeights,
                             const Eigen::MatrixXd &B2Nabla,
                             Polydim::VEM::MCC::VEM_MCC_3D_Velocity_LocalSpace_Data &localSpace) const;

    void ComputeDivergenceCoefficients(const double &polytopeMeasure,
                                       const Eigen::MatrixXd &W2,
                                       Polydim::VEM::MCC::VEM_MCC_3D_Velocity_LocalSpace_Data &localSpace) const;

    void ComputeValuesOnBoundary(const Polydim::VEM::MCC::VEM_MCC_3D_Velocity_ReferenceElement_Data &reference_element_data,
                                 const unsigned int &numFaces,
                                 const std::vector<Eigen::Vector3d> &facesNormals,
                                 const std::vector<bool> &facesNormalDirections,
                                 const std::vector<bool> &faceNormalGlobalDirections,
                                 const std::vector<Eigen::Vector3d> &facesCentroids,
                                 const std::vector<double> &facesAreas,
                                 const std::vector<double> &facesDiameters,
                                 const std::vector<Gedim::Quadrature::QuadratureData> &facesQuadrature,
                                 const Eigen::VectorXd &boundaryQuadratureWeights,
                                 Eigen::MatrixXd &W2,
                                 Eigen::MatrixXd &B2Nabla,
                                 Polydim::VEM::MCC::VEM_MCC_3D_Velocity_LocalSpace_Data &localSpace) const;

    void ComputePolynomialBasisDofs(const double &polytopeMeasure, Polydim::VEM::MCC::VEM_MCC_3D_Velocity_LocalSpace_Data &localSpace) const
    {
        localSpace.Dmatrix = utilities.ComputePolynomialBasisDofs(localSpace.Dimension,
                                                                  polytopeMeasure,
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
    Polydim::VEM::MCC::VEM_MCC_3D_Velocity_LocalSpace_Data CreateLocalSpace(
        const Polydim::VEM::MCC::VEM_MCC_3D_Velocity_ReferenceElement_Data &reference_element_data,
        const Polydim::VEM::MCC::VEM_MCC_3D_Polyhedron_Geometry &polyhedron) const;

    inline Eigen::MatrixXd ComputeDofiDofiStabilizationMatrix(const Polydim::VEM::MCC::VEM_MCC_3D_Velocity_LocalSpace_Data &localSpace,
                                                              const Polydim::VEM::MCC::ProjectionTypes &projectionType) const
    {
        switch (projectionType)
        {
        case Polydim::VEM::MCC::ProjectionTypes::Pi0k:
            return utilities.ComputeDofiDofiStabilizationMatrix(localSpace.Pi0k, localSpace.Measure, localSpace.Dmatrix);
        default:
            throw std::runtime_error("not valid projection type");
        }
    }

    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsValues(const Polydim::VEM::MCC::VEM_MCC_3D_Velocity_LocalSpace_Data &localSpace,
                                                                    const Polydim::VEM::MCC::ProjectionTypes &projectionType) const
    {
        switch (projectionType)
        {
        case Polydim::VEM::MCC::ProjectionTypes::Pi0k:
            return utilities.ComputeBasisFunctionsValues(localSpace.Dimension, localSpace.Pi0k, localSpace.GkVanderInternal);
        default:
            throw std::runtime_error("not valid projectors");
        }
    }

    inline Eigen::MatrixXd ComputeBasisFunctionsDivergenceValues(const Polydim::VEM::MCC::VEM_MCC_3D_Velocity_LocalSpace_Data &localSpace) const
    {
        return localSpace.VanderInternal * localSpace.Vmatrix;
    }

    inline Eigen::MatrixXd ComputePolynomialsValues(const Polydim::VEM::MCC::VEM_MCC_3D_Velocity_LocalSpace_Data &localSpace) const
    {
        return localSpace.VanderInternal;
    }

    inline Eigen::MatrixXd ComputePolynomialsValues(const Polydim::VEM::MCC::VEM_MCC_3D_Velocity_ReferenceElement_Data &reference_element_data,
                                                    const Polydim::VEM::MCC::VEM_MCC_3D_Velocity_LocalSpace_Data &localSpace,
                                                    const Eigen::MatrixXd &points) const
    {
        return monomials3D.Vander(reference_element_data.MonomialsKp1, points, localSpace.Centroid, localSpace.Diameter)
            .leftCols(localSpace.Nk);
    }
};
} // namespace MCC
} // namespace VEM
} // namespace Polydim

#endif
