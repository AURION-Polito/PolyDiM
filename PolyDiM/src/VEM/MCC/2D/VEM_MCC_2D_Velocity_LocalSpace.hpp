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

#ifndef __VEM_MCC_2D_Velocity_LocalSpace_HPP
#define __VEM_MCC_2D_Velocity_LocalSpace_HPP

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

/// \brief Interface class for the velocity space of 2D Mixed Conforming Constant degree Virtual Element Methods \cite
/// secondMixed \cite DaVeiga2016.
class VEM_MCC_2D_Velocity_LocalSpace final : public I_VEM_MCC_2D_Velocity_LocalSpace
{
private:
    Polydim::VEM::MCC::VEM_MCC_Utilities utilities;
    Polydim::Utilities::Monomials_2D monomials;

    void InitializeProjectorsComputation(const Polydim::VEM::MCC::VEM_MCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                         const unsigned int &numEdges,
                                         const Eigen::Vector3d &polygonCentroid,
                                         const double &polygonMeasure,
                                         const double &polygonDiameter,
                                         const Eigen::MatrixXd &internalQuadraturePoints,
                                         const Eigen::VectorXd &internalQuadratureWeights,
                                         const Eigen::MatrixXd &boundaryQuadraturePoints,
                                         Polydim::VEM::MCC::VEM_MCC_2D_Velocity_LocalSpace_Data &localSpace) const;

    void ComputeL2Projectors(const double &polygonMeasure,
                             const Eigen::VectorXd &internalQuadratureWeights,
                             const Eigen::MatrixXd &B2Nabla,
                             Polydim::VEM::MCC::VEM_MCC_2D_Velocity_LocalSpace_Data &localSpace) const;

    void ComputeDivergenceCoefficients(const double &polytopeMeasure,
                                       const Eigen::MatrixXd &W2,
                                       Polydim::VEM::MCC::VEM_MCC_2D_Velocity_LocalSpace_Data &localSpace) const;

    void ComputeValuesOnBoundary(const Eigen::MatrixXd &polytopeVertices,
                                 const Eigen::MatrixXd &edgeNormals,
                                 const std::vector<bool> &edgeDirections,
                                 const Eigen::VectorXd &boundaryQuadratureWeights,
                                 Eigen::MatrixXd &W2,
                                 Eigen::MatrixXd &B2Nabla,
                                 Polydim::VEM::MCC::VEM_MCC_2D_Velocity_LocalSpace_Data &localSpace) const;

    void ComputePolynomialBasisDofs(const double &polytopeMeasure, Polydim::VEM::MCC::VEM_MCC_2D_Velocity_LocalSpace_Data &localSpace) const
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
    Polydim::VEM::MCC::VEM_MCC_2D_Velocity_LocalSpace_Data CreateLocalSpace(
        const Polydim::VEM::MCC::VEM_MCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
        const Polydim::VEM::MCC::VEM_MCC_2D_Polygon_Geometry &polygon) const;

    inline Eigen::MatrixXd ComputeDofiDofiStabilizationMatrix(const Polydim::VEM::MCC::VEM_MCC_2D_Velocity_LocalSpace_Data &localSpace,
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

    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsValues(const Polydim::VEM::MCC::VEM_MCC_2D_Velocity_LocalSpace_Data &localSpace,
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

    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsValues(const Polydim::VEM::MCC::VEM_MCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                                                    const Polydim::VEM::MCC::VEM_MCC_2D_Velocity_LocalSpace_Data &localSpace,
                                                                    const Polydim::VEM::MCC::ProjectionTypes &projectionType,
                                                                    const Eigen::MatrixXd &points) const
    {
        switch(projectionType)
        {

        case ProjectionTypes::Pi0k:
        {
            const unsigned int numQuadrature = points.cols();

            const Eigen::MatrixXd VanderInternalKp1 = monomials.Vander(reference_element_data.MonomialsKp1,
                                                                       points,
                                                                       localSpace.Centroid,
                                                                       localSpace.Diameter);
            const Eigen::MatrixXd VanderInternal = VanderInternalKp1.leftCols(localSpace.Nk);
            Eigen::MatrixXd VanderInternal2k(localSpace.Dimension * localSpace.Nk, localSpace.Dimension * points.cols());
            VanderInternal2k << VanderInternal.transpose(),
                Eigen::MatrixXd::Zero(VanderInternal.cols(), VanderInternal.rows()),
                Eigen::MatrixXd::Zero(VanderInternal.cols(), VanderInternal.rows()),
                VanderInternal.transpose();

            const Eigen::MatrixXd GkNablaVanderInternal = localSpace.TkNabla * VanderInternal2k;
            const Eigen::MatrixXd GkBigOPlusVanderInternal = localSpace.TkBigOPlus * VanderInternal2k;

            Eigen::MatrixXd GkVanderInternal = Eigen::MatrixXd::Zero(localSpace.Dimension * localSpace.Nk,
                                                                     localSpace.Dimension * numQuadrature);

            GkVanderInternal << GkNablaVanderInternal, GkBigOPlusVanderInternal;

            const Eigen::MatrixXd temp = GkVanderInternal.transpose() * localSpace.Pi0k;
            std::vector<Eigen::MatrixXd> result(localSpace.Dimension, Eigen::MatrixXd::Zero(localSpace.Dimension, localSpace.Pi0k.cols()));

            for (unsigned int d = 0; d < localSpace.Dimension; d++)
                result[d] = temp.middleRows(numQuadrature * d, numQuadrature);

            return result;
        }
        default:
            throw std::runtime_error("not valid projector type");
        }

    }

    inline Eigen::MatrixXd ComputeBasisFunctionsDivergenceValues(const Polydim::VEM::MCC::VEM_MCC_2D_Velocity_LocalSpace_Data &localSpace) const
    {
        return localSpace.VanderInternal * localSpace.Vmatrix;
    }

    inline Eigen::MatrixXd ComputePolynomialsValues(const Polydim::VEM::MCC::VEM_MCC_2D_Velocity_LocalSpace_Data &localSpace) const
    {
        return localSpace.VanderInternal;
    }

    inline Eigen::MatrixXd ComputePolynomialsValues(const Polydim::VEM::MCC::VEM_MCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                                    const Polydim::VEM::MCC::VEM_MCC_2D_Velocity_LocalSpace_Data &localSpace,
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
