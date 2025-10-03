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

#ifndef __VEM_PCC_2D_Inertia_LocalSpace_HPP
#define __VEM_PCC_2D_Inertia_LocalSpace_HPP

#include "I_VEM_PCC_2D_LocalSpace.hpp"
#include "Monomials_2D.hpp"

namespace Polydim
{
namespace VEM
{
namespace PCC
{

/// \brief Primal Conforming Constant degree Virtual Element Methods 2D with improvements for high-order \cite
/// Teora2024.
class VEM_PCC_2D_Inertia_LocalSpace final : public I_VEM_PCC_2D_LocalSpace
{
  private:
    VEM_PCC_Utilities<2> utilities;
    Utilities::Monomials_2D monomials;

    void InitializeProjectorsComputation(const Polydim::VEM::PCC::VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                         const Eigen::MatrixXd &polygonVertices,
                                         const Eigen::Vector3d &polygonCentroid,
                                         const double &polygonDiameter,
                                         const Eigen::MatrixXd &internalQuadraturePoints,
                                         const Eigen::VectorXd &internalQuadratureWeights,
                                         const Eigen::MatrixXd &boundaryQuadraturePoints,
                                         Polydim::VEM::PCC::VEM_PCC_2D_LocalSpace_Data &localSpace) const;

    void ComputePiNabla(const Polydim::VEM::PCC::VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                        const double &polygonMeasure,
                        const double &polygonDiameter,
                        const Eigen::VectorXd &internalQuadratureWeights,
                        const Eigen::VectorXd &boundaryQuadratureWeights,
                        const std::vector<Eigen::VectorXd> &boundaryQuadratureWeightsTimesNormal,
                        Polydim::VEM::PCC::VEM_PCC_2D_LocalSpace_Data &localSpace) const;

    void ComputeL2ProjectorsOfDerivatives(const Polydim::VEM::PCC::VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                          const double &polygonMeasure,
                                          const double &polygonDiameter,
                                          const std::vector<Eigen::VectorXd> &boundaryQuadratureWeightsTimesNormal,
                                          Polydim::VEM::PCC::VEM_PCC_2D_LocalSpace_Data &localSpace) const;

    void ComputeL2Projectors(const double &polygonMeasure, Polydim::VEM::PCC::VEM_PCC_2D_LocalSpace_Data &localSpace) const
    {
        utilities.ComputeL2Projectors(polygonMeasure,
                                      localSpace.Order,
                                      localSpace.Nkm1,
                                      localSpace.NumProjectorBasisFunctions,
                                      localSpace.NumInternalBasisFunctions,
                                      localSpace.NumBasisFunctions,
                                      localSpace.Hmatrix,
                                      localSpace.PiNabla,
                                      localSpace.Cmatrix,
                                      localSpace.Pi0km1,
                                      localSpace.Pi0k);
    };

    void ComputeGeometryProperties(const Gedim::GeometryUtilities &geometryUtilities,
                                   const Polydim::Utilities::Inertia_Data &inertia_data,
                                   const PCC::VEM_PCC_2D_Polygon_Geometry &polygon,
                                   PCC::VEM_PCC_2D_Polygon_Geometry &inertia_geometric_data) const;

    void ComputePolynomialsDofs(const double &polytopeMeasure, Polydim::VEM::PCC::VEM_PCC_2D_LocalSpace_Data &localSpace) const;

  public:
    Polydim::VEM::PCC::VEM_PCC_2D_LocalSpace_Data CreateLocalSpace(const Polydim::VEM::PCC::VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                                   const Polydim::VEM::PCC::VEM_PCC_2D_Polygon_Geometry &polygon) const;

    Polydim::VEM::PCC::VEM_PCC_2D_LocalSpace_Data Compute3DUtilities(const Polydim::VEM::PCC::VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                                     const Polydim::VEM::PCC::VEM_PCC_2D_Polygon_Geometry &polygon) const;

    inline Eigen::MatrixXd ComputeDofiDofiStabilizationMatrix(const Polydim::VEM::PCC::VEM_PCC_2D_LocalSpace_Data &localSpace,
                                                              const Polydim::VEM::PCC::ProjectionTypes &projectionType) const
    {
        switch (projectionType)
        {
        case Polydim::VEM::PCC::ProjectionTypes::PiNabla:
            return utilities.ComputeDofiDofiStabilizationMatrix(localSpace.PiNabla, localSpace.constantStiff, localSpace.Dmatrix);
        case Polydim::VEM::PCC::ProjectionTypes::Pi0k:
            return utilities.ComputeDofiDofiStabilizationMatrix(localSpace.Pi0k, localSpace.constantMass, localSpace.Dmatrix);
        default:
            throw std::runtime_error("not valid projection type");
        }
    }

    inline Eigen::MatrixXd ComputeBasisFunctionsValues(const Polydim::VEM::PCC::VEM_PCC_2D_LocalSpace_Data &localSpace,
                                                       const Polydim::VEM::PCC::ProjectionTypes &projectionType) const
    {
        return utilities.ComputeBasisFunctionsValues(projectionType,
                                                     localSpace.Nkm1,
                                                     localSpace.Pi0km1,
                                                     localSpace.Pi0k,
                                                     localSpace.VanderInternal);
    }

    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const Polydim::VEM::PCC::VEM_PCC_2D_LocalSpace_Data &localSpace,
                                                                              const Polydim::VEM::PCC::ProjectionTypes &projectionType) const
    {
        std::vector<Eigen::MatrixXd> basisFunctionsDerivativeValues(localSpace.Dimension);
        const Eigen::MatrixXd FmatrixInvTransp = localSpace.inertia_data.FmatrixInv.transpose();
        switch (projectionType)
        {
        case Polydim::VEM::PCC::ProjectionTypes::PiNabla: {
            for (unsigned short i = 0; i < localSpace.Dimension; ++i)
                basisFunctionsDerivativeValues[i] = localSpace.VanderInternalDerivatives[i] * localSpace.PiNabla;
        }
        break;
        case Polydim::VEM::PCC::ProjectionTypes::Pi0km1Der: {
            for (unsigned short i = 0; i < localSpace.Dimension; ++i)
                basisFunctionsDerivativeValues[i] = localSpace.VanderInternal.leftCols(localSpace.Nkm1) * localSpace.Pi0km1Der[i];
        }
        break;
        default:
            throw std::runtime_error("Unknown projector type");
        }

        std::vector<Eigen::MatrixXd> fmatrixInvTranspTimesBasisFunctionDerivativeValues2D(
            localSpace.Dimension,
            Eigen::MatrixXd::Zero(basisFunctionsDerivativeValues[0].rows(), basisFunctionsDerivativeValues[1].cols()));
        for (unsigned int d1 = 0; d1 < localSpace.Dimension; d1++)
        {
            for (unsigned int d2 = 0; d2 < localSpace.Dimension; d2++)
            {
                fmatrixInvTranspTimesBasisFunctionDerivativeValues2D[d1] +=
                    FmatrixInvTransp(d1, d2) * basisFunctionsDerivativeValues[d2];
            }
        }

        return fmatrixInvTranspTimesBasisFunctionDerivativeValues2D;
    }

    inline Eigen::MatrixXd ComputeBasisFunctionsValues(const Polydim::VEM::PCC::VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                       const Polydim::VEM::PCC::VEM_PCC_2D_LocalSpace_Data &localSpace,
                                                       const Polydim::VEM::PCC::ProjectionTypes &projectionType,
                                                       const Eigen::MatrixXd &points) const
    {
        return utilities.ComputeBasisFunctionsValues(projectionType,
                                                     localSpace.Nkm1,
                                                     localSpace.Pi0km1,
                                                     localSpace.Pi0k,
                                                     ComputePolynomialsValues(reference_element_data, localSpace, points));
    }

    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(
        const Polydim::VEM::PCC::VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
        const Polydim::VEM::PCC::VEM_PCC_2D_LocalSpace_Data &localSpace,
        const Polydim::VEM::PCC::ProjectionTypes &projectionType,
        const Eigen::MatrixXd &points) const
    {
        const Eigen::MatrixXd vander = ComputePolynomialsValues(reference_element_data, localSpace, points);

        std::vector<Eigen::MatrixXd> basisFunctionsDerivativeValues(localSpace.Dimension);
        const Eigen::MatrixXd FmatrixInvTransp = localSpace.inertia_data.FmatrixInv.transpose();
        switch (projectionType)
        {
        case Polydim::VEM::PCC::ProjectionTypes::PiNabla: {
            const std::vector<Eigen::MatrixXd> VanderDerivatives =
                utilities.ComputePolynomialsDerivativeValues(reference_element_data.Monomials,
                                                             monomials,
                                                             localSpace.inertia_polygon.Diameter,
                                                             vander);

            basisFunctionsDerivativeValues.resize(localSpace.Dimension);
            for (unsigned short i = 0; i < localSpace.Dimension; ++i)
                basisFunctionsDerivativeValues[i] = VanderDerivatives[i] * localSpace.PiNabla;
        }
        break;
        case Polydim::VEM::PCC::ProjectionTypes::Pi0km1Der: {
            basisFunctionsDerivativeValues.resize(localSpace.Dimension);
            for (unsigned short i = 0; i < localSpace.Dimension; ++i)
                basisFunctionsDerivativeValues[i] = vander.leftCols(localSpace.Nkm1) * localSpace.Pi0km1Der[i];
        }
        break;
        default:
            throw std::runtime_error("Unknown projector type");
        }

        std::vector<Eigen::MatrixXd> fmatrixInvTranspTimesBasisFunctionDerivativeValues2D(
            localSpace.Dimension,
            Eigen::MatrixXd::Zero(basisFunctionsDerivativeValues[0].rows(), basisFunctionsDerivativeValues[1].cols()));
        for (unsigned int d1 = 0; d1 < localSpace.Dimension; d1++)
        {
            for (unsigned int d2 = 0; d2 < localSpace.Dimension; d2++)
            {
                fmatrixInvTranspTimesBasisFunctionDerivativeValues2D[d1] +=
                    FmatrixInvTransp(d1, d2) * basisFunctionsDerivativeValues[d2];
            }
        }

        return fmatrixInvTranspTimesBasisFunctionDerivativeValues2D;
    }

    inline Eigen::MatrixXd ComputePolynomialsValues(const Polydim::VEM::PCC::VEM_PCC_2D_LocalSpace_Data &localSpace) const
    {
        return localSpace.VanderInternal;
    }

    inline Eigen::MatrixXd ComputePolynomialsValues(const Polydim::VEM::PCC::VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                    const Polydim::VEM::PCC::VEM_PCC_2D_LocalSpace_Data &localSpace,
                                                    const Eigen::MatrixXd &points) const
    {
        const Eigen::MatrixXd referencePoints =
            localSpace.inertia_data.FmatrixInv * (points.colwise() - localSpace.inertia_data.translation);
        return utilities.ComputePolynomialsValues(reference_element_data.Monomials,
                                                  monomials,
                                                  localSpace.inertia_polygon.Centroid,
                                                  localSpace.inertia_polygon.Diameter,
                                                  referencePoints);
    }

    inline std::vector<Eigen::MatrixXd> ComputePolynomialsDerivativeValues(const Polydim::VEM::PCC::VEM_PCC_2D_LocalSpace_Data &localSpace) const
    {
        const std::vector<Eigen::MatrixXd> &polynomialDerivatives = localSpace.VanderInternalDerivatives;
        const Eigen::MatrixXd FmatrixInvTransp = localSpace.inertia_data.FmatrixInv.transpose();

        std::vector<Eigen::MatrixXd> fmatrixInvTranspTimesPolynomialDerivatives(
            localSpace.Dimension,
            Eigen::MatrixXd::Zero(polynomialDerivatives[0].rows(), polynomialDerivatives[1].cols()));
        for (unsigned int d1 = 0; d1 < localSpace.Dimension; d1++)
        {
            for (unsigned int d2 = 0; d2 < localSpace.Dimension; d2++)
            {
                fmatrixInvTranspTimesPolynomialDerivatives[d1] += FmatrixInvTransp(d1, d2) * polynomialDerivatives[d2];
            }
        }
        return fmatrixInvTranspTimesPolynomialDerivatives;
    }

    inline std::vector<Eigen::MatrixXd> ComputePolynomialsDerivativeValues(const Polydim::VEM::PCC::VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                                           const Polydim::VEM::PCC::VEM_PCC_2D_LocalSpace_Data &localSpace,
                                                                           const Eigen::MatrixXd &points) const
    {
        const std::vector<Eigen::MatrixXd> polynomialDerivatives =
            utilities.ComputePolynomialsDerivativeValues(reference_element_data.Monomials,
                                                         monomials,
                                                         localSpace.inertia_polygon.Diameter,
                                                         ComputePolynomialsValues(reference_element_data, localSpace, points));

        const Eigen::MatrixXd FmatrixInvTransp = localSpace.inertia_data.FmatrixInv.transpose();

        std::vector<Eigen::MatrixXd> fmatrixInvTranspTimesPolynomialDerivatives(
            localSpace.Dimension,
            Eigen::MatrixXd::Zero(polynomialDerivatives[0].rows(), polynomialDerivatives[1].cols()));
        for (unsigned int d1 = 0; d1 < localSpace.Dimension; d1++)
        {
            for (unsigned int d2 = 0; d2 < localSpace.Dimension; d2++)
            {
                fmatrixInvTranspTimesPolynomialDerivatives[d1] += FmatrixInvTransp(d1, d2) * polynomialDerivatives[d2];
            }
        }
        return fmatrixInvTranspTimesPolynomialDerivatives;
    }

    inline Eigen::MatrixXd ComputeValuesOnEdge(const Polydim::VEM::PCC::VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                               const Eigen::VectorXd &pointsCurvilinearCoordinates) const
    {
        Eigen::VectorXd edgeInternalPoints;
        if (reference_element_data.Quadrature.ReferenceEdgeDOFsInternalPoints.rows() > 0)
            edgeInternalPoints = reference_element_data.Quadrature.ReferenceEdgeDOFsInternalPoints.row(0).transpose();
        const Eigen::VectorXd edgeBasisCoefficients =
            utilities.ComputeEdgeBasisCoefficients(reference_element_data.Order, edgeInternalPoints);

        return utilities.ComputeValuesOnEdge(edgeInternalPoints.transpose(), reference_element_data.Order, edgeBasisCoefficients, pointsCurvilinearCoordinates);
    }

    Eigen::MatrixXd ComputeBasisFunctionsLaplacianValues(const Polydim::VEM::PCC::VEM_PCC_2D_LocalSpace_Data &localSpace,
                                                         const Polydim::VEM::PCC::ProjectionTypes &projectionType) const
    {
        const Eigen::MatrixXd FmatrixInvTimesFmatrixInvTransp =
            localSpace.inertia_data.FmatrixInv * localSpace.inertia_data.FmatrixInv.transpose();

        switch (projectionType)
        {
        case Polydim::VEM::PCC::ProjectionTypes::Pi0km1Der: {
            Eigen::MatrixXd basisFunctionsLaplacianValues =
                Eigen::MatrixXd::Zero(localSpace.VanderInternalDerivatives[0].rows(), localSpace.NumBasisFunctions);

            for (unsigned short d1 = 0; d1 < localSpace.Dimension; d1++)
                for (unsigned short d2 = 0; d2 < localSpace.Dimension; d2++)
                    basisFunctionsLaplacianValues += FmatrixInvTimesFmatrixInvTransp(d2, d1) *
                                                     localSpace.VanderInternalDerivatives[d1].leftCols(localSpace.Nkm1) *
                                                     localSpace.Pi0km1Der[d2];

            return basisFunctionsLaplacianValues;
        }
        default:
            throw std::runtime_error("Unknown projector type");
        }
    }

    Eigen::MatrixXd ComputeBasisFunctionsLaplacianValues(const Polydim::VEM::PCC::VEM_PCC_2D_ReferenceElement_Data &,
                                                         const Polydim::VEM::PCC::VEM_PCC_2D_LocalSpace_Data &,
                                                         const Polydim::VEM::PCC::ProjectionTypes &,
                                                         const Eigen::MatrixXd &) const
    {
        throw std::runtime_error("Unimplemented method");
    }
    Eigen::MatrixXd ComputePolynomialsLaplacianValues(const Polydim::VEM::PCC::VEM_PCC_2D_ReferenceElement_Data &,
                                                      const Polydim::VEM::PCC::VEM_PCC_2D_LocalSpace_Data &,
                                                      const Eigen::MatrixXd &) const
    {
        throw std::runtime_error("Unimplemented method");
    }
};
} // namespace PCC
} // namespace VEM
} // namespace Polydim

#endif
