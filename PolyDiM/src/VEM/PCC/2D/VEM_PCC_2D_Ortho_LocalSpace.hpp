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

#ifndef __VEM_PCC_2D_Ortho_LocalSpace_HPP
#define __VEM_PCC_2D_Ortho_LocalSpace_HPP

#include "I_VEM_PCC_2D_LocalSpace.hpp"
#include "VEM_Monomials_2D.hpp"

namespace Polydim
{
namespace VEM
{
namespace PCC
{
/// \brief Interface class for Primal Conforming Constant degree 2D Virtual Element Methods \cite Mascotto2018 \cite
/// Teora2024.
class VEM_PCC_2D_Ortho_LocalSpace final : public I_VEM_PCC_2D_LocalSpace
{
  private:
    VEM_PCC_Utilities<2> utilities;
    Utilities::VEM_Monomials_2D monomials;

    void InitializeProjectorsComputation(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                         const Eigen::MatrixXd &polygonVertices,
                                         const Eigen::Vector3d &polygonCentroid,
                                         const double &polygonMeasure,
                                         const double &polygonDiameter,
                                         const Eigen::MatrixXd &internalQuadraturePoints,
                                         const Eigen::VectorXd &internalQuadratureWeights,
                                         const Eigen::MatrixXd &boundaryQuadraturePoints,
                                         VEM_PCC_2D_LocalSpace_Data &localSpace) const;

    void ComputePiNabla(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                        const double &polygonMeasure,
                        const double &polygonDiameter,
                        const Eigen::VectorXd &internalQuadratureWeights,
                        const Eigen::VectorXd &boundaryQuadratureWeights,
                        const std::vector<Eigen::VectorXd> &boundaryQuadratureWeightsTimesNormal,
                        VEM_PCC_2D_LocalSpace_Data &localSpace) const;

    void ComputeL2ProjectorsOfDerivatives(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                          const double &polygonMeasure,
                                          const double &polygonDiameter,
                                          const std::vector<Eigen::VectorXd> &boundaryQuadratureWeightsTimesNormal,
                                          VEM_PCC_2D_LocalSpace_Data &localSpace) const;

    void ComputeL2Projectors(const double &polygonMeasure, VEM_PCC_2D_LocalSpace_Data &localSpace) const
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

  public:
    VEM_PCC_2D_LocalSpace_Data CreateLocalSpace(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                const VEM_PCC_2D_Polygon_Geometry &polygon) const;

    VEM_PCC_2D_LocalSpace_Data Compute3DUtilities(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                  const VEM_PCC_2D_Polygon_Geometry &polygon) const;

    inline Eigen::MatrixXd ComputeDofiDofiStabilizationMatrix(const VEM_PCC_2D_LocalSpace_Data &localSpace,
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

    void ComputePolynomialsDofs(const double &polytopeMeasure, VEM_PCC_2D_LocalSpace_Data &localSpace) const;

    inline Eigen::MatrixXd ComputeBasisFunctionsValues(const VEM_PCC_2D_LocalSpace_Data &localSpace,
                                                       const ProjectionTypes &projectionType) const
    {
        switch (projectionType)
        {
        case ProjectionTypes::Pi0km1:
            return localSpace.VanderInternal.leftCols(localSpace.Nkm1) *
                   localSpace.Qmatrix.topLeftCorner(localSpace.Nkm1, localSpace.Nkm1).transpose() * localSpace.Pi0km1;
        case ProjectionTypes::Pi0k:
            return localSpace.VanderInternal * localSpace.Qmatrix.transpose() * localSpace.Pi0k;
        default:
            throw std::runtime_error("Unsupported projectionType");
        }
    }

    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const VEM_PCC_2D_LocalSpace_Data &localSpace,
                                                                              const ProjectionTypes &projectionType) const
    {
        switch (projectionType)
        {
        case ProjectionTypes::Pi0km1Der: {
            std::vector<Eigen::MatrixXd> basisFunctionsDerivativeValues(localSpace.Dimension);
            for (unsigned short i = 0; i < localSpace.Dimension; ++i)
                basisFunctionsDerivativeValues[i] =
                    localSpace.VanderInternal.leftCols(localSpace.Nkm1) *
                    localSpace.Qmatrix.topLeftCorner(localSpace.Nkm1, localSpace.Nkm1).transpose() * localSpace.Pi0km1Der[i];

            return basisFunctionsDerivativeValues;
        }
        case ProjectionTypes::PiNabla: {
            std::vector<Eigen::MatrixXd> basisFunctionDerivativeValues(localSpace.Dimension);
            for (unsigned short i = 0; i < localSpace.Dimension; ++i)
                basisFunctionDerivativeValues[i] =
                    localSpace.VanderInternalDerivatives[i] * localSpace.Qmatrix.transpose() * localSpace.PiNabla;

            return basisFunctionDerivativeValues;
        }
        default:
            throw std::runtime_error("Unsupported projectionType");
        }
    }

    inline Eigen::MatrixXd ComputeBasisFunctionsValues(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                       const VEM_PCC_2D_LocalSpace_Data &localSpace,
                                                       const ProjectionTypes &projectionType,
                                                       const Eigen::MatrixXd &points) const
    {
        const Eigen::MatrixXd vanderInternal = ComputePolynomialsValues(reference_element_data, localSpace, points);

        switch (projectionType)
        {
        case ProjectionTypes::Pi0km1:
            return vanderInternal.leftCols(localSpace.Nkm1) * localSpace.Pi0km1;
        case ProjectionTypes::Pi0k:
            return vanderInternal * localSpace.Pi0k;
        default:
            throw std::runtime_error("Unsupported projectionType");
        }
    }

    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                                              const VEM_PCC_2D_LocalSpace_Data &localSpace,
                                                                              const ProjectionTypes &projectionType,
                                                                              const Eigen::MatrixXd &points) const
    {
        switch (projectionType)
        {
        case ProjectionTypes::Pi0km1Der: {

            const Eigen::MatrixXd polynomialsValues = ComputePolynomialsValues(reference_element_data, localSpace, points);

            std::vector<Eigen::MatrixXd> basisFunctionsDerivativeValues(localSpace.Dimension);
            for (unsigned short i = 0; i < localSpace.Dimension; ++i)
                basisFunctionsDerivativeValues[i] = polynomialsValues.leftCols(localSpace.Nkm1) * localSpace.Pi0km1Der[i];

            return basisFunctionsDerivativeValues;
        }
        case ProjectionTypes::PiNabla: {

            const std::vector<Eigen::MatrixXd> polynomialsDerivativeValues =
                ComputePolynomialsDerivativeValues(reference_element_data, localSpace, points);

            std::vector<Eigen::MatrixXd> basisFunctionDerivativeValues(localSpace.Dimension);
            for (unsigned short i = 0; i < localSpace.Dimension; ++i)
                basisFunctionDerivativeValues[i] = polynomialsDerivativeValues[i] * localSpace.PiNabla;

            return basisFunctionDerivativeValues;
        }
        default:
            throw std::runtime_error("Unsupported projectionType");
        }
    }

    inline Eigen::MatrixXd ComputePolynomialsValues(const VEM_PCC_2D_LocalSpace_Data &localSpace) const
    {
        return localSpace.VanderInternal * localSpace.Qmatrix.transpose();
    }

    inline Eigen::MatrixXd ComputePolynomialsValues(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                    const VEM_PCC_2D_LocalSpace_Data &localSpace,
                                                    const Eigen::MatrixXd &points) const
    {
        return monomials.Vander(reference_element_data.Monomials, points, localSpace.Centroid, localSpace.Diameter) *
               localSpace.Qmatrix.transpose();
    }

    inline std::vector<Eigen::MatrixXd> ComputePolynomialsDerivativeValues(const VEM_PCC_2D_LocalSpace_Data &localSpace) const
    {
        std::vector<Eigen::MatrixXd> polynomialBasisDerivativeValues(localSpace.Dimension);
        for (unsigned short i = 0; i < localSpace.Dimension; ++i)
            polynomialBasisDerivativeValues[i] = localSpace.VanderInternalDerivatives[i] * localSpace.Qmatrix.transpose();
        return polynomialBasisDerivativeValues;
    }

    inline std::vector<Eigen::MatrixXd> ComputePolynomialsDerivativeValues(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                                           const VEM_PCC_2D_LocalSpace_Data &localSpace,
                                                                           const Eigen::MatrixXd &points) const
    {
        const Eigen::MatrixXd monomialBasisValues =
            monomials.Vander(reference_element_data.Monomials, points, localSpace.Centroid, localSpace.Diameter);
        const std::vector<Eigen::MatrixXd> monomialBasisDerivativeValues =
            monomials.VanderDerivatives(reference_element_data.Monomials, monomialBasisValues, localSpace.Diameter);

        std::vector<Eigen::MatrixXd> polynomialBasisDerivativeValues(localSpace.Dimension);
        for (unsigned short i = 0; i < localSpace.Dimension; ++i)
            polynomialBasisDerivativeValues[i] = monomialBasisDerivativeValues[i] * localSpace.Qmatrix.transpose();
        return polynomialBasisDerivativeValues;
    }

    inline Eigen::MatrixXd ComputeValuesOnEdge(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                               const Eigen::VectorXd &pointsCurvilinearCoordinates) const
    {
        Eigen::VectorXd edgeInternalPoints;
        if (reference_element_data.Quadrature.ReferenceEdgeDOFsInternalPoints.rows() > 0)
            edgeInternalPoints = reference_element_data.Quadrature.ReferenceEdgeDOFsInternalPoints.row(0).transpose();
        const Eigen::VectorXd edgeBasisCoefficients =
            utilities.ComputeEdgeBasisCoefficients(reference_element_data.Order, edgeInternalPoints);

        return utilities.ComputeValuesOnEdge(edgeInternalPoints.transpose(), reference_element_data.Order, edgeBasisCoefficients, pointsCurvilinearCoordinates);
    }

    Eigen::MatrixXd ComputeBasisFunctionsLaplacianValues(const VEM_PCC_2D_LocalSpace_Data &localSpace,
                                                         const ProjectionTypes &projectionType) const
    {
        switch (projectionType)
        {
        case ProjectionTypes::Pi0km1Der: {
            Eigen::MatrixXd basisFunctionsLaplacianValues =
                localSpace.VanderInternalDerivatives[0].leftCols(localSpace.Nkm1) *
                localSpace.Qmatrix.topLeftCorner(localSpace.Nkm1, localSpace.Nkm1).transpose() * localSpace.Pi0km1Der[0];
            for (unsigned int d = 1; d < localSpace.Dimension; ++d)
                basisFunctionsLaplacianValues +=
                    localSpace.VanderInternalDerivatives[d].leftCols(localSpace.Nkm1) *
                    localSpace.Qmatrix.topLeftCorner(localSpace.Nkm1, localSpace.Nkm1).transpose() * localSpace.Pi0km1Der[d];

            return basisFunctionsLaplacianValues;
        }
        default:
            throw std::runtime_error("Unknown projector type");
        }
    }
    Eigen::MatrixXd ComputeBasisFunctionsLaplacianValues(const VEM_PCC_2D_ReferenceElement_Data &,
                                                         const VEM_PCC_2D_LocalSpace_Data &,
                                                         const ProjectionTypes &,
                                                         const Eigen::MatrixXd &) const
    {
        throw std::runtime_error("Unimplemented method");
    }
    Eigen::MatrixXd ComputePolynomialsLaplacianValues(const VEM_PCC_2D_ReferenceElement_Data &,
                                                      const VEM_PCC_2D_LocalSpace_Data &,
                                                      const Eigen::MatrixXd &) const
    {
        throw std::runtime_error("Unimplemented method");
    }
};
} // namespace PCC
} // namespace VEM
} // namespace Polydim

#endif
