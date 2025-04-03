#ifndef __VEM_PCC_3D_Ortho_LocalSpace_HPP
#define __VEM_PCC_3D_Ortho_LocalSpace_HPP

#include "Eigen/Eigen"
#include "I_VEM_PCC_2D_ReferenceElement.hpp"
#include "I_VEM_PCC_3D_LocalSpace.hpp"
#include "I_VEM_PCC_3D_ReferenceElement.hpp"
#include "VEM_Monomials_3D.hpp"
#include "VEM_PCC_2D_LocalSpace_Data.hpp"
#include "VEM_PCC_2D_Ortho_LocalSpace.hpp"
#include "VEM_PCC_3D_LocalSpace_Data.hpp"
#include "VEM_PCC_Utilities.hpp"
#include <vector>

namespace Polydim
{
namespace VEM
{
namespace PCC
{

/// \brief Interface class for Primal Conforming Constant degree 3D Virtual Element Methods \cite DassiMascotto2018
/// \cite Teora2024.
class VEM_PCC_3D_Ortho_LocalSpace final : public I_VEM_PCC_3D_LocalSpace
{
  private:
    VEM_PCC_Utilities<3> utilities;
    Monomials::VEM_Monomials_3D monomials;

    void InitializeProjectorsComputation(const VEM_PCC_3D_ReferenceElement_Data &reference_element_data,
                                         const Eigen::MatrixXd &polyhedronVertices,
                                         const Eigen::MatrixXi &polyhedronEdges,
                                         const std::vector<Eigen::MatrixXi> &polyhedronFaces,
                                         const Eigen::Vector3d &polyhedronCentroid,
                                         const double &polyhedronMeasure,
                                         const double &polyhedronDiameter,
                                         const Eigen::MatrixXd &internalQuadraturePoints,
                                         const Eigen::VectorXd &internalQuadratureWeights,
                                         const Eigen::MatrixXd &boundaryQuadraturePoints,
                                         const Eigen::MatrixXd &edgeInternalQuadraturePoints,
                                         VEM_PCC_3D_LocalSpace_Data &localSpace) const;

    void ComputePiNabla(const VEM_PCC_3D_ReferenceElement_Data &reference_element_data,
                        const double &polyhedronMeasure,
                        const double &polyhedronDiameter,
                        const Eigen::VectorXd &internalQuadratureWeights,
                        const Eigen::VectorXd &boundaryQuadratureWeights,
                        const std::vector<Eigen::VectorXd> &boundaryQuadratureWeightsTimesNormal,
                        VEM_PCC_3D_LocalSpace_Data &localSpace) const;

    void ComputeL2ProjectorsOfDerivatives(const VEM_PCC_3D_ReferenceElement_Data &reference_element_data,
                                          const double &polyhedronMeasure,
                                          const double &polyhedronDiameter,
                                          const std::vector<Eigen::VectorXd> &boundaryQuadratureWeightsTimesNormal,
                                          VEM_PCC_3D_LocalSpace_Data &localSpace) const;

    void ComputeL2Projectors(const double &polyhedronMeasure, VEM_PCC_3D_LocalSpace_Data &localSpace) const
    {
        utilities.ComputeL2Projectors(polyhedronMeasure,
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

    void ComputeFaceProjectors(const VEM_PCC_2D_Ortho_LocalSpace &faceVemValues,
                               const std::vector<Eigen::MatrixXi> &polyhedronFaces,
                               const std::vector<VEM_PCC_2D_Polygon_Geometry> &polygonalFaces,
                               const Eigen::MatrixXd &boundaryQuadraturePoints,
                               const Eigen::VectorXd &boundaryQuadratureWeights,
                               VEM_PCC_3D_LocalSpace_Data &localSpace) const;

    void ComputePolynomialsDofs(const double &polytopeMeasure, VEM_PCC_3D_LocalSpace_Data &localSpace) const;

  public:
    VEM_PCC_3D_LocalSpace_Data CreateLocalSpace(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data_2D,
                                                const VEM_PCC_3D_ReferenceElement_Data &reference_element_data_3D,
                                                const std::vector<VEM_PCC_2D_Polygon_Geometry> &polygonalFaces,
                                                const VEM_PCC_3D_Polyhedron_Geometry &polyhedron) const;

    inline Eigen::MatrixXd ComputeDofiDofiStabilizationMatrix(const VEM_PCC_3D_LocalSpace_Data &localSpace,
                                                              const ProjectionTypes &projectionType) const
    {
        switch (projectionType)
        {
        case ProjectionTypes::PiNabla:
            return utilities.ComputeDofiDofiStabilizationMatrix(localSpace.PiNabla, localSpace.Diameter, localSpace.Dmatrix);
        case ProjectionTypes::Pi0k:
            return utilities.ComputeDofiDofiStabilizationMatrix(localSpace.Pi0k, localSpace.Measure, localSpace.Dmatrix);
        default:
            throw std::runtime_error("not valid projection type");
        }
    }

    inline Eigen::MatrixXd ComputeBasisFunctionsValues(const VEM_PCC_3D_LocalSpace_Data &localSpace,
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

    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const VEM_PCC_3D_LocalSpace_Data &localSpace,
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

    inline Eigen::MatrixXd ComputeBasisFunctionsValues(const VEM_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                       const VEM_PCC_3D_LocalSpace_Data &localSpace,
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

    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const VEM_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                                              const VEM_PCC_3D_LocalSpace_Data &localSpace,
                                                                              const ProjectionTypes &projectionType,
                                                                              const Eigen::MatrixXd &points) const
    {
        const std::vector<Eigen::MatrixXd> polynomialsDerivativeValues =
            ComputePolynomialsDerivativeValues(reference_element_data, localSpace, points);

        switch (projectionType)
        {
        case ProjectionTypes::Pi0km1Der: {
            std::vector<Eigen::MatrixXd> basisFunctionsDerivativeValues(localSpace.Dimension);
            for (unsigned short i = 0; i < localSpace.Dimension; ++i)
                basisFunctionsDerivativeValues[i] = polynomialsDerivativeValues[i] * localSpace.Pi0km1Der[i];

            return basisFunctionsDerivativeValues;
        }
        case ProjectionTypes::PiNabla: {
            std::vector<Eigen::MatrixXd> basisFunctionDerivativeValues(localSpace.Dimension);
            for (unsigned short i = 0; i < localSpace.Dimension; ++i)
                basisFunctionDerivativeValues[i] = polynomialsDerivativeValues[i] * localSpace.PiNabla;

            return basisFunctionDerivativeValues;
        }
        default:
            throw std::runtime_error("Unsupported projectionType");
        }
    }

    inline Eigen::MatrixXd ComputePolynomialsValues(const VEM_PCC_3D_LocalSpace_Data &localSpace) const
    {
        return localSpace.VanderInternal * localSpace.Qmatrix.transpose();
    }

    inline Eigen::MatrixXd ComputePolynomialsValues(const VEM_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                    const VEM_PCC_3D_LocalSpace_Data &localSpace,
                                                    const Eigen::MatrixXd &points) const
    {
        return monomials.Vander(reference_element_data.Monomials, points, localSpace.Centroid, localSpace.Diameter) *
               localSpace.Qmatrix.transpose();
    }

    inline std::vector<Eigen::MatrixXd> ComputePolynomialsDerivativeValues(const VEM_PCC_3D_LocalSpace_Data &localSpace) const
    {
        std::vector<Eigen::MatrixXd> polynomialBasisDerivativeValues(localSpace.Dimension);
        for (unsigned short i = 0; i < localSpace.Dimension; ++i)
            polynomialBasisDerivativeValues[i] = localSpace.VanderInternalDerivatives[i] *
                                                 localSpace.Qmatrix.topLeftCorner(localSpace.Nkm1, localSpace.Nkm1).transpose();
        return polynomialBasisDerivativeValues;
    }

    inline std::vector<Eigen::MatrixXd> ComputePolynomialsDerivativeValues(const VEM_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                                           const VEM_PCC_3D_LocalSpace_Data &localSpace,
                                                                           const Eigen::MatrixXd &points) const
    {
        const Eigen::MatrixXd monomialBasisValues =
            monomials.Vander(reference_element_data.Monomials, points, localSpace.Centroid, localSpace.Diameter);
        const std::vector<Eigen::MatrixXd> monomialBasisDerivativeValues =
            monomials.VanderDerivatives(reference_element_data.Monomials, monomialBasisValues, localSpace.Diameter);

        std::vector<Eigen::MatrixXd> polynomialBasisDerivativeValues(localSpace.Dimension);
        for (unsigned short i = 0; i < localSpace.Dimension; ++i)
            polynomialBasisDerivativeValues[i] =
                monomialBasisDerivativeValues[i] * localSpace.Qmatrix.topLeftCorner(localSpace.Nkm1, localSpace.Nkm1).transpose();
        return polynomialBasisDerivativeValues;
    }

    inline Eigen::MatrixXd ComputeValuesOnEdge(const VEM_PCC_3D_LocalSpace_Data &localSpace,
                                               const Eigen::VectorXd &edgeInternalPoints,
                                               const Eigen::VectorXd &pointsCurvilinearCoordinates) const
    {
        const Eigen::VectorXd edgeBasisCoefficients = utilities.ComputeEdgeBasisCoefficients(localSpace.Order, edgeInternalPoints);
        return utilities.ComputeValuesOnEdge(edgeInternalPoints.transpose(), localSpace.Order, edgeBasisCoefficients, pointsCurvilinearCoordinates);
    }

    Eigen::MatrixXd ComputeBasisFunctionsLaplacianValues(const VEM_PCC_3D_LocalSpace_Data &, const ProjectionTypes &) const
    {
        throw std::runtime_error("Unimplemented method");
    }

    Eigen::MatrixXd ComputeBasisFunctionsLaplacianValues(const VEM_PCC_3D_ReferenceElement_Data &,
                                                         const VEM_PCC_3D_LocalSpace_Data &,
                                                         const ProjectionTypes &,
                                                         const Eigen::MatrixXd &) const
    {
        throw std::runtime_error("Unimplemented method");
    }

    Eigen::MatrixXd ComputePolynomialsLaplacianValues(const VEM_PCC_3D_ReferenceElement_Data &,
                                                      const VEM_PCC_3D_LocalSpace_Data &,
                                                      const Eigen::MatrixXd &) const
    {
        throw std::runtime_error("Unimplemented method");
    }
};
} // namespace PCC
} // namespace VEM
} // namespace Polydim

#endif
