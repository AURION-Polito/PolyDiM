#ifndef __VEM_PCC_3D_Ortho_LocalSpace_HPP
#define __VEM_PCC_3D_Ortho_LocalSpace_HPP

#include "Eigen/Eigen"
#include "VEM_Monomials_2D.hpp"
#include "VEM_PCC_3D_LocalSpace_Data.hpp"
#include "VEM_PCC_3D_ReferenceElement.hpp"
#include "VEM_PCC_Utilities.hpp"
#include <vector>

namespace Polydim
{
namespace VEM
{
namespace PCC
{
class VEM_PCC_3D_Ortho_LocalSpace final
{
private:
    VEM_PCC_Utilities<2> utilities;
    Monomials::VEM_Monomials_2D monomials;

    /// \brief Initialize quantities required for computing projectors.
    /// \details This method computes \ref measure, \ref diameter, \ref
    /// vanderInternal, \ref vanderInternalDerivatives, \ref
    /// internalWeights, \ref vanderBoundary, \ref
    /// vanderBoundaryDerivatives, \ref boundaryWeights, \ref
    /// boundaryWeightsTimesNormal, \ref Hmatrix.
    /// \param geometry The geometry used as domain for the computation of projectors. It has to be
    /// an object of class \ref polyhedron with dimension 2.
    /// \note The following methods have to be called on the geometry before calling this method:
    ///  - \ref polyhedron::Compute2DpolyhedronProperties()
    ///  - \ref polyhedron::ComputePositionPoint()
    ///  - \ref polyhedron::ComputeNormalSign()
    ///  .
    /// Moreover, the method \ref Segment::ComputeNormal() has to be called on each edge of the
    /// geometry.
    void InitializeProjectorsComputation(const VEM_PCC_3D_ReferenceElement_Data& reference_element_data,
                                         const Eigen::MatrixXd& polyhedronVertices,
                                         const Eigen::Vector3d& polyhedronCentroid,
                                         const double& polyhedronDiameter,
                                         const Eigen::MatrixXd& internalQuadraturePoints,
                                         const Eigen::VectorXd& internalQuadratureWeights,
                                         const Eigen::MatrixXd& boundaryQuadraturePoints,
                                         VEM_PCC_3D_LocalSpace_Data& localSpace) const;

    /// \brief Compute matrix \ref piNabla.
    /// \note This requires \ref InitializeProjectorsComputation() to be
    /// called previously.
    void ComputePiNabla(const VEM_PCC_3D_ReferenceElement_Data& reference_element_data,
                        const double& polyhedronMeasure,
                        const double& polyhedronDiameter,
                        const Eigen::VectorXd& internalQuadratureWeights,
                        const Eigen::VectorXd& boundaryQuadratureWeights,
                        const std::vector<Eigen::VectorXd>& boundaryQuadratureWeightsTimesNormal,
                        VEM_PCC_3D_LocalSpace_Data& localSpace) const;

    /// \brief Compute matrices \ref pi0km1Der.
    void ComputeL2ProjectorsOfDerivatives(const VEM_PCC_3D_ReferenceElement_Data& reference_element_data,
                                          const double& polyhedronMeasure,
                                          const double& polyhedronDiameter,
                                          const std::vector<Eigen::VectorXd>& boundaryQuadratureWeightsTimesNormal,
                                          VEM_PCC_3D_LocalSpace_Data& localSpace) const;

    /// \brief Compute matrices \ref pi0km1 and \ref pi0k.
    /// \note This requires \ref ComputePiNabla() to be called previously.
    void ComputeL2Projectors(const double& polyhedronMeasure,
                             VEM_PCC_3D_LocalSpace_Data& localSpace) const
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


    /// \brief Compute the stabilization matrix with PiNabla projector.
    /// \note used for method with stabilization
    /// \return stabilization matrix, size numQuadraturePoints x NumberBasisFunctions()
    inline void ComputeStabilizationMatrix(const double& polyhedronDiameter,
                                           VEM_PCC_3D_LocalSpace_Data& localSpace) const
    {
        localSpace.StabMatrix = utilities.ComputeStabilizationMatrix(localSpace.PiNabla,
                                                                     polyhedronDiameter,
                                                                     localSpace.Dmatrix);
    }
    /// \brief Compute matrix \ref stabMatrix with Pi0k projector.
    /// \note This requires \ref ComputeL2Projectors() to be called previously.
    inline void ComputeStabilizationMatrixPi0k(const double& polyhedronMeasure,
                                               VEM_PCC_3D_LocalSpace_Data& localSpace) const
    {
        localSpace.StabMatrixPi0k = utilities.ComputeStabilizationMatrixPi0k(localSpace.Pi0k,
                                                                             polyhedronMeasure,
                                                                             localSpace.Dmatrix);
    }

    /// Compute the change of basis matrix and the mass matrix of orthogonal polynomial basis
    void ChangeOfBasis(const Eigen::VectorXd& internalQuadratureWeights,
                       VEM_PCC_3D_LocalSpace_Data& localSpace) const;

public:
    VEM_PCC_3D_LocalSpace_Data CreateLocalSpace(const VEM_PCC_3D_ReferenceElement_Data& reference_element_data,
                                                const VEM_PCC_3D_Polyhedron_Geometry& polyhedron) const;


    void ComputePolynomialsDofs(const double& polytopeMeasure,
                                VEM_PCC_3D_LocalSpace_Data& localSpace) const;

    inline Eigen::MatrixXd ComputeBasisFunctionsValues(const VEM_PCC_3D_LocalSpace_Data& localSpace,
                                                       const ProjectionTypes& projectionType) const
    {
        switch (projectionType)
        {
        case ProjectionTypes::Pi0km1:
            return localSpace.VanderInternal.leftCols(localSpace.Nkm1) *
                   localSpace.Qmatrix.topLeftCorner(localSpace.Nkm1,
                                                    localSpace.Nkm1).transpose() *
                   localSpace.Pi0km1;
        case ProjectionTypes::Pi0k:
            return localSpace.VanderInternal *
                   localSpace.Qmatrix.transpose() *
                   localSpace.Pi0k;
        default:
            throw std::runtime_error("Unsupported projectionType");
        }
    }

    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const VEM_PCC_3D_LocalSpace_Data& localSpace,
                                                                              const ProjectionTypes& projectionType) const
    {
        switch (projectionType)
        {
        case ProjectionTypes::Pi0km1Der:
        {
            std::vector<Eigen::MatrixXd> basisFunctionsDerivativeValues(localSpace.Dimension);
            for(unsigned short i = 0; i < localSpace.Dimension; ++i)
                basisFunctionsDerivativeValues[i] = localSpace.VanderInternal.leftCols(localSpace.Nkm1) *
                                                    localSpace.Qmatrix.topLeftCorner(localSpace.Nkm1,
                                                                                     localSpace.Nkm1).transpose() *
                                                    localSpace.Pi0km1Der[i];

            return basisFunctionsDerivativeValues;
        }
        case ProjectionTypes::PiNabla:
        {
            std::vector<Eigen::MatrixXd> basisFunctionDerivativeValues(localSpace.Dimension);
            for(unsigned short i = 0; i < localSpace.Dimension; ++i)
                basisFunctionDerivativeValues[i] = localSpace.VanderInternalDerivatives[i]
                                                   * localSpace.Qmatrix.transpose()
                                                   * localSpace.PiNabla;

            return basisFunctionDerivativeValues;
        }
        default:
            throw std::runtime_error("Unsupported projectionType");
        }
    }


    inline Eigen::MatrixXd ComputeBasisFunctionValues(const VEM_PCC_3D_ReferenceElement_Data& reference_element_data,
                                                      const VEM_PCC_3D_Polyhedron_Geometry& polyhedron,
                                                      const VEM_PCC_3D_LocalSpace_Data& localSpace,
                                                      const ProjectionTypes& projectionType,
                                                      const Eigen::MatrixXd& points) const
    {
        const Eigen::MatrixXd vanderInternal = ComputePolynomialsValues(reference_element_data,
                                                                        localSpace,
                                                                        polyhedron,
                                                                        points);


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

    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionDerivativeValues(const VEM_PCC_3D_ReferenceElement_Data& reference_element_data,
                                                                             const VEM_PCC_3D_Polyhedron_Geometry& polyhedron,
                                                                             const VEM_PCC_3D_LocalSpace_Data& localSpace,
                                                                             const ProjectionTypes& projectionType,
                                                                             const Eigen::MatrixXd& points) const
    {
        const std::vector<Eigen::MatrixXd> polynomialsDerivativeValues = ComputePolynomialsDerivativeValues(reference_element_data,
                                                                                                            localSpace,
                                                                                                            polyhedron,
                                                                                                            points);

        switch (projectionType)
        {
        case ProjectionTypes::Pi0km1Der:
        {
            std::vector<Eigen::MatrixXd> basisFunctionsDerivativeValues(localSpace.Dimension);
            for(unsigned short i = 0; i < localSpace.Dimension; ++i)
                basisFunctionsDerivativeValues[i] = polynomialsDerivativeValues[i] *
                                                    localSpace.Pi0km1Der[i];

            return basisFunctionsDerivativeValues;
        }
        case ProjectionTypes::PiNabla:
        {
            std::vector<Eigen::MatrixXd> basisFunctionDerivativeValues(localSpace.Dimension);
            for(unsigned short i = 0; i < localSpace.Dimension; ++i)
                basisFunctionDerivativeValues[i] = polynomialsDerivativeValues[i]
                                                   * localSpace.PiNabla;

            return basisFunctionDerivativeValues;
        }
        default:
            throw std::runtime_error("Unsupported projectionType");
        }
    }

    inline Eigen::MatrixXd ComputePolynomialsValues(const VEM_PCC_3D_LocalSpace_Data& localSpace) const
    {
        return localSpace.VanderInternal * localSpace.Qmatrix.transpose();
    }

    inline Eigen::MatrixXd ComputePolynomialsValues(const VEM_PCC_3D_ReferenceElement_Data& reference_element_data,
                                                    const VEM_PCC_3D_LocalSpace_Data& localSpace,
                                                    const VEM_PCC_3D_Polyhedron_Geometry& polyhedron,
                                                    const Eigen::MatrixXd& points) const
    {
        return monomials.Vander(reference_element_data.Monomials,
                                points,
                                polyhedron.Centroid,
                                polyhedron.Diameter) *
               localSpace.Qmatrix.transpose();
    }

    inline std::vector<Eigen::MatrixXd> ComputePolynomialsDerivativeValues(const VEM_PCC_3D_LocalSpace_Data& localSpace) const
    {
        std::vector<Eigen::MatrixXd> polynomialBasisDerivativeValues(localSpace.Dimension);
        for(unsigned short i = 0; i < localSpace.Dimension; ++i)
            polynomialBasisDerivativeValues[i] = localSpace.VanderInternalDerivatives[i] *
                                                 localSpace.Qmatrix.topLeftCorner(localSpace.Nkm1,
                                                                                  localSpace.Nkm1).transpose();
        return polynomialBasisDerivativeValues;
    }

    inline std::vector<Eigen::MatrixXd> ComputePolynomialsDerivativeValues(const VEM_PCC_3D_ReferenceElement_Data& reference_element_data,
                                                                           const VEM_PCC_3D_LocalSpace_Data& localSpace,
                                                                           const VEM_PCC_3D_Polyhedron_Geometry& polyhedron,
                                                                           const Eigen::MatrixXd& points) const
    {
        const Eigen::MatrixXd monomialBasisValues = monomials.Vander(reference_element_data.Monomials,
                                                                     points,
                                                                     polyhedron.Centroid,
                                                                     polyhedron.Diameter);
        const std::vector<Eigen::MatrixXd> monomialBasisDerivativeValues = monomials.VanderDerivatives(reference_element_data.Monomials,
                                                                                                       monomialBasisValues,
                                                                                                       polyhedron.Diameter);

        std::vector<Eigen::MatrixXd> polynomialBasisDerivativeValues(localSpace.Dimension);
        for(unsigned short i = 0; i < localSpace.Dimension; ++i)
            polynomialBasisDerivativeValues[i] = monomialBasisDerivativeValues[i] *
                                                 localSpace.Qmatrix.topLeftCorner(localSpace.Nkm1,
                                                                                  localSpace.Nkm1).transpose();
        return polynomialBasisDerivativeValues;
    }

    inline Eigen::MatrixXd ComputeValuesOnEdge(const VEM_PCC_3D_LocalSpace_Data& localSpace,
                                               const Eigen::VectorXd& edgeInternalPoints,
                                               const Eigen::VectorXd& pointsCurvilinearCoordinates) const
    {
        const Eigen::VectorXd edgeBasisCoefficients = utilities.ComputeEdgeBasisCoefficients(localSpace.Order,
                                                                                             edgeInternalPoints);
        return utilities.ComputeValuesOnEdge(edgeInternalPoints.transpose(),
                                             localSpace.Order,
                                             edgeBasisCoefficients,
                                             pointsCurvilinearCoordinates);
    }
};
}
}
}

#endif
