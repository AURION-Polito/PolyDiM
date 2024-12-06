#ifndef __VEM_PCC_3D_LocalSpace_HPP
#define __VEM_PCC_3D_LocalSpace_HPP

#include "Eigen/Eigen"
#include "VEM_Monomials_3D.hpp"
#include "VEM_PCC_2D_LocalSpace.hpp"
#include "VEM_PCC_2D_ReferenceElement.hpp"
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
/// \brief Class used for computing values of basis functions of 2D
/// Primal Conforming Constant degree Virtual Element Methods.
class VEM_PCC_3D_LocalSpace final
{
private:
    VEM_PCC_Utilities<3> utilities;
    Monomials::VEM_Monomials_3D monomials;

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
                                         const Eigen::MatrixXi& polyhedronEdges,
                                         const std::vector<Eigen::MatrixXi>& polyhedronFaces,
                                         const Eigen::Vector3d& polyhedronCentroid,
                                         const double& polyhedronDiameter,
                                         const Eigen::MatrixXd& internalQuadraturePoints,
                                         const Eigen::VectorXd& internalQuadratureWeights,
                                         const Eigen::MatrixXd& boundaryQuadraturePoints,
                                         const Eigen::MatrixXd &edgeInternalQuadraturePoints,
                                         VEM_PCC_3D_LocalSpace_Data &localSpace) const;

    /// \brief Compute matrix \ref piNabla.
    /// \note This requires \ref InitializeProjectorsComputation() to be called previously.
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

    void ComputeFaceProjectors(const VEM_PCC_2D_LocalSpace &faceVemValues,
                               const std::vector<Eigen::MatrixXi> &polyhedronFaces,
                               const std::vector<double> &facesMeasure,
                               const Eigen::MatrixXd &boundaryQuadraturePoints,
                               const Eigen::VectorXd &boundaryQuadratureWeights,
                               VEM_PCC_3D_LocalSpace_Data &localSpace) const;
public:
    VEM_PCC_3D_LocalSpace_Data CreateLocalSpace(const VEM_PCC_2D_ReferenceElement_Data& reference_element_data_2D,
                                                const VEM_PCC_3D_ReferenceElement_Data &reference_element_data_3D,
                                                const VEM_PCC_3D_Polyhedron_Geometry &polyhedron) const;

    /// \brief Compute matrix D: D_{ij} = dof_i(m_j).
    void ComputePolynomialsDofs(const double& polytopeMeasure,
                                VEM_PCC_3D_LocalSpace_Data &localSpace) const;

    inline Eigen::MatrixXd ComputeBasisFunctionsValues(const VEM_PCC_3D_LocalSpace_Data& localSpace,
                                                       const ProjectionTypes& projectionType) const
    {
        return utilities.ComputeBasisFunctionsValues(projectionType,
                                                     localSpace.Nkm1,
                                                     localSpace.Pi0km1,
                                                     localSpace.Pi0k,
                                                     localSpace.VanderInternal);
    }

    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const VEM_PCC_3D_LocalSpace_Data& localSpace,
                                                                              const ProjectionTypes& projectionType) const
    {
        return utilities.ComputeBasisFunctionsDerivativeValues(projectionType,
                                                               localSpace.Nkm1,
                                                               localSpace.VanderInternal,
                                                               localSpace.VanderInternalDerivatives,
                                                               localSpace.PiNabla,
                                                               localSpace.Pi0km1Der);
    }

    inline Eigen::MatrixXd ComputeBasisFunctionsLaplacianValues(const VEM_PCC_3D_LocalSpace_Data& localSpace) const
    {
        return utilities.ComputeBasisFunctionsLaplacianValues(localSpace.Nkm1,
                                                              localSpace.VanderInternalDerivatives,
                                                              localSpace.Pi0km1Der);
    }

    inline Eigen::MatrixXd ComputeBasisFunctionsValues(const VEM_PCC_3D_ReferenceElement_Data& reference_element_data,
                                                       const VEM_PCC_3D_Polyhedron_Geometry& polyhedron,
                                                       const VEM_PCC_3D_LocalSpace_Data& localSpace,
                                                       const ProjectionTypes& projectionType,
                                                       const Eigen::MatrixXd& points) const
    {
        return utilities.ComputeBasisFunctionsValues(projectionType,
                                                     localSpace.Nkm1,
                                                     localSpace.Pi0km1,
                                                     localSpace.Pi0k,
                                                     ComputePolynomialsValues(reference_element_data,
                                                                              polyhedron,
                                                                              points));
    }

    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const VEM_PCC_3D_ReferenceElement_Data& reference_element_data,
                                                                              const VEM_PCC_3D_Polyhedron_Geometry& polyhedron,
                                                                              const VEM_PCC_3D_LocalSpace_Data& localSpace,
                                                                              const ProjectionTypes& projectionType,
                                                                              const Eigen::MatrixXd& points) const
    {
        return utilities.ComputeBasisFunctionsDerivativeValues(projectionType,
                                                               localSpace.Nkm1,
                                                               ComputePolynomialsValues(reference_element_data,
                                                                                        polyhedron,
                                                                                        points),
                                                               ComputePolynomialsDerivativeValues(reference_element_data,
                                                                                                  polyhedron,
                                                                                                  points),
                                                               localSpace.PiNabla,
                                                               localSpace.Pi0km1Der);
    }

    inline Eigen::MatrixXd ComputeBasisFunctionsLaplacianValues(const VEM_PCC_3D_ReferenceElement_Data& reference_element_data,
                                                                const VEM_PCC_3D_Polyhedron_Geometry& polyhedron,
                                                                const VEM_PCC_3D_LocalSpace_Data& localSpace,
                                                                const Eigen::MatrixXd& points) const
    {
        return utilities.ComputeBasisFunctionsLaplacianValues(localSpace.Nkm1,
                                                              localSpace.Pi0km1Der,
                                                              ComputePolynomialsDerivativeValues(reference_element_data,
                                                                                                 polyhedron,
                                                                                                 points));
    }

    inline Eigen::MatrixXd ComputePolynomialsValues(const VEM_PCC_3D_LocalSpace_Data& localSpace) const
    {
        return utilities.ComputePolynomialsValues(localSpace.VanderInternal);
    }

    inline Eigen::MatrixXd ComputePolynomialsValues(const VEM_PCC_3D_ReferenceElement_Data& reference_element_data,
                                                    const VEM_PCC_3D_Polyhedron_Geometry& polyhedron,
                                                    const Eigen::MatrixXd& points) const
    {
        return utilities.ComputePolynomialsValues(reference_element_data.Monomials,
                                                  monomials,
                                                  polyhedron.Centroid,
                                                  polyhedron.Diameter,
                                                  points);
    }

    inline std::vector<Eigen::MatrixXd> ComputePolynomialsDerivativeValues(const VEM_PCC_3D_LocalSpace_Data& localSpace) const
    {
        return utilities.ComputePolynomialsDerivativeValues(localSpace.VanderInternalDerivatives);
    }

    inline std::vector<Eigen::MatrixXd> ComputePolynomialsDerivativeValues(const VEM_PCC_3D_ReferenceElement_Data& reference_element_data,
                                                                           const VEM_PCC_3D_Polyhedron_Geometry& polyhedron,
                                                                           const Eigen::MatrixXd& points) const
    {
        return utilities.ComputePolynomialsDerivativeValues(reference_element_data.Monomials,
                                                            monomials,
                                                            polyhedron.Diameter,
                                                            ComputePolynomialsValues(reference_element_data,
                                                                                     polyhedron,
                                                                                     points));
    }

    inline Eigen::MatrixXd ComputePolynomialsLaplacianValues(const VEM_PCC_3D_ReferenceElement_Data& reference_element_data,
                                                             const VEM_PCC_3D_Polyhedron_Geometry& polyhedron,
                                                             const Eigen::MatrixXd& points) const
    {
        return utilities.ComputePolynomialsLaplacianValues(reference_element_data.Monomials,
                                                           monomials,
                                                           polyhedron.Diameter,
                                                           ComputePolynomialsValues(reference_element_data,
                                                                                    polyhedron,
                                                                                    points));
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
