#ifndef __VEM_PCC_3D_Inertia_LocalSpace_HPP
#define __VEM_PCC_3D_Inertia_LocalSpace_HPP

#include "Eigen/Eigen"
#include "I_VEM_PCC_2D_ReferenceElement.hpp"
#include "I_VEM_PCC_3D_LocalSpace.hpp"
#include "I_VEM_PCC_3D_ReferenceElement.hpp"
#include "VEM_Monomials_3D.hpp"
#include "VEM_PCC_2D_Inertia_LocalSpace.hpp"
#include "VEM_PCC_2D_LocalSpace_Data.hpp"
#include "VEM_PCC_3D_LocalSpace_Data.hpp"
#include "VEM_PCC_Utilities.hpp"
#include <vector>

namespace Polydim
{
namespace VEM
{
namespace PCC
{
/// \brief Class used for computing values of basis functions of 3D
/// Primal Conforming Constant degree Virtual Element Methods.
class VEM_PCC_3D_Inertia_LocalSpace final : public I_VEM_PCC_3D_LocalSpace
{
  private:
    VEM_PCC_Utilities<3> utilities;
    Monomials::VEM_Monomials_3D monomials;

    void InitializeProjectorsComputation(const VEM_PCC_3D_ReferenceElement_Data &reference_element_data,
                                         const Eigen::MatrixXd &polyhedronVertices,
                                         const Eigen::MatrixXi &polyhedronEdges,
                                         const std::vector<Eigen::MatrixXi> &polyhedronFaces,
                                         const Eigen::Vector3d &polyhedronCentroid,
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

    void ComputeFaceProjectors(const VEM_PCC_2D_Inertia_LocalSpace &faceVemValues,
                               const std::vector<Eigen::MatrixXi> &polyhedronFaces,
                               const std::vector<double> &polygonalFaces,
                               const Eigen::MatrixXd &boundaryQuadraturePoints,
                               const Eigen::VectorXd &boundaryQuadratureWeights,
                               VEM_PCC_3D_LocalSpace_Data &localSpace) const;

    void ComputePolynomialsDofs(const double &polytopeMeasure, VEM_PCC_3D_LocalSpace_Data &localSpace) const;

    void InertiaMapping(const Gedim::GeometryUtilities &geometryUtilities,
                        const VEM_PCC_3D_Polyhedron_Geometry &polyhedron,
                        VEM_PCC_3D_Inertia_Data &inertia_data) const;

    void ComputeGeometryProperties(const Gedim::GeometryUtilities &geometryUtilities,
                                   const Eigen::MatrixXi &polyhedronEdges,
                                   const std::vector<Eigen::MatrixXi> &polyhedronFaces,
                                   const std::vector<bool> &polyhedronEdgeDirections,
                                   const std::vector<Eigen::MatrixXd> &mappedFaces3DVertices,
                                   VEM_PCC_3D_Inertia_Data &inertia_data) const;

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
            return utilities.ComputeDofiDofiStabilizationMatrix(localSpace.PiNabla,
                                                                localSpace.inertia_data.constantStiff,
                                                                localSpace.Dmatrix);
        case ProjectionTypes::Pi0k:
            return utilities.ComputeDofiDofiStabilizationMatrix(localSpace.Pi0k,
                                                                localSpace.inertia_data.constantMass,
                                                                localSpace.Dmatrix);
        default:
            throw std::runtime_error("not valid projection type");
        }
    }

    inline Eigen::MatrixXd ComputeBasisFunctionsValues(const VEM_PCC_3D_LocalSpace_Data &localSpace,
                                                       const ProjectionTypes &projectionType) const
    {
        return utilities.ComputeBasisFunctionsValues(projectionType,
                                                     localSpace.Nkm1,
                                                     localSpace.Pi0km1,
                                                     localSpace.Pi0k,
                                                     localSpace.VanderInternal);
    }

    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const VEM_PCC_3D_LocalSpace_Data &localSpace,
                                                                              const ProjectionTypes &projectionType) const
    {
        std::vector<Eigen::MatrixXd> basisFunctionsDerivativeValues(localSpace.Dimension);
        const Eigen::MatrixXd FmatrixInvTransp = localSpace.inertia_data.FmatrixInv.transpose();
        switch (projectionType)
        {
        case ProjectionTypes::PiNabla: {
            basisFunctionsDerivativeValues.resize(localSpace.Dimension);
            for (unsigned short i = 0; i < localSpace.Dimension; ++i)
                basisFunctionsDerivativeValues[i] = localSpace.VanderInternalDerivatives[i] * localSpace.PiNabla;
        }
        break;
        case ProjectionTypes::Pi0km1Der: {
            basisFunctionsDerivativeValues.resize(localSpace.Dimension);
            for (unsigned short i = 0; i < localSpace.Dimension; ++i)
                basisFunctionsDerivativeValues[i] = localSpace.VanderInternal.leftCols(localSpace.Nkm1) * localSpace.Pi0km1Der[i];
        }
        break;
        default:
            throw std::runtime_error("Unknown projector type");
        }

        std::vector<Eigen::MatrixXd> fmatrixInvTranspTimesBasisFunctionDerivativeValues(
            localSpace.Dimension,
            Eigen::MatrixXd::Zero(basisFunctionsDerivativeValues[0].rows(), basisFunctionsDerivativeValues[1].cols()));
        for (unsigned int d1 = 0; d1 < localSpace.Dimension; d1++)
        {
            for (unsigned int d2 = 0; d2 < localSpace.Dimension; d2++)
            {
                fmatrixInvTranspTimesBasisFunctionDerivativeValues[d1] +=
                    FmatrixInvTransp(d1, d2) * basisFunctionsDerivativeValues[d2];
            }
        }

        return fmatrixInvTranspTimesBasisFunctionDerivativeValues;
    }

    inline Eigen::MatrixXd ComputeBasisFunctionsLaplacianValues(const VEM_PCC_3D_LocalSpace_Data &, const ProjectionTypes &) const
    {
        throw std::runtime_error("Unimplemented method");
    }

    inline Eigen::MatrixXd ComputeBasisFunctionsValues(const VEM_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                       const VEM_PCC_3D_LocalSpace_Data &localSpace,
                                                       const ProjectionTypes &projectionType,
                                                       const Eigen::MatrixXd &points) const
    {
        const Eigen::MatrixXd referencePoints =
            localSpace.inertia_data.FmatrixInv * (points.colwise() - localSpace.inertia_data.translation);
        return utilities.ComputeBasisFunctionsValues(projectionType,
                                                     localSpace.Nkm1,
                                                     localSpace.Pi0km1,
                                                     localSpace.Pi0k,
                                                     ComputePolynomialsValues(reference_element_data, localSpace, referencePoints));
    }

    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const VEM_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                                              const VEM_PCC_3D_LocalSpace_Data &localSpace,
                                                                              const ProjectionTypes &projectionType,
                                                                              const Eigen::MatrixXd &points) const
    {
        const Eigen::MatrixXd referencePoints =
            localSpace.inertia_data.FmatrixInv * (points.colwise() - localSpace.inertia_data.translation);

        const Eigen::MatrixXd vander = ComputePolynomialsValues(reference_element_data, localSpace, referencePoints);

        std::vector<Eigen::MatrixXd> basisFunctionsDerivativeValues(localSpace.Dimension);
        const Eigen::MatrixXd FmatrixInvTransp = localSpace.inertia_data.FmatrixInv.transpose();
        switch (projectionType)
        {
        case ProjectionTypes::PiNabla: {
            const std::vector<Eigen::MatrixXd> VanderDerivatives =
                utilities.ComputePolynomialsDerivativeValues(reference_element_data.Monomials,
                                                             monomials,
                                                             localSpace.inertia_data.Diameter,
                                                             vander);

            basisFunctionsDerivativeValues.resize(localSpace.Dimension);
            for (unsigned short i = 0; i < localSpace.Dimension; ++i)
                basisFunctionsDerivativeValues[i] = VanderDerivatives[i] * localSpace.PiNabla;
        }
        break;
        case ProjectionTypes::Pi0km1Der: {
            basisFunctionsDerivativeValues.resize(localSpace.Dimension);
            for (unsigned short i = 0; i < localSpace.Dimension; ++i)
                basisFunctionsDerivativeValues[i] = vander.leftCols(localSpace.Nkm1) * localSpace.Pi0km1Der[i];
        }
        break;
        default:
            throw std::runtime_error("Unknown projector type");
        }

        std::vector<Eigen::MatrixXd> fmatrixInvTranspTimesBasisFunctionDerivativeValues(
            localSpace.Dimension,
            Eigen::MatrixXd::Zero(basisFunctionsDerivativeValues[0].rows(), basisFunctionsDerivativeValues[1].cols()));
        for (unsigned int d1 = 0; d1 < localSpace.Dimension; d1++)
        {
            for (unsigned int d2 = 0; d2 < localSpace.Dimension; d2++)
            {
                fmatrixInvTranspTimesBasisFunctionDerivativeValues[d1] +=
                    FmatrixInvTransp(d1, d2) * basisFunctionsDerivativeValues[d2];
            }
        }

        return fmatrixInvTranspTimesBasisFunctionDerivativeValues;
    }

    inline Eigen::MatrixXd ComputeBasisFunctionsLaplacianValues(const VEM_PCC_3D_ReferenceElement_Data &,
                                                                const VEM_PCC_3D_LocalSpace_Data &,
                                                                const ProjectionTypes &,
                                                                const Eigen::MatrixXd &) const
    {
        throw std::runtime_error("Unimplemented method");
    }

    inline Eigen::MatrixXd ComputePolynomialsValues(const VEM_PCC_3D_LocalSpace_Data &localSpace) const
    {
        return localSpace.VanderInternal;
    }

    inline Eigen::MatrixXd ComputePolynomialsValues(const VEM_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                    const VEM_PCC_3D_LocalSpace_Data &localSpace,
                                                    const Eigen::MatrixXd &points) const
    {
        const Eigen::MatrixXd referencePoints =
            localSpace.inertia_data.FmatrixInv * (points.colwise() - localSpace.inertia_data.translation);
        return utilities.ComputePolynomialsValues(reference_element_data.Monomials,
                                                  monomials,
                                                  localSpace.inertia_data.Centroid,
                                                  localSpace.inertia_data.Diameter,
                                                  referencePoints);
    }

    inline std::vector<Eigen::MatrixXd> ComputePolynomialsDerivativeValues(const VEM_PCC_3D_LocalSpace_Data &localSpace) const
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

    inline std::vector<Eigen::MatrixXd> ComputePolynomialsDerivativeValues(const VEM_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                                           const VEM_PCC_3D_LocalSpace_Data &localSpace,
                                                                           const Eigen::MatrixXd &points) const
    {
        const Eigen::MatrixXd referencePoints =
            localSpace.inertia_data.FmatrixInv * (points.colwise() - localSpace.inertia_data.translation);

        const std::vector<Eigen::MatrixXd> polynomialDerivatives = utilities.ComputePolynomialsDerivativeValues(
            reference_element_data.Monomials,
            monomials,
            localSpace.inertia_data.Diameter,
            ComputePolynomialsValues(reference_element_data, localSpace, referencePoints));

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

    inline Eigen::MatrixXd ComputePolynomialsLaplacianValues(const VEM_PCC_3D_ReferenceElement_Data &,
                                                             const VEM_PCC_3D_LocalSpace_Data &,
                                                             const Eigen::MatrixXd &) const
    {
        throw std::runtime_error("Unimplemented method");
    }

    inline Eigen::MatrixXd ComputeValuesOnEdge(const VEM_PCC_3D_LocalSpace_Data &localSpace,
                                               const Eigen::VectorXd &edgeInternalPoints,
                                               const Eigen::VectorXd &pointsCurvilinearCoordinates) const
    {
        const Eigen::VectorXd edgeBasisCoefficients = utilities.ComputeEdgeBasisCoefficients(localSpace.Order, edgeInternalPoints);
        return utilities.ComputeValuesOnEdge(edgeInternalPoints.transpose(), localSpace.Order, edgeBasisCoefficients, pointsCurvilinearCoordinates);
    }
};
} // namespace PCC
} // namespace VEM
} // namespace Polydim

#endif
