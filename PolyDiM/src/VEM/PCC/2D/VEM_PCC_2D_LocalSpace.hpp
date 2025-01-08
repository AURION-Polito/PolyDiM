#ifndef __VEM_PCC_2D_LocalSpace_HPP
#define __VEM_PCC_2D_LocalSpace_HPP

#include "I_VEM_PCC_2D_LocalSpace.hpp"
#include "VEM_Monomials_2D.hpp"

namespace Polydim
{
namespace VEM
{
namespace PCC
{

/// \brief Class used for computing values of basis functions of 2D
/// Primal Conforming Constant degree Virtual Element Methods.
///
/// Please cite the following article:
///     - <a href="https://doi.org/10.1016/j.matcom.2023.10.003">"Improving high-order VEM stability on badly-shaped
///     elements. Stefano Berrone, Gioana Teora and Fabio Vicini. (2024)"</a>

class VEM_PCC_2D_LocalSpace final : public I_VEM_PCC_2D_LocalSpace
{
  private:
    VEM_PCC_Utilities<2> utilities;
    Monomials::VEM_Monomials_2D monomials;

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

    void ComputePolynomialsDofs(const double &polytopeMeasure, VEM_PCC_2D_LocalSpace_Data &localSpace) const;

    void InitializeE2ProjectorsComputation(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                           const unsigned int &l,
                                           const Eigen::MatrixXd &polygonVertices,
                                           const Eigen::Vector3d &polygonCentroid,
                                           const double &polygonMeasure,
                                           const double &polygonDiameter,
                                           const Eigen::MatrixXd &internalQuadraturePoints,
                                           const Eigen::VectorXd &internalQuadratureWeights,
                                           const Eigen::MatrixXd &internalQuadratureKLPoints,
                                           const Eigen::VectorXd &internalQuadratureKLWeights,
                                           const Eigen::MatrixXd &boundaryQuadraturePoints,
                                           VEM_PCC_2D_LocalSpace_Data &localSpace) const;

    void ComputeL2ProjectorsKL(VEM_PCC_2D_LocalSpace_Data &localSpace) const;

  public:
    VEM_PCC_2D_LocalSpace_Data CreateLocalSpace(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
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

    VEM_PCC_2D_LocalSpace_Data Compute3DUtilities(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                  const VEM_PCC_2D_Polygon_Geometry &polygon) const;

    VEM_PCC_2D_LocalSpace_Data Compute3DUtilities_DF_PCC(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                         const VEM_PCC_2D_Polygon_Geometry &polygon) const;

    inline Eigen::MatrixXd ComputeBasisFunctionsValues(const VEM_PCC_2D_LocalSpace_Data &localSpace,
                                                       const ProjectionTypes &projectionType) const
    {
        if (projectionType == ProjectionTypes::Pi0klm1)
            return localSpace.VanderInternalKL * localSpace.Pi0klm1;

        return utilities.ComputeBasisFunctionsValues(projectionType,
                                                     localSpace.Nkm1,
                                                     localSpace.Pi0km1,
                                                     localSpace.Pi0k,
                                                     localSpace.VanderInternal);
    }

    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const VEM_PCC_2D_LocalSpace_Data &localSpace,
                                                                              const ProjectionTypes &projectionType) const
    {
        return utilities.ComputeBasisFunctionsDerivativeValues(projectionType,
                                                               localSpace.Nkm1,
                                                               localSpace.VanderInternal,
                                                               localSpace.VanderInternalDerivatives,
                                                               localSpace.PiNabla,
                                                               localSpace.Pi0km1Der);
    }

    inline Eigen::MatrixXd ComputeBasisFunctionsValues(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                       const VEM_PCC_2D_LocalSpace_Data &localSpace,
                                                       const ProjectionTypes &projectionType,
                                                       const Eigen::MatrixXd &points) const
    {
        return utilities.ComputeBasisFunctionsValues(projectionType,
                                                     localSpace.Nkm1,
                                                     localSpace.Pi0km1,
                                                     localSpace.Pi0k,
                                                     ComputePolynomialsValues(reference_element_data, localSpace, points));
    }

    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                                              const VEM_PCC_2D_LocalSpace_Data &localSpace,
                                                                              const ProjectionTypes &projectionType,
                                                                              const Eigen::MatrixXd &points) const
    {
        return utilities.ComputeBasisFunctionsDerivativeValues(
            projectionType,
            localSpace.Nkm1,
            ComputePolynomialsValues(reference_element_data, localSpace, points),
            ComputePolynomialsDerivativeValues(reference_element_data, localSpace, points),
            localSpace.PiNabla,
            localSpace.Pi0km1Der);
    }

    inline Eigen::MatrixXd ComputeBasisFunctionsLaplacianValues(const VEM_PCC_2D_LocalSpace_Data &localSpace,
                                                                const ProjectionTypes &projectionType) const
    {
        return utilities.ComputeBasisFunctionsLaplacianValues(projectionType,
                                                              localSpace.Nkm1,
                                                              localSpace.VanderInternalDerivatives,
                                                              localSpace.Pi0km1Der);
    }

    inline Eigen::MatrixXd ComputeBasisFunctionsLaplacianValues(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                                const VEM_PCC_2D_LocalSpace_Data &localSpace,
                                                                const ProjectionTypes &projectionType,
                                                                const Eigen::MatrixXd &points) const
    {
        return utilities.ComputeBasisFunctionsLaplacianValues(
            projectionType,
            localSpace.Nkm1,
            localSpace.Pi0km1Der,
            ComputePolynomialsDerivativeValues(reference_element_data, localSpace, points));
    }

    inline Eigen::MatrixXd ComputePolynomialsValues(const VEM_PCC_2D_LocalSpace_Data &localSpace) const
    {
        return utilities.ComputePolynomialsValues(localSpace.VanderInternal);
    }

    inline Eigen::MatrixXd ComputePolynomialsValues(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                    const VEM_PCC_2D_LocalSpace_Data &localSpace,
                                                    const Eigen::MatrixXd &points) const
    {
        return utilities.ComputePolynomialsValues(reference_element_data.Monomials,
                                                  monomials,
                                                  localSpace.Centroid,
                                                  localSpace.Diameter,
                                                  points);
    }

    inline std::vector<Eigen::MatrixXd> ComputePolynomialsDerivativeValues(const VEM_PCC_2D_LocalSpace_Data &localSpace) const
    {
        return utilities.ComputePolynomialsDerivativeValues(localSpace.VanderInternalDerivatives);
    }
    inline std::vector<Eigen::MatrixXd> ComputePolynomialsDerivativeValues(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                                           const VEM_PCC_2D_LocalSpace_Data &localSpace,
                                                                           const Eigen::MatrixXd &points) const
    {
        return utilities.ComputePolynomialsDerivativeValues(reference_element_data.Monomials,
                                                            monomials,
                                                            localSpace.Diameter,
                                                            ComputePolynomialsValues(reference_element_data, localSpace, points));
    }
    inline Eigen::MatrixXd ComputePolynomialsLaplacianValues(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                             const VEM_PCC_2D_LocalSpace_Data &localSpace,
                                                             const Eigen::MatrixXd &points) const
    {
        return utilities.ComputePolynomialsLaplacianValues(reference_element_data.Monomials,
                                                           monomials,
                                                           localSpace.Diameter,
                                                           ComputePolynomialsValues(reference_element_data, localSpace, points));
    }

    inline Eigen::MatrixXd ComputeValuesOnEdge(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                               const Eigen::VectorXd &pointsCurvilinearCoordinates) const
    {
        Eigen::RowVectorXd edgeInternalPoints;
        if (reference_element_data.Quadrature.ReferenceEdgeDOFsInternalPoints.rows() > 0)
            edgeInternalPoints = reference_element_data.Quadrature.ReferenceEdgeDOFsInternalPoints.row(0);
        const Eigen::VectorXd edgeBasisCoefficients =
            utilities.ComputeEdgeBasisCoefficients(reference_element_data.Order, edgeInternalPoints);

        return utilities.ComputeValuesOnEdge(edgeInternalPoints, reference_element_data.Order, edgeBasisCoefficients, pointsCurvilinearCoordinates);
    }
};
} // namespace PCC
} // namespace VEM
} // namespace Polydim

#endif
