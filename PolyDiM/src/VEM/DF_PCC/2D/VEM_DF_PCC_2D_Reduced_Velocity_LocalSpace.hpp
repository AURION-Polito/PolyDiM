#ifndef __VEM_DF_PCC_2D_Reduced_Velocity_LocalSpace_HPP
#define __VEM_DF_PCC_2D_Reduced_Velocity_LocalSpace_HPP

#include "Eigen/Eigen"
#include "I_VEM_DF_PCC_2D_ReferenceElement.hpp"
#include "I_VEM_DF_PCC_2D_Velocity_LocalSpace.hpp"
#include "VEM_DF_PCC_2D_LocalSpace_Data.hpp"
#include "VEM_DF_PCC_Utilities.hpp"
#include "VEM_Monomials_2D.hpp"
#include <vector>

namespace Polydim
{
namespace VEM
{
namespace DF_PCC
{

/// \brief Class used for computing values of basis functions of 2D
/// Divergence Free Primal Conforming Constant degree Virtual Element Methods.
///
/// Please cite the following article:
///     - <a href="https://doi.org/10.1016/j.matcom.2023.10.003">"Improving high-order VEM stability on badly-shaped
///     elements. Stefano Berrone, Gioana Teora and Fabio Vicini. (2024)"</a>

class VEM_DF_PCC_2D_Reduced_Velocity_LocalSpace final : public I_VEM_DF_PCC_2D_Velocity_LocalSpace
{
  private:
    VEM_DF_PCC_Utilities<2> utilities;
    Monomials::VEM_Monomials_2D monomials;
    Monomials::VEM_GBasis_2D g_basis;

    void InitializeProjectorsComputation(const VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                         const Eigen::MatrixXd &polygonVertices,
                                         const Eigen::Vector3d &polygonCentroid,
                                         const double &polygonMeasure,
                                         const double &polygonDiameter,
                                         const std::vector<bool> &edgeDirections,
                                         const Eigen::MatrixXd &internalQuadraturePoints,
                                         const Eigen::VectorXd &internalQuadratureWeights,
                                         const Eigen::MatrixXd &boundaryQuadraturePoints,
                                         const Eigen::MatrixXd &boundaryDofQuadraturePoints,
                                         const Eigen::MatrixXd &referenceEdgeInternalPoints,
                                         const Eigen::MatrixXd &referenceEdgeDofInternalPoints,
                                         VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace) const;

    void ComputeDivergenceCoefficients(const double &polygonMeasure,
                                       const std::vector<Eigen::VectorXd> &boundaryDofQuadratureWeightsTimesNormal,
                                       VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace) const;

    void ComputePolynomialBasisDofs(const Eigen::VectorXd &internalQuadratureWeights,
                                    VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace) const;

    void ComputeCMatrixkm2(const double &polygonDiameter,
                           const std::vector<Eigen::VectorXd> &boundaryQuadratureWeightsTimesNormal,
                           VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace) const;

    void ComputePiNabla(const VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                        const double &polygonMeasure,
                        const double &polygonDiameter,
                        const Eigen::VectorXd &internalQuadratureWeights,
                        const std::vector<Eigen::VectorXd> &boundaryDofQuadratureWeightsTimesNormal,
                        VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace) const;

    void ComputeL2Projectors(const Eigen::VectorXd &internalQuadratureWeights, VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace) const;

    void ComputeL2ProjectorsOfDerivatives(const VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                          const double &polygonDiameter,
                                          const std::vector<Eigen::VectorXd> &boundaryDofQuadratureWeightsTimesNormal,
                                          VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace) const;

  public:
    virtual ~VEM_DF_PCC_2D_Reduced_Velocity_LocalSpace()
    {
    }

    VEM_DF_PCC_2D_Velocity_LocalSpace_Data CreateLocalSpace(const VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                                            const VEM_DF_PCC_2D_Polygon_Geometry &polygon) const;

    inline Eigen::MatrixXd ComputeDofiDofiStabilizationMatrix(const VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace,
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

    inline Eigen::MatrixXd ComputeValuesOnEdge(const VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                               const Eigen::VectorXd &pointsCurvilinearCoordinates) const
    {
        Eigen::RowVectorXd edgeInternalPoints;
        if (reference_element_data.Quadrature.ReferenceEdgeDOFsInternalPoints.rows() > 0)
            edgeInternalPoints = reference_element_data.Quadrature.ReferenceEdgeDOFsInternalPoints.row(0);
        const Eigen::VectorXd edgeBasisCoefficients =
            utilities.ComputeEdgeBasisCoefficients(reference_element_data.Order, edgeInternalPoints);

        return utilities.ComputeValuesOnEdge(edgeInternalPoints, reference_element_data.Order, edgeBasisCoefficients, pointsCurvilinearCoordinates);
    }

    /// \brief Compute the values of projections of VEM basis functions at the internal quadrature points.
    /// \param localSpace: an object of type \ref VEM::PCC::VEM_DF_PCC_2D_Velocity_LocalSpace_Data which contains local
    /// matrices. \param projectionType: the \ref VEM::PCC::ProjectionTypes reporting the kind of projector used to
    /// access to the point-wise evalution of VE basis functions. \return A matrix of size numQuadrature \f$\times\f$
    /// numDOFs whose columns contain the evaluation of the projection of each basis function at the internal quadrature
    /// points.
    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsValues(const VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace,
                                                                    const ProjectionTypes &projectionType) const
    {
        return utilities.ComputeBasisFunctionsValues(projectionType,
                                                     localSpace.Nkm2,
                                                     localSpace.Pi0km2,
                                                     localSpace.Pi0k,
                                                     localSpace.VanderInternal);
    }

    /// \brief Compute the values of projections of VEM basis function derivatives at the internal quadrature points.
    /// \param localSpace: an object of type \ref VEM::PCC::VEM_DF_PCC_2D_Velocity_LocalSpace_Data which contains local
    /// matrices. \param projectionType: the \ref VEM::PCC::ProjectionTypes reporting the kind of projector used to
    /// access to the point-wise evalution of VE basis function derivatives. \return A vector of 2 matrices of size
    /// numQuadrature \f$\times\f$ numDOFs whose columns contain the evaluation of the projection of each basis function
    /// derivatives with respect x and y, respectively, at the internal quadrature points.
    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace,
                                                                              const ProjectionTypes &projectionType) const
    {
        return utilities.ComputeBasisFunctionsDerivativeValues(projectionType,
                                                               localSpace.Nkm1,
                                                               localSpace.VanderInternal,
                                                               localSpace.VanderInternalDerivatives,
                                                               localSpace.PiNabla,
                                                               localSpace.Pi0km1Der);
    }

    /// \brief Compute the values of projections of VEM basis functions at points.
    /// \param reference_element_data: an object of type \ref VEM::PCC::VEM_DF_PCC_2D_Velocity_ReferenceElement_Data
    /// which contains monomials stuff. \param polygon: an object of type \ref VEM::PCC::VEM_DF_PCC_2D_Polygon_Geometry
    /// which contains the geoemtric properties of the elements. \param localSpace: an object of type \ref
    /// VEM::PCC::VEM_DF_PCC_2D_Velocity_LocalSpace_Data which contains local matrices. \param projectionType: the \ref
    /// VEM::PCC::ProjectionTypes reporting the kind of projector used to access to the point-wise evalution of VE basis
    /// functions. \param points: a matrix 3 \f$\times\f$ numPoints reporting the coordinates of points. \return A
    /// matrix of size numPoints \f$\times\f$ numDOFs whose columns contain the evaluation of the projection of each
    /// basis function at points.
    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsValues(const VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                                                    const VEM_DF_PCC_2D_Polygon_Geometry &polygon,
                                                                    const VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace,
                                                                    const ProjectionTypes &projectionType,
                                                                    const Eigen::MatrixXd &points) const
    {
        return utilities.ComputeBasisFunctionsValues(projectionType,
                                                     localSpace.Nkm2,
                                                     localSpace.Pi0km2,
                                                     localSpace.Pi0k,
                                                     ComputePolynomialsValues(reference_element_data, polygon, points));
    }

    /// \brief Compute the values of projections of VEM basis function derivatives at points.
    /// \param reference_element_data: an object of type \ref VEM::PCC::VEM_DF_PCC_2D_Velocity_ReferenceElement_Data
    /// which contains monomials stuff. \param polygon: an object of type \ref VEM::PCC::VEM_DF_PCC_2D_Polygon_Geometry
    /// which contains the geoemtric properties of the elements. \param localSpace: an object of type \ref
    /// VEM::PCC::VEM_DF_PCC_2D_Velocity_LocalSpace_Data which contains local matrices. \param projectionType: the \ref
    /// VEM::PCC::ProjectionTypes reporting the kind of projector used to access to the point-wise evalution of VE basis
    /// function derivatives. \param points: a matrix 3 \f$\times\f$ numPoints reporting the coordinates of points.
    /// \return A vector of 2 matrices of size numPoints \f$\times\f$ numDOFs whose columns contain the evaluation of
    /// the projection of each basis function derivatives with respect x and y, respectively, at points.
    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                                                              const VEM_DF_PCC_2D_Polygon_Geometry &polygon,
                                                                              const VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace,
                                                                              const ProjectionTypes &projectionType,
                                                                              const Eigen::MatrixXd &points) const
    {
        return utilities.ComputeBasisFunctionsDerivativeValues(
            projectionType,
            localSpace.Nkm1,
            ComputePolynomialsValues(reference_element_data, polygon, points),
            ComputePolynomialsDerivativeValues(reference_element_data, polygon, points),
            localSpace.PiNabla,
            localSpace.Pi0km1Der);
    }

    /// \brief Compute the values of monomial basis functions at the internal quadrature points.
    /// \param localSpace: an object of type \ref VEM::PCC::VEM_DF_PCC_2D_Velocity_LocalSpace_Data which contains local
    /// matrices. \return A matrix of size numQuadrature \f$\times\f$ numMonomials whose columns contain the evaluation
    /// of monomials at the internal quadrature points.
    inline Eigen::MatrixXd ComputePolynomialsValues(const VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace) const
    {
        return utilities.ComputePolynomialsValues(localSpace.VanderInternal);
    }

    /// \brief Compute the values of monomial basis functions at points.
    /// \param reference_element_data: an object of type \ref VEM::PCC::VEM_DF_PCC_2D_Velocity_ReferenceElement_Data
    /// which contains monomials stuff. \param polygon: an object of type \ref VEM::PCC::VEM_DF_PCC_2D_Polygon_Geometry
    /// which contains the geoemtric properties of the elements. \param points: a matrix 3 \f$\times\f$ numPoints
    /// reporting the coordinates of points. \return A matrix of size numPoints \f$\times\f$ numMonomials whose columns
    /// contain the evaluation of monomials at points.
    inline Eigen::MatrixXd ComputePolynomialsValues(const VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                                    const VEM_DF_PCC_2D_Polygon_Geometry &polygon,
                                                    const Eigen::MatrixXd &points) const
    {
        return utilities.ComputePolynomialsValues(reference_element_data.Monomials, monomials, polygon.Centroid, polygon.Diameter, points);
    }

    /// \brief Compute the values of monomial basis function derivatives at the internal quadrature points.
    /// \param localSpace: an object of type \ref VEM::PCC::VEM_DF_PCC_2D_Velocity_LocalSpace_Data which contains local
    /// matrices. \return A vector of two matrices of size numQuadrature \f$\times\f$ numMonomials whose columns contain
    /// the evaluation of monomials derivatives with respect x and y, respectively, at the internal quadrature points
    inline std::vector<Eigen::MatrixXd> ComputePolynomialsDerivativeValues(const VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace) const
    {
        return utilities.ComputePolynomialsDerivativeValues(localSpace.VanderInternalDerivatives);
    }

    /// \brief Compute the values of monomial basis functions at points.
    /// \param reference_element_data: an object of type \ref VEM::PCC::VEM_DF_PCC_2D_Velocity_ReferenceElement_Data
    /// which contains monomials stuff. \param polygon: an object of type \ref VEM::PCC::VEM_DF_PCC_2D_Polygon_Geometry
    /// which contains the geoemtric properties of the elements. \param points: a matrix 3 \f$\times\f$ numPoints
    /// reporting the coordinates of points. \return A vector of two matrices of size numPoints \f$\times\f$
    /// numMonomials whose columns contain the evaluation of monomials derivatives with respect x and y, respectively,
    /// at points
    inline std::vector<Eigen::MatrixXd> ComputePolynomialsDerivativeValues(const VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                                                           const VEM_DF_PCC_2D_Polygon_Geometry &polygon,
                                                                           const Eigen::MatrixXd &points) const
    {
        return utilities.ComputePolynomialsDerivativeValues(reference_element_data.Monomials,
                                                            monomials,
                                                            polygon.Diameter,
                                                            ComputePolynomialsValues(reference_element_data, polygon, points));
    }

    /// \brief Compute the values of monomials laplacian at points.
    /// \param reference_element_data: an object of type \ref VEM::PCC::VEM_DF_PCC_2D_Velocity_ReferenceElement_Data
    /// which contains monomials stuff. \param polygon: an object of type \ref VEM::PCC::VEM_DF_PCC_2D_Polygon_Geometry
    /// which contains the geoemtric properties of the elements. \param points: a matrix 3 \f$\times\f$ numPoints
    /// reporting the coordinates of points. \return A matrix of size numPoints \f$\times\f$ numMonomials whose columns
    /// contain the evaluation of monomials at points.
    inline Eigen::MatrixXd ComputePolynomialsLaplacianValues(const VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                                             const VEM_DF_PCC_2D_Polygon_Geometry &polygon,
                                                             const Eigen::MatrixXd &points) const
    {
        return utilities.ComputePolynomialsLaplacianValues(reference_element_data.Monomials,
                                                           monomials,
                                                           polygon.Diameter,
                                                           ComputePolynomialsValues(reference_element_data, polygon, points));
    }

    inline Eigen::MatrixXd ComputeBasisFunctionsDivergenceValues(const VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace) const
    {
        return utilities.ComputeBasisFunctionsDivergenceValues(1, localSpace.VanderInternal, localSpace.Vmatrix);
    }
};
} // namespace DF_PCC
} // namespace VEM
} // namespace Polydim

#endif
