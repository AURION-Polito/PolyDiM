#ifndef __VEM_PCC_2D_LocalSpace_HPP
#define __VEM_PCC_2D_LocalSpace_HPP

#include "Eigen/Eigen"
#include "VEM_Monomials_2D.hpp"
#include "VEM_PCC_2D_LocalSpace_Data.hpp"
#include "VEM_PCC_2D_ReferenceElement.hpp"
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
///
/// Please cite the following article:
///     - <a href="https://doi.org/10.1016/j.matcom.2023.10.003">"Improving high-order VEM stability on badly-shaped
///     elements. Stefano Berrone, Gioana Teora and Fabio Vicini. (2024)"</a>

class VEM_PCC_2D_LocalSpace final
{
  private:
    VEM_PCC_Utilities<2> utilities;
    Monomials::VEM_Monomials_2D monomials;

    void InitializeProjectorsComputation(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                         const Eigen::MatrixXd &polygonVertices,
                                         const Eigen::Vector3d &polygonCentroid,
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

    inline void ComputeStabilizationMatrix(const double &polygonDiameter, VEM_PCC_2D_LocalSpace_Data &localSpace) const
    {
        localSpace.StabMatrix =
            utilities.ComputeStabilizationMatrix(localSpace.PiNabla, polygonDiameter, localSpace.Dmatrix);
    }

    inline void ComputeStabilizationMatrixPi0k(const double &polygonMeasure,
                                               VEM_PCC_2D_LocalSpace_Data &localSpace) const
    {
        localSpace.StabMatrixPi0k =
            utilities.ComputeStabilizationMatrixPi0k(localSpace.Pi0k, polygonMeasure, localSpace.Dmatrix);
    }

    void InitializeE2ProjectorsComputation(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                           const unsigned int &l,
                                           const Eigen::MatrixXd &polygonVertices,
                                           const Eigen::Vector3d &polygonCentroid,
                                           const double &polygonDiameter,
                                           const Eigen::MatrixXd &internalQuadraturePoints,
                                           const Eigen::VectorXd &internalQuadratureWeights,
                                           const Eigen::MatrixXd &internalQuadratureKLPoints,
                                           const Eigen::VectorXd &internalQuadratureKLWeights,
                                           const Eigen::MatrixXd &boundaryQuadraturePoints,
                                           VEM_PCC_2D_LocalSpace_Data &localSpace) const;

    void ComputeL2ProjectorsKL(VEM_PCC_2D_LocalSpace_Data &localSpace) const;

  public:
    /// \brief Create and Initialize all the variables contained in \ref VEM::PCC::VEM_PCC_2D_LocalSpace_Data
    /// \param reference_element_data: an object of type \ref VEM::PCC::VEM_PCC_2D_ReferenceElement_Data which contains
    /// monomials, quadrature and the number of degrees of freedom, counting in order DOFS associated with vertices,
    /// edges and internal values.
    /// \param polygon: an object of type \ref VEM::PCC::VEM_PCC_2D_Polygon_Geometry which
    /// contains the geoemtric properties of the elements.
    /// \return An object of type \ref VEM::PCC::VEM_PCC_2D_LocalSpace_Data.
    VEM_PCC_2D_LocalSpace_Data CreateLocalSpace(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                const VEM_PCC_2D_Polygon_Geometry &polygon) const;

    /// \brief Create and Initialize variables contained in \ref VEM::PCC::VEM_PCC_2D_LocalSpace_Data which are used by
    /// \ref VEM::PCC::VEM_PCC_3D_LocalSpace. \param reference_element_data: an object of type \ref
    /// VEM::PCC::VEM_PCC_2D_ReferenceElement_Data which contains monomials, quadrature and the number of degrees of
    /// freedom, counting in order DOFS associated with vertices, edges and internal values. \param polygon: an object
    /// of type \ref VEM::PCC::VEM_PCC_2D_Polygon_Geometry which contains the geoemtric properties of the elements.
    /// \return An object of type \ref VEM::PCC::VEM_PCC_2D_LocalSpace_Data.
    VEM_PCC_2D_LocalSpace_Data Compute3DUtilities(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                  const VEM_PCC_2D_Polygon_Geometry &polygon) const;

    VEM_PCC_2D_LocalSpace_Data Compute3DUtilities_DF_PCC(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                         const VEM_PCC_2D_Polygon_Geometry &polygon) const;

    /// \brief Compute the values of projections of VEM basis functions at the internal quadrature points.
    /// \param localSpace: an object of type \ref VEM::PCC::VEM_PCC_2D_LocalSpace_Data which contains local matrices.
    /// \param projectionType: the \ref VEM::PCC::ProjectionTypes reporting the kind of projector used to access to the
    /// point-wise evalution of VE basis functions. \return A matrix of size numQuadrature \f$\times\f$ numDOFs whose
    /// columns contain the evaluation of the projection of each basis function at the internal quadrature points.
    inline Eigen::MatrixXd ComputeBasisFunctionsValues(const VEM_PCC_2D_LocalSpace_Data &localSpace,
                                                       const ProjectionTypes &projectionType) const
    {
        if (projectionType == ProjectionTypes::Pi0klm1)
            return localSpace.VanderInternalKL * localSpace.Pi0klm1;

        return utilities.ComputeBasisFunctionsValues(
            projectionType, localSpace.Nkm1, localSpace.Pi0km1, localSpace.Pi0k, localSpace.VanderInternal);
    }

    /// \brief Compute the values of projections of VEM basis function derivatives at the internal quadrature points.
    /// \param localSpace: an object of type \ref VEM::PCC::VEM_PCC_2D_LocalSpace_Data which contains local matrices.
    /// \param projectionType: the \ref VEM::PCC::ProjectionTypes reporting the kind of projector used to access to the
    /// point-wise evalution of VE basis function derivatives. \return A vector of 2 matrices of size numQuadrature
    /// \f$\times\f$ numDOFs whose columns contain the evaluation of the projection of each basis function derivatives
    /// with respect x and y, respectively, at the internal quadrature points.
    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(
        const VEM_PCC_2D_LocalSpace_Data &localSpace, const ProjectionTypes &projectionType) const
    {
        return utilities.ComputeBasisFunctionsDerivativeValues(projectionType,
                                                               localSpace.Nkm1,
                                                               localSpace.VanderInternal,
                                                               localSpace.VanderInternalDerivatives,
                                                               localSpace.PiNabla,
                                                               localSpace.Pi0km1Der);
    }

    /// \brief Compute the values of VEM basis function laplacian at the internal quadrature points, here approximated
    /// using \ref VEM::PCC::ProjectionTypes::Pi0km1Der. \param localSpace: an object of type \ref
    /// VEM::PCC::VEM_PCC_2D_LocalSpace_Data which contains local matrices. \return A matrix of size numQuadrature
    /// \f$\times\f$ numDOFs whose columns contain the evaluation of the approximated laplacian at the internal
    /// quadrature points.
    inline Eigen::MatrixXd ComputeBasisFunctionsLaplacianValues(const VEM_PCC_2D_LocalSpace_Data &localSpace) const
    {
        return utilities.ComputeBasisFunctionsLaplacianValues(
            localSpace.Nkm1, localSpace.VanderInternalDerivatives, localSpace.Pi0km1Der);
    }

    /// \brief Compute the values of projections of VEM basis functions at points.
    /// \param reference_element_data: an object of type \ref VEM::PCC::VEM_PCC_2D_ReferenceElement_Data which contains
    /// monomials stuff. \param localSpace: an object of type \ref VEM::PCC::VEM_PCC_2D_LocalSpace_Data which contains
    /// local matrices. \param projectionType: the \ref VEM::PCC::ProjectionTypes reporting the kind of projector used
    /// to access to the point-wise evalution of VE basis functions. \param points: a matrix 3 \f$\times\f$ numPoints
    /// reporting the coordinates of points. \return A matrix of size numPoints \f$\times\f$ numDOFs whose columns
    /// contain the evaluation of the projection of each basis function at points.
    inline Eigen::MatrixXd ComputeBasisFunctionsValues(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                       const VEM_PCC_2D_LocalSpace_Data &localSpace,
                                                       const ProjectionTypes &projectionType,
                                                       const Eigen::MatrixXd &points) const
    {
        return utilities.ComputeBasisFunctionsValues(
            projectionType,
            localSpace.Nkm1,
            localSpace.Pi0km1,
            localSpace.Pi0k,
            ComputePolynomialsValues(reference_element_data, localSpace, points));
    }

    /// \brief Compute the values of projections of VEM basis function derivatives at points.
    /// \param reference_element_data: an object of type \ref VEM::PCC::VEM_PCC_2D_ReferenceElement_Data which contains
    /// monomials stuff. \param localSpace: an object of type \ref VEM::PCC::VEM_PCC_2D_LocalSpace_Data which contains
    /// local matrices. \param projectionType: the \ref VEM::PCC::ProjectionTypes reporting the kind of projector used
    /// to access to the point-wise evalution of VE basis function derivatives. \param points: a matrix 3 \f$\times\f$
    /// numPoints reporting the coordinates of points. \return A vector of 2 matrices of size numPoints \f$\times\f$
    /// numDOFs whose columns contain the evaluation of the projection of each basis function derivatives with respect x
    /// and y, respectively, at points.
    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(
        const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
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

    /// \brief Compute the values of VEM basis function laplacian at points, here approximated using \ref
    /// VEM::PCC::ProjectionTypes::Pi0km1Der. \param reference_element_data: an object of type \ref
    /// VEM::PCC::VEM_PCC_2D_ReferenceElement_Data which contains monomials stuff. \param localSpace: an object of type
    /// \ref VEM::PCC::VEM_PCC_2D_LocalSpace_Data which contains local matrices. \param points: a matrix 3 \f$\times\f$
    /// numPoints reporting the coordinates of points. \return A matrix of size numPoints \f$\times\f$ numDOFs whose
    /// columns contain the evaluation of the approximated laplacian at points.
    inline Eigen::MatrixXd ComputeBasisFunctionsLaplacianValues(
        const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
        const VEM_PCC_2D_LocalSpace_Data &localSpace,
        const Eigen::MatrixXd &points) const
    {
        return utilities.ComputeBasisFunctionsLaplacianValues(
            localSpace.Nkm1,
            localSpace.Pi0km1Der,
            ComputePolynomialsDerivativeValues(reference_element_data, localSpace, points));
    }

    /// \brief Compute the values of monomial basis functions at the internal quadrature points.
    /// \param localSpace: an object of type \ref VEM::PCC::VEM_PCC_2D_LocalSpace_Data which contains local matrices.
    /// \return A matrix of size numQuadrature \f$\times\f$ numMonomials whose columns contain the evaluation of
    /// monomials at the internal quadrature points.
    inline Eigen::MatrixXd ComputePolynomialsValues(const VEM_PCC_2D_LocalSpace_Data &localSpace) const
    {
        return utilities.ComputePolynomialsValues(localSpace.VanderInternal);
    }

    /// \brief Compute the values of monomial basis functions at points.
    /// \param reference_element_data: an object of type \ref VEM::PCC::VEM_PCC_2D_ReferenceElement_Data which contains
    /// monomials stuff. \param points: a matrix 3 \f$\times\f$ numPoints reporting the coordinates of points. \return A
    /// matrix of size numPoints \f$\times\f$ numMonomials whose columns contain the evaluation of monomials at points.
    inline Eigen::MatrixXd ComputePolynomialsValues(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                    const VEM_PCC_2D_LocalSpace_Data &localSpace,
                                                    const Eigen::MatrixXd &points) const
    {
        return utilities.ComputePolynomialsValues(
            reference_element_data.Monomials, monomials, localSpace.Centroid, localSpace.Diameter, points);
    }

    /// \brief Compute the values of monomial basis function derivatives at the internal quadrature points.
    /// \param localSpace: an object of type \ref VEM::PCC::VEM_PCC_2D_LocalSpace_Data which contains local matrices.
    /// \return A vector of two matrices of size numQuadrature \f$\times\f$ numMonomials whose columns contain the
    /// evaluation of monomials derivatives with respect x and y, respectively, at the internal quadrature points
    inline std::vector<Eigen::MatrixXd> ComputePolynomialsDerivativeValues(
        const VEM_PCC_2D_LocalSpace_Data &localSpace) const
    {
        return utilities.ComputePolynomialsDerivativeValues(localSpace.VanderInternalDerivatives);
    }

    /// \brief Compute the values of monomial basis functions at points.
    /// \param reference_element_data: an object of type \ref VEM::PCC::VEM_PCC_2D_ReferenceElement_Data which contains
    /// monomials stuff. \param points: a matrix 3 \f$\times\f$ numPoints reporting the coordinates of points. \return A
    /// vector of two matrices of size numPoints \f$\times\f$ numMonomials whose columns contain the evaluation of
    /// monomials derivatives with respect x and y, respectively, at points
    inline std::vector<Eigen::MatrixXd> ComputePolynomialsDerivativeValues(
        const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
        const VEM_PCC_2D_LocalSpace_Data &localSpace,
        const Eigen::MatrixXd &points) const
    {
        return utilities.ComputePolynomialsDerivativeValues(
            reference_element_data.Monomials,
            monomials,
            localSpace.Diameter,
            ComputePolynomialsValues(reference_element_data, localSpace, points));
    }

    /// \brief Compute the values of monomials laplacian at points.
    /// \param reference_element_data: an object of type \ref VEM::PCC::VEM_PCC_2D_ReferenceElement_Data which contains
    /// monomials stuff. \param points: a matrix 3 \f$\times\f$ numPoints reporting the coordinates of points. \return A
    /// matrix of size numPoints \f$\times\f$ numMonomials whose columns contain the evaluation of monomials at points.
    inline Eigen::MatrixXd ComputePolynomialsLaplacianValues(
        const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
        const VEM_PCC_2D_LocalSpace_Data &localSpace,
        const Eigen::MatrixXd &points) const
    {
        return utilities.ComputePolynomialsLaplacianValues(
            reference_element_data.Monomials,
            monomials,
            localSpace.Diameter,
            ComputePolynomialsValues(reference_element_data, localSpace, points));
    }

    inline Eigen::MatrixXd ComputeValuesOnEdge(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                               const Eigen::VectorXd &pointsCurvilinearCoordinates) const
    {
        Eigen::RowVectorXd edgeInternalPoints;
        if (reference_element_data.Quadrature.ReferenceSegmentInternalPoints.rows() > 0)
            edgeInternalPoints = reference_element_data.Quadrature.ReferenceSegmentInternalPoints.row(0);
        const Eigen::VectorXd edgeBasisCoefficients =
            utilities.ComputeEdgeBasisCoefficients(reference_element_data.Order, edgeInternalPoints);

        return utilities.ComputeValuesOnEdge(
            edgeInternalPoints, reference_element_data.Order, edgeBasisCoefficients, pointsCurvilinearCoordinates);
    }
};
} // namespace PCC
} // namespace VEM
} // namespace Polydim

#endif
