#ifndef __VEM_PCC_2D_Inertia_LocalSpace_HPP
#define __VEM_PCC_2D_Inertia_LocalSpace_HPP

#include "GeometryUtilities.hpp"
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

class VEM_PCC_2D_Inertia_LocalSpace final : public I_VEM_PCC_2D_LocalSpace
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

    void InertiaMapping(const VEM_PCC_2D_Polygon_Geometry &polygon, VEM_PCC_2D_Inertia_Data &inertia_data) const;

    void ComputeGeometryProperties(const Gedim::GeometryUtilities &geometryUtilities,
                                   const std::vector<bool> &polygonEdgeDirections,
                                   const std::vector<Eigen::Matrix3d> &polygonTriangulation,
                                   VEM_PCC_2D_Inertia_Data &data) const;

public:
    /// \brief Create and Initialize all the variables contained in \ref VEM::PCC::VEM_PCC_2D_LocalSpace_Data
    /// \param reference_element_data: an object of type \ref VEM::PCC::VEM_PCC_2D_ReferenceElement_Data which contains
    /// monomials, quadrature and the number of degrees of freedom, counting in order DOFS associated with vertices,
    /// edges and internal values. \param polygon: an object of type \ref VEM::PCC::VEM_PCC_2D_Polygon_Geometry which
    /// contains the geoemtric properties of the elements. \return An object of type \ref
    /// VEM::PCC::VEM_PCC_2D_LocalSpace_Data.
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

    /// \brief Compute the values of projections of VEM basis functions at the internal quadrature points.
    /// \param localSpace: an object of type \ref VEM::PCC::VEM_PCC_2D_LocalSpace_Data which contains local matrices.
    /// \param projectionType: the \ref VEM::PCC::ProjectionTypes reporting the kind of projector used to access to the
    /// point-wise evalution of VE basis functions. \return A matrix of size numQuadrature \f$\times\f$ numDOFs whose
    /// columns contain the evaluation of the projection of each basis function at the internal quadrature points.
    inline Eigen::MatrixXd ComputeBasisFunctionsValues(const VEM_PCC_2D_LocalSpace_Data &localSpace,
                                                       const ProjectionTypes &projectionType) const
    {
        return utilities.ComputeBasisFunctionsValues(projectionType,
                                                     localSpace.Nkm1,
                                                     localSpace.Pi0km1,
                                                     localSpace.Pi0k,
                                                     localSpace.VanderInternal);
    }

    /// \brief Compute the values of projections of VEM basis function derivatives at the internal quadrature points.
    /// \param localSpace: an object of type \ref VEM::PCC::VEM_PCC_2D_LocalSpace_Data which contains local matrices.
    /// \param projectionType: the \ref VEM::PCC::ProjectionTypes reporting the kind of projector used to access to the
    /// point-wise evalution of VE basis function derivatives. \return A vector of 2 matrices of size numQuadrature
    /// \f$\times\f$ numDOFs whose columns contain the evaluation of the projection of each basis function derivatives
    /// with respect x and y, respectively, at the internal quadrature points.
    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const VEM_PCC_2D_LocalSpace_Data &localSpace,
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
                basisFunctionsDerivativeValues[i] =
                    localSpace.VanderInternal.leftCols(localSpace.Nkm1) * localSpace.Pi0km1Der[i];
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
        const Eigen::MatrixXd referencePoints = localSpace.inertia_data.FmatrixInv
                                                * (points.colwise() - localSpace.inertia_data.translation);
        return utilities.ComputeBasisFunctionsValues(projectionType,
                                                     localSpace.Nkm1,
                                                     localSpace.Pi0km1,
                                                     localSpace.Pi0k,
                                                     ComputePolynomialsValues(reference_element_data,
                                                                              localSpace,
                                                                              referencePoints));
    }

    /// \brief Compute the values of projections of VEM basis function derivatives at points.
    /// \param reference_element_data: an object of type \ref VEM::PCC::VEM_PCC_2D_ReferenceElement_Data which contains
    /// monomials stuff. \param localSpace: an object of type \ref VEM::PCC::VEM_PCC_2D_LocalSpace_Data which contains
    /// local matrices. \param projectionType: the \ref VEM::PCC::ProjectionTypes reporting the kind of projector used
    /// to access to the point-wise evalution of VE basis function derivatives. \param points: a matrix 3 \f$\times\f$
    /// numPoints reporting the coordinates of points. \return A vector of 2 matrices of size numPoints \f$\times\f$
    /// numDOFs whose columns contain the evaluation of the projection of each basis function derivatives with respect x
    /// and y, respectively, at points.
    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                                              const VEM_PCC_2D_LocalSpace_Data &localSpace,
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
            const std::vector<Eigen::MatrixXd> VanderDerivatives = utilities.ComputePolynomialsDerivativeValues(
                reference_element_data.Monomials, monomials, localSpace.Diameter, vander);

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

    /// \brief Compute the values of monomial basis functions at the internal quadrature points.
    /// \param localSpace: an object of type \ref VEM::PCC::VEM_PCC_2D_LocalSpace_Data which contains local matrices.
    /// \return A matrix of size numQuadrature \f$\times\f$ numMonomials whose columns contain the evaluation of
    /// monomials at the internal quadrature points.
    inline Eigen::MatrixXd ComputePolynomialsValues(const VEM_PCC_2D_LocalSpace_Data &localSpace) const
    {
        return localSpace.VanderInternal;
    }

    /// \brief Compute the values of monomial basis functions at points.
    /// \param reference_element_data: an object of type \ref VEM::PCC::VEM_PCC_2D_ReferenceElement_Data which contains
    /// monomials stuff. \param points: a matrix 3 \f$\times\f$ numPoints reporting the coordinates of points. \return A
    /// matrix of size numPoints \f$\times\f$ numMonomials whose columns contain the evaluation of monomials at points.
    inline Eigen::MatrixXd ComputePolynomialsValues(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                    const VEM_PCC_2D_LocalSpace_Data &localSpace,
                                                    const Eigen::MatrixXd &points) const
    {
        const Eigen::MatrixXd referencePoints =
            localSpace.inertia_data.FmatrixInv * (points.colwise() - localSpace.inertia_data.translation);
        return utilities.ComputePolynomialsValues(
            reference_element_data.Monomials, monomials, localSpace.Centroid, localSpace.Diameter, referencePoints);
    }

    /// \brief Compute the values of monomial basis function derivatives at the internal quadrature points.
    /// \param localSpace: an object of type \ref VEM::PCC::VEM_PCC_2D_LocalSpace_Data which contains local matrices.
    /// \return A vector of two matrices of size numQuadrature \f$\times\f$ numMonomials whose columns contain the
    /// evaluation of monomials derivatives with respect x and y, respectively, at the internal quadrature points
    inline std::vector<Eigen::MatrixXd> ComputePolynomialsDerivativeValues(
        const VEM_PCC_2D_LocalSpace_Data &localSpace) const
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
        const Eigen::MatrixXd referencePoints =
            localSpace.inertia_data.FmatrixInv * (points.colwise() - localSpace.inertia_data.translation);

        const std::vector<Eigen::MatrixXd> polynomialDerivatives = utilities.ComputePolynomialsDerivativeValues(
            reference_element_data.Monomials,
            monomials,
            localSpace.Diameter,
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

    inline Eigen::MatrixXd ComputeValuesOnEdge(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                               const Eigen::VectorXd &pointsCurvilinearCoordinates) const
    {
        Eigen::VectorXd edgeInternalPoints;
        if (reference_element_data.Quadrature.ReferenceSegmentInternalPoints.rows() > 0)
            edgeInternalPoints = reference_element_data.Quadrature.ReferenceSegmentInternalPoints.row(0).transpose();
        const Eigen::VectorXd edgeBasisCoefficients =
            utilities.ComputeEdgeBasisCoefficients(reference_element_data.Order, edgeInternalPoints);

        return utilities.ComputeValuesOnEdge(edgeInternalPoints.transpose(),
                                             reference_element_data.Order,
                                             edgeBasisCoefficients,
                                             pointsCurvilinearCoordinates);
    }

    Eigen::MatrixXd ComputeBasisFunctionsLaplacianValues(const VEM_PCC_2D_LocalSpace_Data &) const
    {
        throw std::runtime_error("Unimplemented method");
    }
    Eigen::MatrixXd ComputeBasisFunctionsLaplacianValues(const VEM_PCC_2D_ReferenceElement_Data &,
                                                         const VEM_PCC_2D_LocalSpace_Data &,
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
