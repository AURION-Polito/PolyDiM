#ifndef __I_VEM_PCC_2D_LocalSpace_HPP
#define __I_VEM_PCC_2D_LocalSpace_HPP

#include "Eigen/Eigen"
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

class I_VEM_PCC_2D_LocalSpace
{
  public:
    virtual ~I_VEM_PCC_2D_LocalSpace()
    {
    }

    /// \brief Create and Initialize all the variables contained in \ref VEM::PCC::VEM_PCC_2D_LocalSpace_Data
    /// \param reference_element_data: an object of type \ref VEM::PCC::VEM_PCC_2D_ReferenceElement_Data which contains
    /// monomials, quadrature and the number of degrees of freedom, counting in order DOFS associated with vertices,
    /// edges and internal values.
    /// \param polygon: an object of type \ref VEM::PCC::VEM_PCC_2D_Polygon_Geometry which
    /// contains the geoemtric properties of the elements.
    /// \return An object of type \ref VEM::PCC::VEM_PCC_2D_LocalSpace_Data.
    virtual VEM_PCC_2D_LocalSpace_Data CreateLocalSpace(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                        const VEM_PCC_2D_Polygon_Geometry &polygon) const = 0;

    /// \brief Compute the values of projections of VEM basis functions at the internal quadrature points.
    /// \param localSpace: an object of type \ref VEM::PCC::VEM_PCC_2D_LocalSpace_Data which contains local matrices.
    /// \param projectionType: the \ref VEM::PCC::ProjectionTypes reporting the kind of projector used to access to the
    /// point-wise evalution of VE basis functions. \return A matrix of size numQuadrature \f$\times\f$ numDOFs whose
    /// columns contain the evaluation of the projection of each basis function at the internal quadrature points.
    virtual Eigen::MatrixXd ComputeBasisFunctionsValues(const VEM_PCC_2D_LocalSpace_Data &localSpace,
                                                        const ProjectionTypes &projectionType) const = 0;

    /// \brief Compute the values of projections of VEM basis function derivatives at the internal quadrature points.
    /// \param localSpace: an object of type \ref VEM::PCC::VEM_PCC_2D_LocalSpace_Data which contains local matrices.
    /// \param projectionType: the \ref VEM::PCC::ProjectionTypes reporting the kind of projector used to access to the
    /// point-wise evalution of VE basis function derivatives. \return A vector of 2 matrices of size numQuadrature
    /// \f$\times\f$ numDOFs whose columns contain the evaluation of the projection of each basis function derivatives
    /// with respect x and y, respectively, at the internal quadrature points.
    virtual std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(
        const VEM_PCC_2D_LocalSpace_Data &localSpace, const ProjectionTypes &projectionType) const = 0;

    /// \brief Compute the values of VEM basis function laplacian at the internal quadrature points, here approximated
    /// using \ref VEM::PCC::ProjectionTypes::Pi0km1Der. \param localSpace: an object of type \ref
    /// VEM::PCC::VEM_PCC_2D_LocalSpace_Data which contains local matrices. \return A matrix of size numQuadrature
    /// \f$\times\f$ numDOFs whose columns contain the evaluation of the approximated laplacian at the internal
    /// quadrature points.
    virtual Eigen::MatrixXd ComputeBasisFunctionsLaplacianValues(
        const VEM_PCC_2D_LocalSpace_Data &localSpace) const = 0;

    /// \brief Compute the values of projections of VEM basis functions at points.
    /// \param reference_element_data: an object of type \ref VEM::PCC::VEM_PCC_2D_ReferenceElement_Data which contains
    /// monomials stuff. \param localSpace: an object of type \ref VEM::PCC::VEM_PCC_2D_LocalSpace_Data which contains
    /// local matrices. \param projectionType: the \ref VEM::PCC::ProjectionTypes reporting the kind of projector used
    /// to access to the point-wise evalution of VE basis functions. \param points: a matrix 3 \f$\times\f$ numPoints
    /// reporting the coordinates of points. \return A matrix of size numPoints \f$\times\f$ numDOFs whose columns
    /// contain the evaluation of the projection of each basis function at points.
    virtual Eigen::MatrixXd ComputeBasisFunctionsValues(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                        const VEM_PCC_2D_LocalSpace_Data &localSpace,
                                                        const ProjectionTypes &projectionType,
                                                        const Eigen::MatrixXd &points) const = 0;

    /// \brief Compute the values of projections of VEM basis function derivatives at points.
    /// \param reference_element_data: an object of type \ref VEM::PCC::VEM_PCC_2D_ReferenceElement_Data which contains
    /// monomials stuff. \param localSpace: an object of type \ref VEM::PCC::VEM_PCC_2D_LocalSpace_Data which contains
    /// local matrices. \param projectionType: the \ref VEM::PCC::ProjectionTypes reporting the kind of projector used
    /// to access to the point-wise evalution of VE basis function derivatives. \param points: a matrix 3 \f$\times\f$
    /// numPoints reporting the coordinates of points. \return A vector of 2 matrices of size numPoints \f$\times\f$
    /// numDOFs whose columns contain the evaluation of the projection of each basis function derivatives with respect x
    /// and y, respectively, at points.
    virtual std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(
        const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
        const VEM_PCC_2D_LocalSpace_Data &localSpace,
        const ProjectionTypes &projectionType,
        const Eigen::MatrixXd &points) const = 0;

    /// \brief Compute the values of VEM basis function laplacian at points, here approximated using \ref
    /// VEM::PCC::ProjectionTypes::Pi0km1Der. \param reference_element_data: an object of type \ref
    /// VEM::PCC::VEM_PCC_2D_ReferenceElement_Data which contains monomials stuff. \param localSpace: an object of type
    /// \ref VEM::PCC::VEM_PCC_2D_LocalSpace_Data which contains local matrices. \param points: a matrix 3 \f$\times\f$
    /// numPoints reporting the coordinates of points. \return A matrix of size numPoints \f$\times\f$ numDOFs whose
    /// columns contain the evaluation of the approximated laplacian at points.
    virtual Eigen::MatrixXd ComputeBasisFunctionsLaplacianValues(
        const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
        const VEM_PCC_2D_LocalSpace_Data &localSpace,
        const Eigen::MatrixXd &points) const = 0;

    /// \brief Compute the values of monomial basis functions at the internal quadrature points.
    /// \param localSpace: an object of type \ref VEM::PCC::VEM_PCC_2D_LocalSpace_Data which contains local matrices.
    /// \return A matrix of size numQuadrature \f$\times\f$ numMonomials whose columns contain the evaluation of
    /// monomials at the internal quadrature points.
    virtual Eigen::MatrixXd ComputePolynomialsValues(const VEM_PCC_2D_LocalSpace_Data &localSpace) const = 0;

    /// \brief Compute the values of monomial basis functions at points.
    /// \param reference_element_data: an object of type \ref VEM::PCC::VEM_PCC_2D_ReferenceElement_Data which contains
    /// monomials stuff. \param points: a matrix 3 \f$\times\f$ numPoints reporting the coordinates of points. \return A
    /// matrix of size numPoints \f$\times\f$ numMonomials whose columns contain the evaluation of monomials at points.
    virtual Eigen::MatrixXd ComputePolynomialsValues(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                     const VEM_PCC_2D_LocalSpace_Data &localSpace,
                                                     const Eigen::MatrixXd &points) const = 0;

    /// \brief Compute the values of monomial basis function derivatives at the internal quadrature points.
    /// \param localSpace: an object of type \ref VEM::PCC::VEM_PCC_2D_LocalSpace_Data which contains local matrices.
    /// \return A vector of two matrices of size numQuadrature \f$\times\f$ numMonomials whose columns contain the
    /// evaluation of monomials derivatives with respect x and y, respectively, at the internal quadrature points
    virtual std::vector<Eigen::MatrixXd> ComputePolynomialsDerivativeValues(
        const VEM_PCC_2D_LocalSpace_Data &localSpace) const = 0;

    /// \brief Compute the values of monomial basis functions at points.
    /// \param reference_element_data: an object of type \ref VEM::PCC::VEM_PCC_2D_ReferenceElement_Data which contains
    /// monomials stuff. \param points: a matrix 3 \f$\times\f$ numPoints reporting the coordinates of points. \return A
    /// vector of two matrices of size numPoints \f$\times\f$ numMonomials whose columns contain the evaluation of
    /// monomials derivatives with respect x and y, respectively, at points
    virtual std::vector<Eigen::MatrixXd> ComputePolynomialsDerivativeValues(
        const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
        const VEM_PCC_2D_LocalSpace_Data &localSpace,
        const Eigen::MatrixXd &points) const = 0;

    /// \brief Compute the values of monomials laplacian at points.
    /// \param reference_element_data: an object of type \ref VEM::PCC::VEM_PCC_2D_ReferenceElement_Data which contains
    /// monomials stuff. \param points: a matrix 3 \f$\times\f$ numPoints reporting the coordinates of points. \return A
    /// matrix of size numPoints \f$\times\f$ numMonomials whose columns contain the evaluation of monomials at points.
    virtual Eigen::MatrixXd ComputePolynomialsLaplacianValues(
        const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
        const VEM_PCC_2D_LocalSpace_Data &localSpace,
        const Eigen::MatrixXd &points) const = 0;

    virtual Eigen::MatrixXd ComputeValuesOnEdge(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data,
                                                const Eigen::VectorXd &pointsCurvilinearCoordinates) const = 0;
};
} // namespace PCC
} // namespace VEM
} // namespace Polydim

#endif
