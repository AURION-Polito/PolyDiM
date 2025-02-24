#ifndef __VEM_PCC_Utilities_HPP
#define __VEM_PCC_Utilities_HPP

#include "Eigen/Eigen"
#include "VEM_Monomials_Data.hpp"

namespace Polydim
{
namespace VEM
{
namespace PCC
{
/** @file */

/// \brief Enumeration for Projector Types
enum struct ProjectionTypes
{
    Pi0km1 = 0,    ///< \f$\Pi^0_{order-1}\f$ projection to project basis
    Pi0k = 1,      ///< \f$\Pi^0_{order}\f$ projection to project basis
    PiNabla = 2,   ///< \f$\Pi^{\nabla}_{order-1}\f$ projection to project basis gradient
    Pi0km1Der = 3, ///< \f$\Pi^{0}_{order-1}\f$ projection to project basis gradient
    Pi0klm1 = 4
};

/// \brief Base class for computing values of basis functions of Primal Conforming Constant degree
/// Virtual Element Methods.
/// \copyright See top level LICENSE file for details.
template <unsigned short dimension> struct VEM_PCC_Utilities final
{
    /// Compute the Edge basis coefficients
    Eigen::VectorXd ComputeEdgeBasisCoefficients(const unsigned int &order, const Eigen::VectorXd &edgeInternalPoints) const;

    /// \brief Compute matrices of coefficients related to \ref VEM::PCC::ProjectionTypes::Pi0km1 and \ref
    /// VEM::PCC::ProjectionTypes::Pi0k. \param[in] measure: a double representing the measure of the element \param[in]
    /// order: the vem order \param[in] Nkm1: number of polynomials of degree order-1 \param[in] Nk: number of
    /// polynomials of degree order \param[in] NumInternalBasisFunctions: number of internal basis functions \param[in]
    /// NumBasisFunctions: number of internal functions \param[in] Hmatrix: mass matrix of polynomials \param[in]
    /// PiNabla: \ref VEM::PCC::ProjectionTypes::PiNabla coefficients \param[out] Cmatrix: the right-hand side to define
    /// \ref VEM::PCC::ProjectionTypes::Pi0k coefficients \param[out] Pi0km1: \ref VEM::PCC::ProjectionTypes::Pi0km1
    /// coefficients \param[out] Pi0k: \ref VEM::PCC::ProjectionTypes::Pi0k coefficients
    void ComputeL2Projectors(const double &measure,
                             const unsigned int &order,
                             const unsigned int &Nkm1,
                             const unsigned int &Nk,
                             const unsigned int &NumInternalBasisFunctions,
                             const unsigned int &NumBasisFunctions,
                             const Eigen::MatrixXd &Hmatrix,
                             const Eigen::MatrixXd &PiNabla,
                             Eigen::MatrixXd &Cmatrix,
                             Eigen::MatrixXd &Pi0km1,
                             Eigen::MatrixXd &Pi0k) const;

    /// \brief Compute the values of the polynomial projection of
    /// derivatives of basis functions at internal quadrature points on
    /// the geometry.
    /// \param basisFunctionsDerivativeValues The vector of matrices of
    /// values. Its length equals \ref Dimension(). Each column of each
    /// matrix will contain the values of a basis function's projected
    /// derivative at internal quadrature points.
    std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const ProjectionTypes &projectionType,
                                                                       const unsigned int &Nkm1,
                                                                       const Eigen::MatrixXd &vanderInternal,
                                                                       const std::vector<Eigen::MatrixXd> &vanderInternalDerivatives,
                                                                       const Eigen::MatrixXd &piNabla,
                                                                       const std::vector<Eigen::MatrixXd> &pi0km1Der) const
    {
        switch (projectionType)
        {
        case ProjectionTypes::PiNabla: {
            std::vector<Eigen::MatrixXd> basisFunctionsDerivativeValues(dimension);

            for (unsigned short i = 0; i < dimension; ++i)
                basisFunctionsDerivativeValues[i] = vanderInternalDerivatives[i] * piNabla;

            return basisFunctionsDerivativeValues;
        }
        case ProjectionTypes::Pi0km1Der: {
            std::vector<Eigen::MatrixXd> basisFunctionsDerivativeValues(dimension);

            for (unsigned short i = 0; i < dimension; ++i)
                basisFunctionsDerivativeValues[i] = vanderInternal.leftCols(Nkm1) * pi0km1Der[i];

            return basisFunctionsDerivativeValues;
        }
        default:
            throw std::runtime_error("Unknown projector type");
        }
    }

    /// \brief Compute the values of the polynomial projection of
    /// the laplacian of basis functions at internal quadrature
    /// points on the geometry.
    /// \details The projection of the laplacian considered is \f$
    /// \widetilde{\Delta}\varphi = \nabla
    /// \cdot \left(\Pi^0_{\mathrm{order}-1} \nabla \varphi \right)\f$.
    /// \param basisFunctionLaplacianValues The matrix of
    /// values. Each column will contain the values
    /// of the projected laplacian of a basis function at internal quadrature points.
    /// \sa \ref ComputeInternalQuadratureWeights().
    Eigen::MatrixXd ComputeBasisFunctionsLaplacianValues(const ProjectionTypes &projectionType,
                                                         const unsigned int &Nkm1,
                                                         const std::vector<Eigen::MatrixXd> &vanderInternalDerivatives,
                                                         const std::vector<Eigen::MatrixXd> &pi0km1Der) const
    {
        switch (projectionType)
        {
        case ProjectionTypes::Pi0km1Der: {
            Eigen::MatrixXd basisFunctionsLaplacianValues = vanderInternalDerivatives[0].leftCols(Nkm1) * pi0km1Der[0];
            for (unsigned int d = 1; d < dimension; ++d)
                basisFunctionsLaplacianValues += vanderInternalDerivatives[d].leftCols(Nkm1) * pi0km1Der[d];

            return basisFunctionsLaplacianValues;
        }
        default:
            throw std::runtime_error("Unknown projector type");
        }
    }

    /// \brief Compute the values of the polynomial projection of basis functions at internal
    /// quadrature points on the geometry.
    /// \param basisFunctionsValues The matrix of values. Each column will contain the values of a
    /// basis function's projection at internal quadrature points.
    inline Eigen::MatrixXd ComputeBasisFunctionsValues(const ProjectionTypes &projectionType,
                                                       const unsigned int &Nkm1,
                                                       const Eigen::MatrixXd &pi0km1,
                                                       const Eigen::MatrixXd &pi0k,
                                                       const Eigen::MatrixXd &vanderInternal) const
    {
        switch (projectionType)
        {
        case ProjectionTypes::Pi0km1:
            return vanderInternal.leftCols(Nkm1) * pi0km1;
        case ProjectionTypes::Pi0k:
            return vanderInternal * pi0k;
        default:
            throw std::runtime_error("Unknown projector type");
        }
    }

    /// \brief Compute the values of the basis functions of the polynomial basis of projectors at
    /// internal quadrature points on the geometry.
    /// a basis function at internal quadrature points.
    inline Eigen::MatrixXd ComputePolynomialsValues(const Eigen::MatrixXd &vanderInternal) const
    {
        return vanderInternal;
    }

    /// \brief Compute the values of the basis functions of the polynomial basis of projectors at
    /// given points on the geometry.
    /// \param points The points at which to evaluate the basis functions.
    /// a basis function at the given points.
    template <typename VEM_MonomialType>
    inline Eigen::MatrixXd ComputePolynomialsValues(const Monomials::VEM_Monomials_Data &data,
                                                    const VEM_MonomialType &monomials,
                                                    const Eigen::Vector3d &centroid,
                                                    const double &diameter,
                                                    const Eigen::MatrixXd &points) const
    {
        return monomials.Vander(data, points, centroid, diameter);
    }

    /// \brief Compute the values of the derivatives of the basis functions of the polynomial
    /// basis of projectors at internal quadrature points on the geometry.
    /// values of a basis function's derivative at internal quadrature points.
    inline std::vector<Eigen::MatrixXd> ComputePolynomialsDerivativeValues(const std::vector<Eigen::MatrixXd> &vanderInternalDerivatives) const
    {
        return vanderInternalDerivatives;
    }

    /// \brief Compute the values of the derivatives of the basis functions of the polynomial
    /// basis of projectors at given points on the geometry.
    /// \param points The points at which to evaluate the derivatives.
    /// \ref Dimension(). Each column of each matrix will contain the values of a basis function's
    /// derivative at the given points.
    template <typename VEM_MonomialType>
    inline std::vector<Eigen::MatrixXd> ComputePolynomialsDerivativeValues(const Monomials::VEM_Monomials_Data &data,
                                                                           const VEM_MonomialType &monomials,
                                                                           const double &diameter,
                                                                           const Eigen::MatrixXd &vander) const
    {
        return monomials.VanderDerivatives(data, vander, diameter);
    }

    /// \brief Compute the values of the laplacian of the basis functions of the polynomial basis
    /// of projectors at given points.
    /// \param points The points at which to evaluate the laplacian.
    /// values of a basis function's laplacian at the given points.
    template <typename VEM_MonomialType>
    inline Eigen::MatrixXd ComputePolynomialsLaplacianValues(const Monomials::VEM_Monomials_Data &data,
                                                             const VEM_MonomialType &monomials,
                                                             const double &diameter,
                                                             const Eigen::MatrixXd &vander) const
    {
        return monomials.VanderLaplacian(data, vander, diameter);
    }

    /// \brief Compute basis function values at given points on an edge.
    /// \details
    /// The first two columns of the output contain the values of the basis functions relative to
    /// the edge extrema, the others contain the values of basis functions relative to edge
    /// internal dofs.
    /// \param edgeInternalPoints reference points internal to each edge, also used as quadrature points.
    /// \param pointsCurvilinearCoordinates Curvilinear coordinates of the points, expressed in
    /// curvilinear coordinates in the interval [0,1].
    /// \param values Matrix containing on each column the values of a local basis function at the
    /// given points, as described above.
    /// \note This is compatible with the 1D points returned by quadrature rules.
    Eigen::MatrixXd ComputeValuesOnEdge(const Eigen::RowVectorXd &edgeInternalPoints,
                                        const unsigned int &order,
                                        const Eigen::VectorXd &edgeBasisCoefficients,
                                        const Eigen::VectorXd &pointsCurvilinearCoordinates) const;

    /// \brief Compute the dofi-dofi stabilization matrix using \ref VEM::PCC::ProjectionTypes projector.
    /// \param piNabla: the projection coefficients related to \ref VEM::PCC::ProjectionTypes
    /// \param diameter: a double representing the diameter of the element
    /// \param DMatrix: matrix containing the polynomials degrees of freedom
    /// \return The dofi-dofi stabilization matrix.
    Eigen::MatrixXd ComputeDofiDofiStabilizationMatrix(const Eigen::MatrixXd &projector,
                                                       const double &coefficient,
                                                       const Eigen::MatrixXd &Dmatrix) const;
};
} // namespace PCC
} // namespace VEM
} // namespace Polydim

#endif
