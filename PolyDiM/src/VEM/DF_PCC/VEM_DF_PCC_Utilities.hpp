#ifndef __VEM_DF_PCC_Utilities_HPP
#define __VEM_DF_PCC_Utilities_HPP

#include "Eigen/Eigen"
#include "VEM_Monomials_Data.hpp"

namespace Polydim
{
namespace VEM
{
namespace DF_PCC
{
/** @file */

/// \brief Enumeration for Projector Types
enum struct ProjectionTypes
{
    Pi0km2 = 0,   ///< \f$\Pi^0_{order-1}\f$ projection to project basis
    Pi0k = 1,     ///< \f$\Pi^0_{order}\f$ projection to project basis
    PiNabla = 2,  ///< \f$\Pi^{\nabla}_{order-1}\f$ projection to project basis gradient
    Pi0km1Der = 3 ///< \f$\Pi^{0}_{order-1}\f$ projection to project basis gradient
};

/// \brief Base class for computing values of basis functions of Primal Conforming Constant degree
/// Virtual Element Methods.
/// \copyright See top level LICENSE file for details.
template <unsigned short dimension> struct VEM_DF_PCC_Utilities final
{
    /// Compute the Edge basis coefficients
    Eigen::VectorXd ComputeEdgeBasisCoefficients(const unsigned int &order, const Eigen::VectorXd &edgeInternalPoints) const;

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
                                                                       const std::vector<Eigen::MatrixXd> &piNabla,
                                                                       const std::vector<Eigen::MatrixXd> &pi0km1Der) const
    {
        switch (projectionType)
        {
        case ProjectionTypes::PiNabla: {
            std::vector<Eigen::MatrixXd> basisFunctionDerivativeValues(dimension * dimension);

            for (unsigned short j = 0; j < dimension; ++j)
                for (unsigned short i = 0; i < dimension; ++i)
                    basisFunctionDerivativeValues[dimension * j + i] = vanderInternalDerivatives[i] * piNabla[j];

            return basisFunctionDerivativeValues;
        }
        case ProjectionTypes::Pi0km1Der: {
            std::vector<Eigen::MatrixXd> basisFunctionDerivativeValues(dimension * dimension);

            for (unsigned short j = 0; j < dimension; ++j)
                for (unsigned short i = 0; i < dimension; ++i)
                    basisFunctionDerivativeValues[dimension * j + i] =
                        vanderInternal.leftCols(Nkm1) * pi0km1Der[dimension * j + i];

            return basisFunctionDerivativeValues;
        }
        default:
            throw std::runtime_error("Unknown projector type");
        }
    }

    inline Eigen::MatrixXd ComputeBasisFunctionsDivergenceValues(const unsigned int &Nkm1,
                                                                 const Eigen::MatrixXd &vanderInternal,
                                                                 const Eigen::MatrixXd &vmatrix) const
    {
        return vanderInternal.leftCols(Nkm1) * vmatrix;
    }

    /// \brief Compute the values of the polynomial projection of basis functions at internal
    /// quadrature points on the geometry.
    /// \param basisFunctionsValues The matrix of values. Each column will contain the values of a
    /// basis function's projection at internal quadrature points.
    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsValues(const ProjectionTypes &projectionType,
                                                                    const unsigned int &Nkm2,
                                                                    const std::vector<Eigen::MatrixXd> &pi0km2,
                                                                    const std::vector<Eigen::MatrixXd> &pi0k,
                                                                    const Eigen::MatrixXd &vanderInternal) const
    {
        std::vector<Eigen::MatrixXd> basisFunctionValues(dimension);
        switch (projectionType)
        {
        case ProjectionTypes::Pi0km2:
            for (unsigned short i = 0; i < dimension; ++i)
                basisFunctionValues[i] = vanderInternal.leftCols(Nkm2) * pi0km2[i];
            break;
        case ProjectionTypes::Pi0k:
            for (unsigned short i = 0; i < dimension; ++i)
                basisFunctionValues[i] = vanderInternal * pi0k[i];
            break;
        default:
            throw std::runtime_error("Unknown projector type");
        }
        return basisFunctionValues;
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

    Eigen::MatrixXd ComputeDofiDofiStabilizationMatrix(const std::vector<Eigen::MatrixXd> &projector,
                                                       const double &coefficient,
                                                       const std::vector<Eigen::MatrixXd> &dmatrix) const;
};
} // namespace DF_PCC
} // namespace VEM
} // namespace Polydim

#endif
