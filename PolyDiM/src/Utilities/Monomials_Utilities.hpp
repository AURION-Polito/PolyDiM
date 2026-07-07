// _LICENSE_HEADER_
//
// Copyright (C) 2019 - 2025.
// Terms register on the GPL-3.0 license.
//
// This file can be redistributed and/or modified under the license terms.
//
// See top level LICENSE file for more details.
//
// This file can be used citing references in CITATION.cff file.

#ifndef __Monomials_Utilities_HPP
#define __Monomials_Utilities_HPP

#include "LAPACK_utilities.hpp"
#include "Monomials_Data.hpp"

namespace Polydim
{
namespace Utilities
{
template <unsigned short dimension> struct Monomials_Utilities final
{
    /// @brief Collect the monomial exponents into a single matrix.
    ///
    /// @param data Monomial basis data.
    /// @return A \f$\text{dimension} \times N_{\text{mon}}\f$ integer matrix whose
    ///         \f$m\f$-th column is the exponent multi-index \f$\alpha\f$ of the \f$m\f$-th monomial.
    Eigen::MatrixXi Exponents(const Polydim::Utilities::Monomials_Data &data) const
    {
        Eigen::MatrixXi exponents(dimension, data.NumMonomials);

        for (unsigned int m = 0; m < data.NumMonomials; m++)
            exponents.col(m) << data.Exponents[m];

        return exponents;
    }

    /// @brief Evaluate the Vandermonde matrix of the scaled monomial basis.
    ///
    /// Builds the matrix whose \f$(p, m)\f$ entry is the \f$m\f$-th scaled monomial
    /// evaluated at the \f$p\f$-th point,
    /// \f$\prod_d \left(\frac{x_d - x_{E,d}}{h_E}\right)^{\alpha_d}\f$. Per-coordinate
    /// integer powers are precomputed once (@c VanderPartial) and combined according
    /// to each monomial's exponent multi-index.
    ///
    /// @param data     Monomial basis data (degree, exponents, number of monomials).
    /// @param points   Evaluation points, one per column.
    /// @param centroid Element centroid \f$\boldsymbol{x}_E\f$.
    /// @param diam     Element diameter \f$h_E\f$ used to scale the coordinates.
    /// @return A \f$N_{\text{points}} \times N_{\text{mon}}\f$ Vandermonde matrix;
    ///         a column of ones is returned when the basis reduces to the constant.
    Eigen::MatrixXd Vander(const Polydim::Utilities::Monomials_Data &data,
                           const Eigen::MatrixXd &points,
                           const Eigen::Vector3d &centroid,
                           const double &diam) const
    {
        Eigen::MatrixXd vander;
        const unsigned int numPoints = points.cols();
        if (data.NumMonomials > 1)
        {
            // VanderPartial[i]'s rows contain (x-x_E)^i/h_E^i,
            // (y-y_E)^i/h_E^i and (possibly) (z-z_E)^i/h_E^i respectively.
            // Size is dimension x numPoints.
            std::vector<Eigen::MatrixXd> VanderPartial(data.PolynomialDegree + 1, Eigen::MatrixXd(dimension, numPoints));
            double inverseDiam = 1.0 / diam;
            VanderPartial[0].setOnes(dimension, numPoints);
            VanderPartial[1] = (points.colwise() - centroid) * inverseDiam;

            for (unsigned int i = 2; i <= data.PolynomialDegree; i++)
                VanderPartial[i] = VanderPartial[i - 1].cwiseProduct(VanderPartial[1]);

            vander.resize(numPoints, data.NumMonomials);
            vander.col(0).setOnes();
            for (unsigned int i = 1; i < data.NumMonomials; ++i)
            {
                const Eigen::VectorXi expo = data.Exponents[i];

                vander.col(i) = (VanderPartial[expo[0]].row(0)).transpose();
                if (dimension > 1)
                    vander.col(i) = vander.col(i).cwiseProduct(VanderPartial[expo[1]].row(1).transpose());
                if (dimension > 2)
                    vander.col(i) = vander.col(i).cwiseProduct(VanderPartial[expo[2]].row(2).transpose());
            }
        }
        else
            vander.setOnes(numPoints, 1);

        return vander;
    }

    /// @brief Evaluate the Vandermonde matrices of the first-order partial derivatives.
    ///
    /// For each spatial direction \f$i\f$, returns the Vandermonde matrix of
    /// \f$\partial_{x_i} m_\alpha\f$, reusing the already-evaluated monomial values
    /// in @p Vander. Each derivative maps to a lower-degree monomial (via
    /// @c DerivativeIndices) scaled by the derivative coefficient and by the chain-rule
    /// factor \f$1/h_E\f$; directions with no contribution yield a zero column.
    ///
    /// @tparam MonomialType Monomial type exposing @c DerivativeIndices and @c DerivativeMatrix.
    /// @param data      Monomial basis data.
    /// @param monomials Monomial object providing the derivative maps.
    /// @param Vander    Vandermonde matrix of the monomials at the same points (see Vander()).
    /// @param diam      Element diameter \f$h_E\f$.
    /// @return A vector of @c dimension matrices, the \f$i\f$-th being the Vandermonde
    ///         matrix of \f$\partial_{x_i}\f$ of the basis (same shape as @p Vander).
    template <typename MonomialType>
    std::vector<Eigen::MatrixXd> VanderDerivatives(const Polydim::Utilities::Monomials_Data &data,
                                                   const MonomialType &monomials,
                                                   const Eigen::MatrixXd &Vander,
                                                   const double &diam) const
    {
        std::vector<Eigen::MatrixXd> vanderDerivatives;
        vanderDerivatives.resize(dimension);
        for (unsigned int i = 0; i < dimension; i++)
        {
            vanderDerivatives[i].resizeLike(Vander);
            vanderDerivatives[i].col(0).setZero();
        }
        if (data.NumMonomials > 1)
        {
            double inverseDiam = 1.0 / diam;
            for (unsigned int k = 1; k < data.NumMonomials; k++)
            {
                std::vector<int> derIndices = monomials.DerivativeIndices(data, k);
                for (unsigned int i = 0; i < dimension; i++)
                {
                    if (derIndices[i] >= 0)
                        vanderDerivatives[i].col(k) =
                            inverseDiam * monomials.DerivativeMatrix(data, i)(k, derIndices[i]) * Vander.col(derIndices[i]);
                    else
                        vanderDerivatives[i].col(k).setZero();
                }
            }
        }

        return vanderDerivatives;
    }

    /// @brief Evaluate the Vandermonde matrix of the monomial Laplacian.
    ///
    /// Returns the matrix whose \f$k\f$-th column holds \f$\Delta m_k\f$ evaluated at
    /// the points, assembled from the precomputed @c data.Laplacian coefficients and
    /// the second-derivative maps (@c SecondDerivativeIndices), scaled by
    /// \f$1/h_E^2\f$. The constant and linear monomials, whose Laplacian vanishes,
    /// give zero columns.
    ///
    /// @tparam MonomialType Monomial type exposing @c SecondDerivativeIndices.
    /// @param data      Monomial basis data (must provide the @c Laplacian matrix).
    /// @param monomials Monomial object providing the second-derivative maps.
    /// @param Vander    Vandermonde matrix of the monomials at the same points (see Vander()).
    /// @param diam      Element diameter \f$h_E\f$.
    /// @return A matrix, shaped like @p Vander, holding \f$\Delta\f$ of the basis.
    template <typename MonomialType>
    Eigen::MatrixXd VanderLaplacian(const Polydim::Utilities::Monomials_Data &data,
                                    const MonomialType &monomials,
                                    const Eigen::MatrixXd &Vander,
                                    const double &diam) const
    {
        Eigen::MatrixXd vanderLaplacian;

        vanderLaplacian.resizeLike(Vander);
        vanderLaplacian.block(0, 0, Vander.rows(), 3).setZero();
        Eigen::MatrixXd laplacian = data.Laplacian;

        if (data.NumMonomials > 3)
        {
            const double inverseDiamSqrd = 1.0 / (diam * diam);
            for (unsigned int k = 3; k < data.NumMonomials; k++)
            {
                std::vector<int> secondDerIndices = monomials.SecondDerivativeIndices(data, k);
                if (secondDerIndices[0] >= 0)
                    vanderLaplacian.col(k) =
                        inverseDiamSqrd * laplacian(k, secondDerIndices[0]) * Vander.col(secondDerIndices[0]);
                else
                    vanderLaplacian.col(k).setZero();
                for (unsigned int i = 1; i < dimension; i++)
                {
                    if (secondDerIndices[i] >= 0)
                        vanderLaplacian.col(k) +=
                            inverseDiamSqrd * laplacian(k, secondDerIndices[i]) * Vander.col(secondDerIndices[i]);
                }
            }
        }

        return vanderLaplacian;
    }

    /// @brief \f$L^2(E)\f$-orthonormalize the monomial basis via modified Gram–Schmidt.
    ///
    /// Orthonormalizes the scaled monomials with respect to the weighted
    /// \f$L^2(E)\f$ inner product defined by the quadrature @p weights, using a
    /// modified Gram–Schmidt factorization followed by a re-orthogonalization step
    /// for numerical stability. The routine returns the change-of-basis matrices
    /// mapping between the monomial and the orthonormal basis.
    ///
    /// @param weights   Quadrature weights defining the \f$L^2(E)\f$ inner product.
    /// @param Vander     Vandermonde matrix of the monomials at the quadrature points.
    /// @param[out] Hmatrix    Mass matrix in the orthonormal basis, \f$Q_2^\top Q_2\f$ (close to the identity).
    /// @param[out] QmatrixInv Lower-triangular inverse change-of-basis matrix, \f$(R_2 R_1)^\top\f$.
    /// @param[out] Qmatrix    Change-of-basis matrix to the orthonormal basis (inverse of @p QmatrixInv).
    void MGSOrthonormalize(const Eigen::VectorXd &weights,
                           const Eigen::MatrixXd &Vander,
                           Eigen::MatrixXd &Hmatrix,
                           Eigen::MatrixXd &QmatrixInv,
                           Eigen::MatrixXd &Qmatrix) const
    {
        Eigen::MatrixXd Q1;
        Eigen::MatrixXd R1;
        LAPACK_utilities::MGS(Vander, Q1, R1);

        // L2(E)-re-orthogonalization process
        Eigen::MatrixXd Q2;
        Eigen::MatrixXd R2;
        LAPACK_utilities::MGS(weights.array().sqrt().matrix().asDiagonal() * Q1, Q2, R2);

        Hmatrix = Q2.transpose() * Q2;

        QmatrixInv = (R2 * R1).transpose();
        LAPACK_utilities::inverseTri(QmatrixInv, Qmatrix, 'L', 'N');
    }
};
} // namespace Utilities
} // namespace Polydim

#endif
