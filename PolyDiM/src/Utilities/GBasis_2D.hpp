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

#ifndef __GBasis_2D_HPP
#define __GBasis_2D_HPP

#include "GBasis_Data.hpp"
#include "Monomials_2D.hpp"

namespace Polydim
{
namespace Utilities
{
class GBasis_2D final
{
private:
    /// @brief 2D scalar monomial basis \f$\{m_\alpha\}\f$ used to build the vector basis.
    Polydim::Utilities::Monomials_2D monomials;


    /// @brief Compute the decomposition indices of a single vector monomial onto the G-basis.
    ///
    /// Given the multi-index @p expo of a scalar monomial, returns the list of
    /// index pairs \f$(i,j)\f$ that describe how the associated vector-valued
    /// monomial is expressed in terms of the \f$\mathcal{G}^\nabla\f$ /
    /// \f$\mathcal{G}^\oplus\f$ decomposition stored in @p data. Each entry is an
    /// `Eigen::Vector2i` holding the two component indices of the decomposition.
    ///
    /// @param data Precomputed G-basis data (dimensions, exponent map, vector decomposition).
    /// @param expo Exponent multi-index \f$\alpha = (\alpha_1, \alpha_2)\f$ of the scalar monomial.
    /// @return Index pairs describing the decomposition of the vector monomial.
    /// @note Internal helper invoked by Compute(); not part of the public interface.
    std::vector<Eigen::Vector2i> VectorDecompositionIndices(const Polydim::Utilities::GBasis_Data &data,
                                                            const Eigen::VectorXi &expo) const;

public:
    /// @brief Build the vector polynomial G-basis of a given degree.
    ///
    /// Assembles all quantities needed to represent the vector polynomial space
    /// \f$[\mathbb{P}_k]^2\f$ through its decomposition
    /// \f$\mathcal{G}_k^\nabla \oplus \mathcal{G}_k^\oplus\f$, where
    /// \f$\mathcal{G}_k^\nabla = \nabla \mathbb{P}_{k+1}\f$ and
    /// \f$\mathcal{G}_k^\oplus = \boldsymbol{x}^{\perp}\,\mathbb{P}_{k-1}\f$
    /// with \f$\boldsymbol{x}^{\perp} = (x_2, -x_1)\f$.
    ///
    /// @param polynomial_degree Maximum polynomial degree \f$k\f$ of the basis.
    /// @return A fully populated Polydim::Utilities::GBasis_Data structure, including the
    ///         space dimensions (\f$N_k\f$, \f$N_{k-1}\f$, \f$N_{k+1}\f$, and the sizes of the
    ///         \f$\mathcal{G}^\oplus\f$ and \f$\mathcal{G}^\nabla\f$ blocks), the exponent
    ///         matrix, and the vector-decomposition maps.
    Polydim::Utilities::GBasis_Data Compute(const unsigned int polynomial_degree);

    /// @brief Evaluate the Vandermonde matrix of the \f$\mathcal{G}_k^\oplus\f$ block.
    ///
    /// Builds the two Cartesian components of the Vandermonde matrix of
    /// \f$\mathcal{G}_k^\oplus = \boldsymbol{x}^{\perp}\,\mathbb{P}_{k-1}\f$,
    /// obtained by multiplying each scalar monomial of \f$\mathbb{P}_{k-1}\f$ by
    /// \f$\boldsymbol{x}^{\perp} = (x_2, -x_1)\f$. Concretely, the first
    /// \f$N_{k-1}\f$ columns of @p vander (the \f$\mathbb{P}_{k-1}\f$ monomials)
    /// are scaled column-wise by \f$x_2\f$ for the first component and by
    /// \f$-x_1\f$ for the second, where \f$x_1\f$ and \f$x_2\f$ are the linear
    /// monomials stored in columns 1 and 2 of @p vander.
    ///
    /// @param data   Precomputed G-basis data; only @c Nkm1 (\f$N_{k-1}\f$) is used here.
    /// @param vander Vandermonde matrix of the scalar monomial basis evaluated at the
    ///               evaluation points (one row per point); column 1 holds \f$x_1\f$,
    ///               column 2 holds \f$x_2\f$, and the leading \f$N_{k-1}\f$ columns hold
    ///               the \f$\mathbb{P}_{k-1}\f$ monomials.
    /// @return A vector of two `Eigen::MatrixXd` (x- and y-components), each of size
    ///         @c vander.rows() \f$\times\f$ @c data.Nkm1, forming the Vandermonde matrix
    ///         of \f$\mathcal{G}_k^\oplus\f$.
    std::vector<Eigen::MatrixXd> VanderGBigOPlus(const Polydim::Utilities::GBasis_Data &data, const Eigen::MatrixXd &vander) const
    {
        std::vector<Eigen::MatrixXd> vanderGBigOPlus(2, Eigen::MatrixXd::Zero(vander.rows(), data.Nkm1));
        vanderGBigOPlus[0] = vander.leftCols(data.Nkm1).array().colwise() * vander.col(2).array();
        vanderGBigOPlus[1] = vander.leftCols(data.Nkm1).array().colwise() * (-1.0 * vander.col(1).array());
        return vanderGBigOPlus;
    };
};
} // namespace Utilities
} // namespace Polydim

#endif
