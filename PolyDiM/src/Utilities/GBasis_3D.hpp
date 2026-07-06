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

#ifndef __GBasis_3D_HPP
#define __GBasis_3D_HPP

#include "GBasis_Data.hpp"
#include "Monomials_3D.hpp"

namespace Polydim
{
namespace Utilities
{
class GBasis_3D final
{
  private:
    /// @brief 3D scalar monomial basis \f$\{m_\alpha\}\f$ used to build the vector basis.
    Polydim::Utilities::Monomials_3D monomials;

    /// @brief Compute the decomposition indices of a single vector monomial onto the G-basis.
    ///
    /// Given the multi-index @p expo of a scalar monomial, returns, for each of the
    /// three Cartesian components, the four column indices addressing the internal
    /// storage blocks of the decomposition: the \f$\mathcal{G}^\nabla\f$ block
    /// (index 0) and the three \f$\mathcal{G}^\oplus\f$ sub-blocks (indices 1–3,
    /// corresponding respectively to the first group of size @c DimFirstBasis and to
    /// the two \f$\mathbb{P}_{k-1}\f$ blocks). Each entry is an `Eigen::Vector4i`.
    ///
    /// @param data Precomputed G-basis data (dimensions, exponent map, vector decomposition).
    /// @param expo Exponent multi-index \f$\alpha = (\alpha_1, \alpha_2, \alpha_3)\f$ of the scalar monomial.
    /// @return One `Eigen::Vector4i` per Cartesian component, holding the four block indices.
    /// @note Internal helper invoked by Compute(); not part of the public interface.
    std::vector<Eigen::Vector4i> VectorDecompositionIndices(const Polydim::Utilities::GBasis_Data &data,
                                                            const Eigen::VectorXi &expo) const;

  public:
    /// @brief Build the vector polynomial G-basis of a given degree in 3D.
    ///
    /// Assembles all quantities needed to represent the vector polynomial space
    /// \f$[\mathbb{P}_k]^3\f$ through its decomposition
    /// \f$\mathcal{G}_k^\nabla \oplus \mathcal{G}_k^\oplus\f$, where
    /// \f$\mathcal{G}_k^\nabla = \nabla \mathbb{P}_{k+1}\f$ and
    /// \f$\mathcal{G}_k^\oplus = \boldsymbol{x} \times [\mathbb{P}_{k-1}]^3\f$.
    /// The \f$\mathcal{G}^\oplus\f$ block is stored in three parts of widths
    /// @c DimFirstBasis, @c Nkm1 and @c Nkm1 (total \f$\dim\mathcal{G}_k^\oplus = \f$
    /// @c DimFirstBasis @c + @c 2*Nkm1).
    ///
    /// @param polynomial_degree Maximum polynomial degree \f$k\f$ of the basis.
    /// @return A fully populated Polydim::Utilities::GBasis_Data structure (space
    ///         dimensions, exponent matrix, and the four-block vector-decomposition maps).
    Polydim::Utilities::GBasis_Data Compute(const unsigned int polynomial_degree);

    /// @brief Evaluate the Vandermonde matrix of the \f$\mathcal{G}_k^\oplus\f$ block.
    ///
    /// Builds the three Cartesian components (x, y, z) of the Vandermonde matrix of
    /// \f$\mathcal{G}_k^\oplus = \boldsymbol{x} \times [\mathbb{P}_{k-1}]^3\f$,
    /// evaluated at the points whose scalar monomial values are provided in @p vander.
    ///
    /// @param data   Precomputed G-basis data (@c DimFirstBasis and @c Nkm1 are used here).
    /// @param vander Vandermonde matrix of the scalar monomial basis evaluated at the
    ///               evaluation points (one row per point), from which the linear
    ///               monomials \f$x_1, x_2, x_3\f$ and the \f$\mathbb{P}_{k-1}\f$
    ///               monomials are extracted.
    /// @return A vector of three `Eigen::MatrixXd` (x-, y- and z-components), each of size
    ///         @c vander.rows() \f$\times\f$ (@c data.DimFirstBasis @c + @c 2*data.Nkm1),
    ///         forming the Vandermonde matrix of \f$\mathcal{G}_k^\oplus\f$.
    std::vector<Eigen::MatrixXd> VanderGBigOPlus(const Polydim::Utilities::GBasis_Data &data, const Eigen::MatrixXd &vander) const;

    /// @brief Assemble the vector decomposition into a two-block (\f$\nabla\f$, \f$\oplus\f$) form.
    ///
    /// Repackages the internal four-block storage of @p data into a compact
    /// per-component representation. For each of the three Cartesian components,
    /// the result holds two matrices: the gradient block
    /// \f$\mathcal{G}_k^\nabla\f$ (size \f$N_k \times N_{k+1}\f$) and the full
    /// \f$\mathcal{G}_k^\oplus\f$ block, obtained by horizontally concatenating the
    /// three stored sub-blocks into a single \f$N_k \times (\text{DimFirstBasis} + 2N_{k-1})\f$ matrix.
    ///
    /// @param data Precomputed G-basis data containing the four-block @c VectorDecomposition.
    /// @return A `std::vector` of three components, each a pair of `Eigen::MatrixXd`:
    ///         index 0 is the \f$\nabla\f$ block, index 1 is the assembled \f$\oplus\f$ block.
    std::vector<std::vector<Eigen::MatrixXd>> VectorDecomposition(const Polydim::Utilities::GBasis_Data &data) const
    {
        std::vector<std::vector<Eigen::MatrixXd>> result(3);
        for (unsigned int i = 0; i < data.Dimension; i++)
        {
            result[i].resize(2);
            result[i][0] = data.VectorDecomposition[i][0];
            result[i][1].setZero(data.Nk, 2 * data.Nkm1 + data.DimFirstBasis);
            result[i][1] << data.VectorDecomposition[i][1], data.VectorDecomposition[i][2], data.VectorDecomposition[i][3];
        }
        return result;
    }
};
} // namespace Utilities
} // namespace Polydim

#endif
