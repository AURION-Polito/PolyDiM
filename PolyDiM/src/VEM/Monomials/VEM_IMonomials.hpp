#ifndef __VEM_Monomials_IMonomials_HPP
#define __VEM_Monomials_IMonomials_HPP

#include <vector>
#include "Eigen/Eigen"

namespace Polydim
{
  namespace VEM
  {
    struct VEM_Monomials_Data
    {
        unsigned int Order; ///< monomial space order
        unsigned int Dimension; ///< The geometric dimension
        unsigned int NumMonomials; ///< Number of monomials in the basis.
        std::vector<Eigen::VectorXi> Exponents; ///< Table of exponents of each monomial.
        std::vector<Eigen::MatrixXd> DerivativeMatrices; ///< Matrices used to compute derivatives of monomials.
        Eigen::MatrixXd Laplacian; ///< Matrix used to compute the laplacian of monomials.
    };

    ///  \brief Class used to manage scaled monomials in VEM classes.
    ///  \details These scaled monomials are used as basis for local
    ///  polynomial projections in VEM classes. The generic element of the
    ///  basis on a given polytope \f$E\f$ of dimension 2 or 3 has the
    ///  following analytical expression: \f[ m(\mathbf{x}) =
    ///  \frac{\prod_{i=1}^{d}(x_i -
    ///  x_{E,i})^{\alpha_i}}{h_E^{\sum_{i=1}^d\alpha_i}} \f] where \f$d\f$
    ///  is the dimension of the polytope, \f$\mathbf{x}_E\f$ is a
    ///  reference point on \f$E\f$ (usually the centroid) and \f$h_E\f$ is
    ///  the diameter of \f$E\f$.
    ///  \copyright See top level LICENSE file for details.
    class VEM_IMonomials
    {
      public:
        virtual ~VEM_IMonomials() {}

        /// \brief Compute the monomials data
        /// \param order the polynomial space order
        virtual VEM_Monomials_Data Compute(const unsigned int order) const = 0;
        /// \brief Get the dimension of the polynomial basis.
        /// \returns A const reference to \ref dimension.
        virtual unsigned int Dimension() const = 0;
        /// \brief Get the polynomial degree of the basis.
        /// \returns A const reference to \ref polynomialDegree.
        virtual unsigned int PolynomialDegree() const = 0;
        /// \brief Get the number of monomials.
        /// \returns A const reference to \ref numMonomials.
        virtual unsigned int NumMonomials() const = 0;
        /// \brief Get the table of the exponents of the monomial basis.
        /// \returns A const reference to \ref exponents.
        virtual Eigen::VectorXi Exponents( const int& _index ) const = 0;
        virtual Eigen::MatrixXi Exponents() const = 0;
        /// \param The required derivative.
        /// \returns A const reference to \ref derivativeMatrices[i].
        virtual Eigen::MatrixXd DerivativeMatrix(const unsigned int& i) const = 0;
        /// \returns A const reference to the x-derivative matrix.
        virtual Eigen::MatrixXd D_x() const = 0;
        /// \returns A const reference to the y-derivative matrix.
        virtual Eigen::MatrixXd D_y() const = 0;
        /// \returns A const reference to the z-derivative matrix.
        virtual Eigen::MatrixXd D_z() const = 0;
        /// \returns A const reference to the laplacian matrix.
        virtual Eigen::MatrixXd Lapl() const = 0;

        /// \brief Get the index of the monomial with input exponents.
        /// \param _exponents The exponents of the requested monomial.
        /// \returns The index of the monomial, or -1 if the exponents are not valid.
        virtual int Index( const Eigen::VectorXi& _exponents) const = 0;
        /// \brief Compute the indices of the two monomials that are multiples
        /// of the derivatives of a monomial.
        /// \details If \ref dimension is 3 and the input monomial is
        /// \f[m(\mathbf{x}) =
        /// \frac{(\mathbf{x}-\mathbf{x}_E)^{\mathbf{\alpha}}}{h_E^{|\alpha|}},
        /// \alpha = (\alpha_x,\alpha_y,\alpha_z)\,, \f] the method
        /// returns the indices of \f$m_x\f$, \f$m_y\f$ and \f$m_z\f$ such
        /// that \f[\frac{\partial m}{\partial x} = \frac{\alpha_x}{h_E}
        /// m_x \,, \frac{\partial m}{\partial y} = \frac{\alpha_y}{h_E}
        /// m_y \,, \frac{\partial m}{\partial z} = \frac{\alpha_z}{h_E}
        /// m_z \,.\f] If \ref dimension is 2, only the indices of \f$m_x\f$
        /// and \f$m_y\f$ are returned. One or more indices can be -1,
        /// indicating that the corresponding derivative is 0.
        /// \param index The monomial to be derived.
        /// \param derivativeIndices the resulting indices
        /// \returns A vector of
        /// size \ref dimension containing the indices of \f$m_x\f$,
        /// \f$m_y\f$ and (if \ref dimension is 3), \f$m_z\f$.
        /// \sa \ref SecondDerivativeIndices().
        virtual std::vector<int> DerivativeIndices(const unsigned int& index) const = 0;
        /// \brief Compute the indices of the two monomials that are multiples
        /// of the second derivatives of a monomial.
        /// \details If \ref dimension is 3 and the input monomial is
        /// \f[m(\mathbf{x}) =
        /// \frac{(\mathbf{x}-\mathbf{x}_E)^{\mathbf{\alpha}}}{h_E^{|\alpha|}},
        /// \alpha = (\alpha_x,\alpha_y,\alpha_z)\,, \f] the method
        /// returns the indices of \f$m_x\f$, \f$m_y\f$ and \f$m_z\f$ such
        /// that \f[\frac{\partial^2 m}{\partial x^2} =
        /// \frac{\alpha_x(\alpha_x-1)}{h_E^2} m_x\,, \frac{\partial^2
        /// m}{\partial y^2} = \frac{\alpha_y(\alpha_y-1)}{h_E^2} m_y\,,
        /// \frac{\partial^2 m}{\partial z^2} =
        /// \frac{\alpha_z(\alpha_z-1)}{h_E^2} m_z \,. \f] If \ref
        /// dimension is 2, only the indices of \f$m_x\f$ and \f$m_y\f$
        /// are returned. One or more indices can be -1, indicating that
        /// the corresponding derivative is 0.
        /// \param index The monomial to be derived.
        /// \param secondDerivativeIndices the resulting indices
        /// \returns A vector of size \ref dimension containing the
        /// indices of \f$m_x\f$, \f$m_y\f$ and (if \ref dimension is 3),
        /// \f$m_z\f$.
        /// \sa \ref DerivativeIndices().
        virtual std::vector<int> SecondDerivativeIndices(const unsigned int& index) const = 0;

        /// \brief Compute the Vandermonde matrix of the monomial basis in
        /// given points.
        /// \details Each column of the matrix contains the values of a
        /// monomial at the input points.
        /// \param points The vector of coordinates of the points.
        /// \param centroid The reference point \f$\mathbf{x}_E\f$ (see Detailed description).
        /// \param diam The diameter of the polytope (see Detailed description).
        /// \param Vander The matrix to be filled.
        virtual Eigen::MatrixXd Vander(const std::vector<Eigen::VectorXd>& points,
                                       const Eigen::VectorXd& centroid,
                                       const double& diam) const = 0;

        /// \brief Compute the Vandermonde matrix of the monomial basis in
        /// given points.
        /// \details Each column of the matrix contains the values of a
        /// monomial at the input points.
        /// \param points The matrix of coordinates of the points. Each
        /// column contains the coordinates of a point.
        /// \param centroid The reference point \f$\mathbf{x}_E\f$ (see
        /// Detailed description).
        /// \param diam The diameter of the polytope (see Detailed
        /// description).
        /// \param Vander The matrix to be filled.
        virtual Eigen::MatrixXd Vander(const Eigen::MatrixXd& points,
                                       const Eigen::Vector3d& centroid,
                                       const double& diam) const = 0;

        /// \brief Compute the Vandermonde matrices of the derivatives of
        /// the monomial basis in given points.
        /// \details Each column of the output matrices contains the
        /// values of a monomial derivative at the input points.
        /// \param Vander The Vandermonde matrix of monomials at given
        /// points.
        /// \param diam The diameter of the polytope (see Detailed
        /// description).
        /// \param vanderDerivatives The vector of matrices to be filled.
        /// \sa \ref Vander().
        virtual std::vector<Eigen::MatrixXd> VanderDerivatives(const Eigen::MatrixXd& Vander,
                                                               const double& diam) const = 0;

        /// \brief Compute the Vandermonde matrix of the laplacian of the
        /// monomial basis in given points.
        /// \details Each column of the matrix contains the values of a
        /// monomial's laplacian at the input points.
        /// \param Vander The Vandermonde matrix of monomials at given
        /// points.
        /// \param diam The diameter of the polytope (see Detailed
        /// description).
        /// \param vanderLaplacian The matrix to be filled.
        /// \sa \ref Vander().
        virtual Eigen::MatrixXd VanderLaplacian(const Eigen::MatrixXd& Vander,
                                                const double& diam) const = 0;
    };
  }
}

#endif
