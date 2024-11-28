#ifndef __VEM_Monomials_VEM_Monomials_2D_HPP
#define __VEM_Monomials_VEM_Monomials_2D_HPP

#include "VEM_Monomials_Utilities.hpp"

namespace Polydim
{
  namespace VEM
  {
    class VEM_Monomials_2D final : public VEM_IMonomials
    {
      private:
        VEM_Monomials_Data data;
        VEM_Monomials_Utilities utilities;

      public:
        VEM_Monomials_2D();
        virtual ~VEM_Monomials_2D() {}

        /// \brief Get the dimension of the polynomial basis.
        /// \returns A const reference to \ref dimension.
        inline unsigned int Dimension() const { return data.Dimension; }
        /// \brief Get the polynomial degree of the basis.
        /// \returns A const reference to \ref polynomialDegree.
        inline unsigned int PolynomialDegree() const { return config.PolynomialDegree(); }
        /// \brief Get the number of monomials.
        /// \returns A const reference to \ref numMonomials.
        inline unsigned int NumMonomials() const { return data.NumMonomials; }
        /// \brief Get the table of the exponents of the monomial basis.
        /// \returns A const reference to \ref exponents.
        Eigen::VectorXi Exponents( const int& _index ) const { return data.Exponents[_index]; }
        Eigen::MatrixXi Exponents() const;
        /// \param The required derivative.
        /// \returns A const reference to \ref derivativeMatrices[i].
        Eigen::MatrixXd DerivativeMatrix(const unsigned int& i) const { return data.DerivativeMatrices[i]; }
        /// \returns A const reference to the x-derivative matrix.
        Eigen::MatrixXd D_x() const { return data.DerivativeMatrices[0]; }
        /// \returns A const reference to the y-derivative matrix.
        Eigen::MatrixXd D_y() const { return data.DerivativeMatrices[1]; }
        /// \returns A const reference to the z-derivative matrix.
        Eigen::MatrixXd D_z() const { return data.DerivativeMatrices[2]; }
        /// \returns A const reference to the laplacian matrix.
        Eigen::MatrixXd Lapl() const { return data.Laplacian; }

        int Index( const Eigen::VectorXi& _exponents ) const;
        vector<int> DerivativeIndices(const unsigned int& index) const;
        vector<int> SecondDerivativeIndices(const unsigned int& index) const;
        Eigen::MatrixXd Vander(const vector<Eigen::VectorXd>& points,
                               const Eigen::VectorXd& centroid,
                               const double& diam) const
        {
          return utilities.Vander(points,
                                  centroid,
                                  diam);
        }
        Eigen::MatrixXd Vander(const Eigen::MatrixXd& points,
                               const Eigen::Vector3d& centroid,
                               const double& diam) const
        {
          return utilities.Vander(points,
                                  centroid,
                                  diam);
        }
        vector<Eigen::MatrixXd> VanderDerivatives(const Eigen::MatrixXd& vander,
                                                  const double& diam) const
        {
          return utilities.VanderDerivatives(vander,
                                             diam);
        }
        Eigen::MatrixXd VanderLaplacian(const Eigen::MatrixXd& vander,
                                        const double& diam) const
        {
          return utilities.VanderLaplacian(vander,
                                           diam);
        }
    };
  }

#endif
