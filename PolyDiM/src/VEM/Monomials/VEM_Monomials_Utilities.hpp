#ifndef __VEM_Monomials_Monomials_Utilities_HPP
#define __VEM_Monomials_Monomials_Utilities_HPP

#include "VEM_IMonomials.hpp"

namespace Polydim
{
  namespace VEM
  {
    struct VEM_Monomials_Data
    {
        /// \brief The dimension of the polytope.
        /// \note The value can be 1, 2 or 3.
        unsigned int Dimension;
        unsigned int NumMonomials; ///< Number of monomials in the basis.
        std::vector<Eigen::VectorXi> Exponents; ///< Table of exponents of each monomial.
        std::vector<Eigen::MatrixXd> DerivativeMatrices; ///< Matrices used to compute derivatives of monomials.
        Eigen::MatrixXd Laplacian; ///< Matrix used to compute the laplacian of monomials.
    };

    /// Class for VEM_Monomial base implementation
    class VEM_Monomials_Utilities final
    {
      private:
        const VEM_IMonomials& vemMonomials;

      public:
        VEM_Monomials_Utilities(const VEM_IMonomials& vemMonomials); ///< Class constructor.
        ~VEM_Monomials_Utilities() {} ///< Class destructor.


        Eigen::MatrixXd Vander(const std::vector<Eigen::VectorXd>& points,
                               const Eigen::VectorXd& centroid,
                               const double& diam) const;
        Eigen::MatrixXd Vander(const Eigen::MatrixXd& points,
                               const Eigen::Vector3d& centroid,
                               const double& diam) const;
        std::vector<Eigen::MatrixXd> VanderDerivatives(const Eigen::MatrixXd& Vander,
                                                       const double& diam) const;
        Eigen::MatrixXd VanderLaplacian(const Eigen::MatrixXd& Vander,
                                        const double& diam) const;
    };
  }
}

#endif
