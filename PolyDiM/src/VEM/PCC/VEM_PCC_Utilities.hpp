#ifndef __VEM_PCC_Utilities_HPP
#define __VEM_PCC_Utilities_HPP

#include "Eigen/Eigen"

namespace Polydim
{
  namespace VEM
  {
    namespace PCC
    {
      /// \brief Base class for computing values of basis functions of Primal Conforming Constant degree
      /// Virtual Element Methods.
      /// \copyright See top level LICENSE file for details.
      template<unsigned short dimension>
      class VEM_ValuesUtilities_PCC
      {
        private:
        public:
          VEM_ValuesUtilities_PCC();
          virtual ~VEM_ValuesUtilities_PCC() {}

          /// Compute the Edge basis coefficients
          Eigen::VectorXd ComputeEdgeBasisCoefficients(const unsigned int& order,
                                                       const Eigen::VectorXd& edgeInternalPoints) const;

          /// \brief Compute matrix \ref stabMatrix with PiNabla projector.
          /// \note This requires \ref ComputePiNabla() to be called previously.
          /// \return MainApplication::Output::Success if the computation was successful.
          Eigen::MatrixXd ComputeStabilizationMatrix(const Eigen::MatrixXd& piNabla,
                                                     const double& diameter,
                                                     const Eigen::MatrixXd& DMatrix) const;
          /// \brief Compute matrix \ref stabMatrix with Pi0k projector.
          /// \note This requires \ref ComputeL2Projectors() to be called previously.
          /// \return MainApplication::Output::Success if the computation was successful.
          Eigen::MatrixXd ComputeStabilizationMatrixPi0k(const Eigen::MatrixXd& pi0k,
                                                         const double& measure,
                                                         const Eigen::MatrixXd& DMatrix) const;
          /// \brief Compute matrices \ref pi0km1 and \ref pi0k.
          /// \return MainApplication::Output::Success if the computation was successful.
          /// \note This requires \ref ComputePiNabla() to be called previously.
          void ComputeL2Projectors(const unsigned int& numberProjectorBasisFunctions,
                                   const unsigned int& numberBasisFunctions,
                                   const unsigned int& numberInternalBasisFunctions,
                                   const unsigned int& order,
                                   const Eigen::MatrixXd& piNabla,
                                   const Eigen::MatrixXd& Hmatrix,
                                   const double& measure,
                                   const unsigned int& Nkm1,
                                   const Eigen::LLT<Eigen::MatrixXd>& H_km1_LLT,
                                   Eigen::MatrixXd& pi0km1,
                                   Eigen::MatrixXd& pi0k) const;

          /// \brief Compute the values of the polynomial projection of
          /// derivatives of basis functions at internal quadrature points on
          /// the geometry.
          /// \param basisFunctionsDerivativeValues The vector of matrices of
          /// values. Its length equals \ref Dimension(). Each column of each
          /// matrix will contain the values of a basis function's projected
          /// derivative at internal quadrature points.
          std::vector<Eigen::MatrixXd> ComputeBasisFunctionDerivativeValues(const unsigned int& Nkm1,
                                                                            const Eigen::MatrixXd& vanderInternal,
                                                                            const std::vector<Eigen::MatrixXd>& pi0km1Der) const
          {
            std::vector<Eigen::MatrixXd> basisFunctionsDerivativeValues;

            basisFunctionsDerivativeValues.resize(dimension);
            for(unsigned short i = 0; i < dimension; ++i)
              basisFunctionsDerivativeValues[i] = vanderInternal.leftCols(Nkm1) *
                                                  pi0km1Der[i];

            return basisFunctionsDerivativeValues;
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
          Eigen::MatrixXd ComputeBasisFunctionLaplacianValues(const unsigned int& Nkm1,
                                                              const std::vector<Eigen::MatrixXd>& vanderInternalDerivatives,
                                                              const std::vector<Eigen::MatrixXd>& pi0km1Der) const
          {
            Eigen::MatrixXd basisFunctionLaplacianValues;

            basisFunctionLaplacianValues = vanderInternalDerivatives[0].leftCols(Nkm1) * pi0km1Der[0];
            for(unsigned int d = 1; d < dimension; ++d)
              basisFunctionLaplacianValues += vanderInternalDerivatives[d].leftCols(Nkm1) *
                                              pi0km1Der[d];

            return basisFunctionLaplacianValues;
          }

          /// \brief Compute internal quadrature weights on the geometry used to compute projectors.
          /// \details These weights correspond to the quadrature points used in the methods to get
          /// values that do not accept points as inputs, that is:
          /// - \ref ComputeProjectedBasisFunctionValues() [1/2]
          /// - \ref ComputeProjectedBasisFunctionDerivativeValues() [1/2]
          /// - \ref ComputePolynomialBasisValues() [1/2]
          /// - \ref ComputePolynomialBasisDerivativeValues() [1/2]
          /// - \ref ComputePolynomialBasisLaplacianValues() [1/2]
          /// \param weightsVector The vector to be filled.
          /// \sa \ref ComputeInternalQuadraturePoints().
          inline Eigen::VectorXd ComputeInternalQuadratureWeights(const Eigen::VectorXd& internalWeights) const
          {
            return internalWeights;
          }

          /// \brief Compute internal quadrature points on the geometry used to compute projectors.
          /// \details These points are the quadrature points used in the methods to get values that do
          /// not accept points as inputs, that is:
          /// - \ref ComputeProjectedBasisFunctionValues() [1/2]
          /// - \ref ComputeProjectedBasisFunctionDerivativeValues() [1/2]
          /// - \ref ComputePolynomialBasisValues() [1/2]
          /// - \ref ComputePolynomialBasisDerivativeValues() [1/2]
          /// - \ref ComputePolynomialBasisLaplacianValues() [1/2]
          /// \param pointsMatrix The matrix to be filled. The size will be \ref Dimension() x \ref
          /// NumberInternalQuadraturePoints().
          /// \sa \ref ComputeInternalQuadratureWeights().
          inline Eigen::MatrixXd ComputeInternalQuadraturePoints(const Eigen::MatrixXd& internalQuadraturePoints) const
          {
            return internalQuadraturePoints;
          }

          /// \brief Compute the values of the polynomial projection of basis functions at internal
          /// quadrature points on the geometry.
          /// \details To control the projection used, call either \ref UsePi0km1() or \ref UsePi0k().
          /// \param basisFunctionValues The matrix of values. Each column will contain the values of a
          /// basis function's projection at internal quadrature points.
          inline Eigen::MatrixXd ComputeBasisFunctionValues(const bool& usePi0km1,
                                                            const unsigned int& Nkm1,
                                                            const Eigen::MatrixXd& pi0km1,
                                                            const Eigen::MatrixXd& pi0k,
                                                            const Eigen::MatrixXd& vanderInternal) const
          {
            if (usePi0km1)
              return vanderInternal.leftCols(Nkm1)*pi0km1;
            else
              return vanderInternal*pi0k;
          }

          /// \brief Compute the values of the polynomial projection of basis functions at given points
          /// on the geometry.
          /// \details To control the projection used, call either \ref UsePi0km1() or \ref UsePi0k().
          /// \param points The points at which to evaluate projections.
          /// \param basisFunctionValues The matrix of values. Each column will contain the values of a
          /// basis function's projection at the given points.
          /// \return MainApplication::Output::Success if the computation was successful.
          Eigen::MatrixXd ComputeBasisFunctionValues(const bool& usePi0km1,
                                                     const unsigned int& Nkm1,
                                                     const Eigen::MatrixXd& pi0km1,
                                                     const Eigen::MatrixXd& pi0k,
                                                     const Eigen::MatrixXd& vanderInternal,
                                                     const Eigen::MatrixXd& polynomialBasisValues) const;

          /// \brief Compute the values of the polynomial projection of derivatives of basis functions
          /// at given points on the geometry.
          /// \param points The points at which to evaluate projections.
          /// \param basisFunctionDerivativeValues The vector of matrices of values. Its length equals
          /// \ref Dimension(). Each column of each matrix will contain the values of a basis function's
          /// projected derivative at the given points.
          std::vector<Eigen::MatrixXd> ComputeBasisFunctionDerivativeValues(const unsigned int& Nkm1,
                                                                            const std::vector<Eigen::MatrixXd>& pi0km1Der,
                                                                            const Eigen::MatrixXd& polynomialBasisValues) const;

          /// \brief Compute the values of the polynomial projection of the laplacian of basis functions
          /// at given points on the geometry.
          /// \details The projection of the laplacian considered is \f$ \widetilde{\Delta}\varphi =
          /// \nabla \cdot \left(\Pi^0_{\mathrm{order}-1} \nabla \varphi \right)\f$.
          /// \param points The points at which to evaluate projections.
          /// \param basisFunctionValues The matrix of values. Each column will contain the values of a
          /// basis function's projected laplacian at the given points.
          Eigen::MatrixXd ComputeBasisFunctionLaplacianValuesOnPoints(const unsigned int& Nkm1,
                                                                      const std::vector<Eigen::MatrixXd>& pi0km1Der,
                                                                      const std::vector<Eigen::MatrixXd>& polynomialBasisDerivativeValues) const;

          /// \brief Compute the values of the basis functions of the polynomial basis of projectors at
          /// internal quadrature points on the geometry.
          /// \param polynomialBasisValues The matrix of values. Each column will contain the values of
          /// a basis function at internal quadrature points.
          inline Eigen::MatrixXd ComputePolynomialBasisValues(const Eigen::MatrixXd& vanderInternal) const
          {
            return vanderInternal;
          }

          /// \brief Compute the values of the basis functions of the polynomial basis of projectors at
          /// given points on the geometry.
          /// \param points The points at which to evaluate the basis functions.
          /// \param polynomialBasisValues The matrix of values. Each column will contain the values of
          /// a basis function at the given points.
          inline Eigen::MatrixXd ComputePolynomialBasisValues(const VEM_IMonomials& monomials,
                                                              const Eigen::Vector3d& centroid,
                                                              const double& diameter,
                                                              const Eigen::MatrixXd& points) const
          {
            return monomials.Vander(points, centroid, diameter);
          }

          /// \brief Compute the values of the derivatives of the basis functions of the polynomial
          /// basis of projectors at internal quadrature points on the geometry.
          /// \param polynomialBasisDerivativeValues The matrix of values. Each column will contain the
          /// values of a basis function's derivative at internal quadrature points.
          inline std::vector<Eigen::MatrixXd> ComputePolynomialBasisDerivativeValues(const std::vector<Eigen::MatrixXd>& vanderInternalDerivatives) const
          {
            return vanderInternalDerivatives;
          }

          /// \brief Compute the values of the derivatives of the basis functions of the polynomial
          /// basis of projectors at given points on the geometry.
          /// \param points The points at which to evaluate the derivatives.
          /// \param polynomialBasisDerivativeValues The vector of matrices of values. Its length equals
          /// \ref Dimension(). Each column of each matrix will contain the values of a basis function's
          /// derivative at the given points.
          inline std::vector<Eigen::MatrixXd> ComputePolynomialBasisDerivativeValues(const VEM_IMonomials& monomials,
                                                                                     const double& diameter,
                                                                                     const Eigen::MatrixXd& vander) const
          {
            return monomials.VanderDerivatives(vander,
                                               diameter);
          }

          /// \brief Compute the values of the laplacian of the basis functions of the polynomial basis
          /// of projectors at given points.
          /// \param points The points at which to evaluate the laplacian.
          /// \param polynomialBasisLaplacianValues The matrix of values. Each column will contain the
          /// values of a basis function's laplacian at the given points.
          inline Eigen::MatrixXd ComputePolynomialBasisLaplacianValues(const VEM_IMonomials& monomials,
                                                                       const double& diameter,
                                                                       const Eigen::MatrixXd& vander) const
          {
            return monomials.VanderLaplacian(vander,
                                             diameter);
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
          Eigen::MatrixXd ComputeValuesOnEdge(const Eigen::RowVectorXd& edgeInternalPoints,
                                              const unsigned int& order,
                                              const Eigen::VectorXd& edgeBasisCoefficients,
                                              const Eigen::VectorXd& pointsCurvilinearCoordinates) const;

      };
    }
  }
}

#endif
