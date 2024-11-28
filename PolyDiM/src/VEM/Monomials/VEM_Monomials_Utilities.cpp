#include "VEM_Monomials_Utilities.hpp"

using namespace std;
using namespace Eigen;

namespace Polydim
{
  namespace VEM
  {
    //****************************************************************************
    VEM_Monomials_Utilities::VEM_Monomials_Utilities(const VEM_IMonomials& vemMonomials) :
      vemMonomials(vemMonomials)
    {
    }
    // ***************************************************************************
    MatrixXd VEM_Monomials_Utilities::Vander(const vector<VectorXd>& points,
                                             const VectorXd& centroid,
                                             const double& diam) const
    {
      MatrixXd vander;
      const unsigned int dimension = vemMonomials.Dimension();
      const unsigned int numMonomials = vemMonomials.NumMonomials();
      const unsigned int nQ = points.size();
      vector<MatrixXd> VanderPartial;
      if (numMonomials > 1)
      {
        const unsigned int polynomialDegree = vemMonomials.PolynomialDegree();
        const double inverseDiam = 1.0 / diam;
        VanderPartial.resize(dimension);
        for(unsigned int i = 0; i < dimension; i++)
        {
          VanderPartial[i].resize(nQ, polynomialDegree + 1);
          // first column is set to one
          VanderPartial[i].col(0).setConstant(1);
          // second column has (x-xc)/he
          for(unsigned int k = 0; k < nQ; k++)
            VanderPartial[i](k, 1) = (points[k](i) - centroid(i)) * inverseDiam;
          // other columns are computed by multiplication of VanderPartial columns
          for(unsigned int k = 2; k <= polynomialDegree; k++)
            VanderPartial[i].col(k) = VanderPartial[i].col(k-1).cwiseProduct(VanderPartial[i].col(1));
        }
      }

      // Compute Vander using VanderPartial
      vander.resize(nQ, numMonomials);
      vander.col(0).setConstant(1);
      for(unsigned int k = 1; k < numMonomials; k++)
      {
        const VectorXi expo = vemMonomials.Exponents(k);
        vander.col(k) = VanderPartial[0].col(expo(0)).cwiseProduct(VanderPartial[1].col(expo(1)));
        if(dimension == 3)
          vander.col(k) = vander.col(k).cwiseProduct(VanderPartial[2].col(expo(2)));
      }
      return vander;
    }
    //****************************************************************************
    MatrixXd VEM_Monomials_Utilities::Vander(const MatrixXd& points,
                                             const Vector3d& centroid,
                                             const double& diam) const
    {
      MatrixXd vander;
      const unsigned int numPoints = points.cols();
      if(vemMonomials.NumMonomials()>1)
      {
        // VanderPartial[i]'s rows contain (x-x_E)^i/h_E^i,
        // (y-y_E)^i/h_E^i and (possibly) (z-z_E)^i/h_E^i respectively.
        // Size is dimension x numPoints.
        vector<MatrixXd> VanderPartial(vemMonomials.PolynomialDegree()+1, MatrixXd(vemMonomials.Dimension(), numPoints));
        double inverseDiam = 1.0/diam;
        VanderPartial[0].setOnes(vemMonomials.Dimension(), numPoints);
        VanderPartial[1] = (points.colwise() - centroid)*inverseDiam;

        for(unsigned int i = 2; i <= vemMonomials.PolynomialDegree(); i++)
          VanderPartial[i] = VanderPartial[i-1].cwiseProduct(VanderPartial[1]);

        vander.resize(numPoints, vemMonomials.NumMonomials());
        vander.col(0).setOnes();
        for(unsigned int i = 1; i < vemMonomials.NumMonomials(); ++i)
        {
          const VectorXi expo = vemMonomials.Exponents(i);

          vander.col(i) = (VanderPartial[expo[0]].row(0)).transpose();
          if (vemMonomials.Dimension() > 1)
            vander.col(i) = vander.col(i)
                            .cwiseProduct(VanderPartial[expo[1]].row(1).transpose());
          if (vemMonomials.Dimension() > 2)
            vander.col(i) = vander.col(i)
                            .cwiseProduct(VanderPartial[expo[2]].row(2).transpose());
        }
      }
      else
        vander.setOnes(numPoints,1);

      return vander;
    }
    //****************************************************************************
    vector<MatrixXd> VEM_Monomials_Utilities::VanderDerivatives(const MatrixXd& Vander,
                                                                const double& diam) const
    {
      vector<MatrixXd> vanderDerivatives;
      vanderDerivatives.resize(vemMonomials.Dimension());
      for(unsigned int i=0; i< vemMonomials.Dimension(); i++)
      {
        vanderDerivatives[i].resizeLike(Vander);
        vanderDerivatives[i].col(0).setZero();
      }
      if(vemMonomials.NumMonomials()>1)
      {
        double inverseDiam = 1.0/diam;
        for(unsigned int k=1; k < vemMonomials.NumMonomials(); k++)
        {
          vector<int> derIndices = vemMonomials.DerivativeIndices(k);
          for(unsigned int i=0; i<vemMonomials.Dimension(); i++)
          {
            if(derIndices[i]>=0)
              vanderDerivatives[i].col(k) = inverseDiam *
                                            vemMonomials.DerivativeMatrix(i)(k, derIndices[i]) *
                                            Vander.col(derIndices[i]);
            else
              vanderDerivatives[i].col(k).setZero();
          }
        }
      }

      return vanderDerivatives;
    }
    //****************************************************************************
    MatrixXd VEM_Monomials_Utilities::VanderLaplacian(const MatrixXd& vander,
                                                      const double& diam) const
    {
      MatrixXd vanderLaplacian;

      vanderLaplacian.resizeLike(vander);
      vanderLaplacian.block(0, 0, vander.rows(), 3).setZero();
      MatrixXd laplacian = vemMonomials.Lapl();

      if(vemMonomials.NumMonomials()>3)
      {
        const double inverseDiamSqrd = 1.0/(diam*diam);
        for(unsigned int k = 3; k < vemMonomials.NumMonomials(); k++)
        {
          vector<int> secondDerIndices = vemMonomials.SecondDerivativeIndices(k);
          if(secondDerIndices[0]>=0)
            vanderLaplacian.col(k) = inverseDiamSqrd *
                                     laplacian(k, secondDerIndices[0]) *
                vander.col(secondDerIndices[0]);
          else
            vanderLaplacian.col(k).setZero();
          for(unsigned int i=1; i<vemMonomials.Dimension(); i++)
          {
            if(secondDerIndices[i]>=0)
              vanderLaplacian.col(k) += inverseDiamSqrd *
                                        laplacian(k, secondDerIndices[i])*vander.col(secondDerIndices[i]);
          }
        }
      }

      return vanderLaplacian;
    }
    //****************************************************************************
  }
}
