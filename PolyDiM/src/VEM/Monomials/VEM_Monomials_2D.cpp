#include "VEM_Monomials_2D.hpp"

using namespace std;
using namespace Eigen;

namespace Polydim
{
  namespace VEM
  {
    //****************************************************************************
    VEM_Monomials_2D::VEM_Monomials_2D() :
      utilities(*this)
    {
      data.Dimension = 2;

      data.NumMonomials = (config.PolynomialDegree() + 1) * (config.PolynomialDegree() + 2) * 0.5;
      data.Exponents.resize(data.NumMonomials) ;
      data.Exponents[0].setZero(data.Dimension);
      data.DerivativeMatrices.resize(data.Dimension);
      for(unsigned int i=0; i< data.Dimension; i++)
        data.DerivativeMatrices[i].setZero(data.NumMonomials, data.NumMonomials);
      data.Laplacian.setZero(data.NumMonomials, data.NumMonomials);
      for (unsigned int i=1; i < data.NumMonomials; i++)
      {
        data.Exponents[i].resize(data.Dimension);
        if (data.Exponents[i-1](0)==0)
          data.Exponents[i] << data.Exponents[i-1](1)+1 , 0;
        else
          data.Exponents[i] << data.Exponents[i-1](0)-1 , data.Exponents[i-1](1)+1;
        vector<int> derIndices = DerivativeIndices(i);
        vector<int> secondDerIndices = SecondDerivativeIndices(i);
        for(unsigned int j = 0; j < data.Dimension; ++j)
        {
          if(derIndices[j]>=0)
            data.DerivativeMatrices[j](i,derIndices[j]) = data.Exponents[i](j);
          if(secondDerIndices[j]>=0)
            data.Laplacian(i,secondDerIndices[j]) = data.Exponents[i](j)*(data.Exponents[i](j)-1);
        }
      }
    }
    //****************************************************************************
    MatrixXi VEM_Monomials_2D::Exponents() const
    {
      MatrixXi exponents(2, data.NumMonomials);

      for (unsigned int m = 0; m < data.NumMonomials; m++)
        exponents.col(m)<< data.Exponents[m];

      return exponents;
    }
    //****************************************************************************
    int VEM_Monomials_2D::Index( const VectorXi& _exponents ) const
    {
      int out = -1;
      if(_exponents[0]>=0 && _exponents[1]>=0)
      {
        unsigned int sumExponents = _exponents.sum();
        out = sumExponents*(sumExponents+1)/2+_exponents(1);
      }
      return out;
    }
    //****************************************************************************
    vector<int> VEM_Monomials_2D::DerivativeIndices(const unsigned int& index) const
    {
      vector<int> derivativeIndices;

      const VectorXi& expo = data.Exponents[index];
      const int xDerIndex = index - expo.sum();
      derivativeIndices.resize(2);
      derivativeIndices[0] = (expo(0)>0) ? xDerIndex : -1;
      derivativeIndices[1] = (expo(1)>0) ? xDerIndex-1 : -1;

      return derivativeIndices;
    }
    //****************************************************************************
    vector<int> VEM_Monomials_2D::SecondDerivativeIndices(const unsigned int& index) const
    {
      vector<int> secondDerivativeIndices;
      const VectorXi& expo = data.Exponents[index];
      const int xxDerIndex = index + 1 - 2 * data.Exponents[index].sum();

      secondDerivativeIndices.resize(2);
      secondDerivativeIndices[0] = (expo(0)>1) ? xxDerIndex : -1;
      secondDerivativeIndices[1] = (expo(1)>1) ? xxDerIndex-2 : -1;
      return secondDerivativeIndices;
    }
    //****************************************************************************
  }
}
