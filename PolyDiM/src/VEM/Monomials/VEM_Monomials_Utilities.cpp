#include "VEM_Monomials_Utilities.hpp"

using namespace std;
using namespace Eigen;

namespace Polydim
{
namespace VEM
{
namespace Monomials
{
template struct VEM_Monomials_Utilities<1>;
template struct VEM_Monomials_Utilities<2>;
template struct VEM_Monomials_Utilities<3>;
//****************************************************************************
template <unsigned short dimension>
MatrixXi VEM_Monomials_Utilities<dimension>::Exponents(const VEM_Monomials_Data &data) const
{
    MatrixXi exponents(dimension, data.NumMonomials);

    for (unsigned int m = 0; m < data.NumMonomials; m++)
        exponents.col(m) << data.Exponents[m];

    return exponents;
}
//****************************************************************************
template <unsigned short dimension>
MatrixXd VEM_Monomials_Utilities<dimension>::Vander(const VEM_Monomials_Data &data,
                                                    const MatrixXd &points,
                                                    const Vector3d &centroid,
                                                    const double &diam) const
{
    MatrixXd vander;
    const unsigned int numPoints = points.cols();
    if (data.NumMonomials > 1)
    {
        // VanderPartial[i]'s rows contain (x-x_E)^i/h_E^i,
        // (y-y_E)^i/h_E^i and (possibly) (z-z_E)^i/h_E^i respectively.
        // Size is dimension x numPoints.
        vector<MatrixXd> VanderPartial(data.PolynomialDegree + 1, MatrixXd(dimension, numPoints));
        double inverseDiam = 1.0 / diam;
        VanderPartial[0].setOnes(dimension, numPoints);
        VanderPartial[1] = (points.colwise() - centroid) * inverseDiam;

        for (unsigned int i = 2; i <= data.PolynomialDegree; i++)
            VanderPartial[i] = VanderPartial[i - 1].cwiseProduct(VanderPartial[1]);

        vander.resize(numPoints, data.NumMonomials);
        vander.col(0).setOnes();
        for (unsigned int i = 1; i < data.NumMonomials; ++i)
        {
            const VectorXi expo = data.Exponents[i];

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
//****************************************************************************
} // namespace Monomials
} // namespace VEM
} // namespace Polydim
