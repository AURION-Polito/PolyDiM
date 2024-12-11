#ifndef __VEM_GBasis_VEM_GBasis_3D_HPP
#define __VEM_GBasis_VEM_GBasis_3D_HPP

#include "VEM_GBasis_Data.hpp"
#include "VEM_Monomials_3D.hpp"

namespace Polydim
{
namespace VEM
{
namespace Monomials
{
class VEM_GBasis_3D final
{
private:
    VEM_Monomials_3D monomials;

    std::vector<Eigen::Vector4i> VectorDecompositionIndices(const VEM_GBasis_Data& data,
                                                            const Eigen::VectorXi& expo) const;

public:
    VEM_GBasis_Data Compute(const unsigned int polynomial_degree);

    std::vector<Eigen::MatrixXd> VanderGBigOPlus(const VEM_GBasis_Data& data,
                                                 const Eigen::MatrixXd &vander) const;

    std::vector<std::vector<Eigen::MatrixXd>> VectorDecomposition(const VEM_GBasis_Data& data) const
    {
        std::vector<std::vector<Eigen::MatrixXd>> result(3);
        for(unsigned int i = 0; i < data.Dimension; i++)
        {
            result[i].resize(2);
            result[i][0] = data.VectorDecomposition[i][0];
            result[i][1].setZero(data.Nk, 2 * data.Nkm1 + data.DimFirstBasis);
            result[i][1] << data.VectorDecomposition[i][1], data.VectorDecomposition[i][2], data.VectorDecomposition[i][3];
        }
        return result;
    }

};
}
}
}

#endif
