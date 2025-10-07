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

#ifndef __VEM_DF_PCC_Utilities_HPP
#define __VEM_DF_PCC_Utilities_HPP

#include "Eigen/Eigen"
#include "Gedim_Macro.hpp"
#include "Monomials_Data.hpp"
#include "lagrange_1D.hpp"

namespace Polydim
{
namespace VEM
{
namespace DF_PCC
{

enum struct ProjectionTypes
{
    Pi0km2 = 0,
    Pi0k = 1,
    PiNabla = 2,
    Pi0km1Der = 3
};

struct VEM_DF_PCC_Utilities final
{
    Eigen::VectorXd ComputeEdgeBasisCoefficients(const unsigned int &order, const Eigen::VectorXd &edgeInternalPoints) const
    {
        // Compute basis function coefficients on the generic edge.
        Eigen::VectorXd interpolation_points_x(order + 1);
        interpolation_points_x << 0.0, 1.0, edgeInternalPoints;
        return Polydim::Interpolation::Lagrange::Lagrange_1D_coefficients(interpolation_points_x);
    }

    std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const unsigned int dimension,
                                                                       const Polydim::VEM::DF_PCC::ProjectionTypes &projectionType,
                                                                       const unsigned int &Nkm1,
                                                                       const Eigen::MatrixXd &vanderInternal,
                                                                       const std::vector<Eigen::MatrixXd> &vanderInternalDerivatives,
                                                                       const std::vector<Eigen::MatrixXd> &piNabla,
                                                                       const std::vector<Eigen::MatrixXd> &pi0km1Der) const
    {
        switch (projectionType)
        {
        case Polydim::VEM::DF_PCC::ProjectionTypes::PiNabla: {
            std::vector<Eigen::MatrixXd> basisFunctionDerivativeValues(dimension * dimension);

            for (unsigned short j = 0; j < dimension; ++j)
                for (unsigned short i = 0; i < dimension; ++i)
                    basisFunctionDerivativeValues[dimension * j + i] = vanderInternalDerivatives[i] * piNabla[j];

            return basisFunctionDerivativeValues;
        }
        case Polydim::VEM::DF_PCC::ProjectionTypes::Pi0km1Der: {
            std::vector<Eigen::MatrixXd> basisFunctionDerivativeValues(dimension * dimension);

            for (unsigned short j = 0; j < dimension; ++j)
                for (unsigned short i = 0; i < dimension; ++i)
                    basisFunctionDerivativeValues[dimension * j + i] =
                        vanderInternal.leftCols(Nkm1) * pi0km1Der[dimension * j + i];

            return basisFunctionDerivativeValues;
        }
        default:
            throw std::runtime_error("Unknown projector type");
        }
    }

    inline Eigen::MatrixXd ComputeBasisFunctionsDivergenceValues(const unsigned int &Nkm1,
                                                                 const Eigen::MatrixXd &vanderInternal,
                                                                 const Eigen::MatrixXd &vmatrix) const
    {
        return vanderInternal.leftCols(Nkm1) * vmatrix;
    }

    inline std::vector<Eigen::MatrixXd> ComputeBasisFunctionsValues(const unsigned int dimension,
                                                                    const Polydim::VEM::DF_PCC::ProjectionTypes &projectionType,
                                                                    const unsigned int &Nkm2,
                                                                    const std::vector<Eigen::MatrixXd> &pi0km2,
                                                                    const std::vector<Eigen::MatrixXd> &pi0k,
                                                                    const Eigen::MatrixXd &vanderInternal) const
    {
        std::vector<Eigen::MatrixXd> basisFunctionValues(dimension);
        switch (projectionType)
        {
        case Polydim::VEM::DF_PCC::ProjectionTypes::Pi0km2:
            for (unsigned short i = 0; i < dimension; ++i)
                basisFunctionValues[i] = vanderInternal.leftCols(Nkm2) * pi0km2[i];
            break;
        case Polydim::VEM::DF_PCC::ProjectionTypes::Pi0k:
            for (unsigned short i = 0; i < dimension; ++i)
                basisFunctionValues[i] = vanderInternal * pi0k[i];
            break;
        default:
            throw std::runtime_error("Unknown projector type");
        }
        return basisFunctionValues;
    }

#if PYBIND == 1
    template <typename MonomialType>
    inline Eigen::MatrixXd ComputePolynomialsValues(const Eigen::MatrixXd &vanderInternal, const MonomialType &monomials) const
    {
        return vanderInternal;
    }
#else
    inline Eigen::MatrixXd ComputePolynomialsValues(const Eigen::MatrixXd &vanderInternal) const
    {
        return vanderInternal;
    }
#endif

    template <typename MonomialType>
    inline Eigen::MatrixXd ComputePolynomialsValues(const Polydim::Utilities::Monomials_Data &data,
                                                    const MonomialType &monomials,
                                                    const Eigen::Vector3d &centroid,
                                                    const double &diameter,
                                                    const Eigen::MatrixXd &points) const
    {
        return monomials.Vander(data, points, centroid, diameter);
    }

#if PYBIND == 1
    template <typename MonomialType>
    inline std::vector<Eigen::MatrixXd> ComputePolynomialsDerivativeValues(const std::vector<Eigen::MatrixXd> &vanderInternalDerivatives,
                                                                           const MonomialType &monomials) const
    {
        return vanderInternalDerivatives;
    }
#else
    inline std::vector<Eigen::MatrixXd> ComputePolynomialsDerivativeValues(const std::vector<Eigen::MatrixXd> &vanderInternalDerivatives) const
    {
        return vanderInternalDerivatives;
    }
#endif

    template <typename MonomialType>
    inline std::vector<Eigen::MatrixXd> ComputePolynomialsDerivativeValues(const Polydim::Utilities::Monomials_Data &data,
                                                                           const MonomialType &monomials,
                                                                           const double &diameter,
                                                                           const Eigen::MatrixXd &vander) const
    {
        return monomials.VanderDerivatives(data, vander, diameter);
    }

    template <typename MonomialType>
    inline Eigen::MatrixXd ComputePolynomialsLaplacianValues(const Polydim::Utilities::Monomials_Data &data,
                                                             const MonomialType &monomials,
                                                             const double &diameter,
                                                             const Eigen::MatrixXd &vander) const
    {
        return monomials.VanderLaplacian(data, vander, diameter);
    }

    Eigen::MatrixXd ComputeValuesOnEdge(const Eigen::RowVectorXd &edgeInternalPoints,
                                        const unsigned int &order,
                                        const Eigen::VectorXd &edgeBasisCoefficients,
                                        const Eigen::VectorXd &pointsCurvilinearCoordinates) const
    {
        Eigen::VectorXd interpolation_points_x(order + 1);
        interpolation_points_x << 0.0, 1.0, edgeInternalPoints.transpose();
        return Polydim::Interpolation::Lagrange::Lagrange_1D_values(interpolation_points_x, edgeBasisCoefficients, pointsCurvilinearCoordinates);
    }

    Eigen::MatrixXd ComputeDofiDofiStabilizationMatrix(const unsigned int dimension,
                                                       const std::vector<Eigen::MatrixXd> &projector,
                                                       const double &coefficient,
                                                       const std::vector<Eigen::MatrixXd> &dmatrix) const
    {
        Eigen::MatrixXd staBmatrix = dmatrix[0] * projector[0];

        for (unsigned int d = 1; d < dimension; d++)
            staBmatrix += dmatrix[d] * projector[d];

        staBmatrix.diagonal().array() -= 1;

        // staBmatrix = (\Pi^{0,dofs}_order - I)^T * (\Pi^{0,dofs}_order - I).
        staBmatrix = coefficient * staBmatrix.transpose() * staBmatrix;

        return staBmatrix;
    }
};
} // namespace DF_PCC
} // namespace VEM
} // namespace Polydim

#endif
