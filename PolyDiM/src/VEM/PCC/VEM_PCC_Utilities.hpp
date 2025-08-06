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

#ifndef __VEM_PCC_Utilities_HPP
#define __VEM_PCC_Utilities_HPP

#include "Eigen/Eigen"
#include "VEM_Monomials_Data.hpp"
#include "lagrange_1D.hpp"

namespace Polydim
{
namespace VEM
{
namespace PCC
{

enum struct ProjectionTypes
{
    Pi0km1 = 0,
    Pi0k = 1,
    PiNabla = 2,
    Pi0km1Der = 3,
    Pi0klm1 = 4
};

template <unsigned short dimension> struct VEM_PCC_Utilities final
{
    Eigen::VectorXd ComputeEdgeBasisCoefficients(const unsigned int &order, const Eigen::VectorXd &edgeInternalPoints) const
    {
        // Compute basis function coefficients on the generic edge.
        Eigen::VectorXd interpolation_points_x(order + 1);
        interpolation_points_x << 0.0, 1.0, edgeInternalPoints;
        return Interpolation::Lagrange::Lagrange_1D_coefficients(interpolation_points_x);
    }

    void ComputeL2Projectors(const double &measure,
                             const unsigned int &order,
                             const unsigned int &Nkm1,
                             const unsigned int &Nk,
                             const unsigned int &NumInternalBasisFunctions,
                             const unsigned int &NumBasisFunctions,
                             const Eigen::MatrixXd &Hmatrix,
                             const Eigen::MatrixXd &PiNabla,
                             Eigen::MatrixXd &Cmatrix,
                             Eigen::MatrixXd &Pi0km1,
                             Eigen::MatrixXd &Pi0k) const
    {
        Cmatrix = Eigen::MatrixXd::Zero(Nk, NumBasisFunctions);
        // \int_E \Pi^\nabla_order \phi_j · m_i for m_i of degree > order-2 (enhancement property).
        Cmatrix.bottomRows(Nk - NumInternalBasisFunctions) = Hmatrix.bottomRows(Nk - NumInternalBasisFunctions) * PiNabla;

        if (order > 1)
        {
            Cmatrix.topLeftCorner(NumInternalBasisFunctions, NumBasisFunctions - NumInternalBasisFunctions).setZero();
            //\int_E \phi_j · m_i = measure*\delta_{ij} for m_i of degree <= order-2 (internal dofs).
            Cmatrix.topRightCorner(NumInternalBasisFunctions, NumInternalBasisFunctions) =
                measure * Eigen::MatrixXd::Identity(NumInternalBasisFunctions, NumInternalBasisFunctions);
        }

        Pi0km1 = Hmatrix.topLeftCorner(Nkm1, Nkm1).llt().solve(Cmatrix.topRows(Nkm1));
        Pi0k = Hmatrix.llt().solve(Cmatrix);
    }

    std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const ProjectionTypes &projectionType,
                                                                       const unsigned int &Nkm1,
                                                                       const Eigen::MatrixXd &vanderInternal,
                                                                       const std::vector<Eigen::MatrixXd> &vanderInternalDerivatives,
                                                                       const Eigen::MatrixXd &piNabla,
                                                                       const std::vector<Eigen::MatrixXd> &pi0km1Der) const
    {
        switch (projectionType)
        {
        case ProjectionTypes::PiNabla: {
            std::vector<Eigen::MatrixXd> basisFunctionsDerivativeValues(dimension);

            for (unsigned short i = 0; i < dimension; ++i)
                basisFunctionsDerivativeValues[i] = vanderInternalDerivatives[i] * piNabla;

            return basisFunctionsDerivativeValues;
        }
        case ProjectionTypes::Pi0km1Der: {
            std::vector<Eigen::MatrixXd> basisFunctionsDerivativeValues(dimension);

            for (unsigned short i = 0; i < dimension; ++i)
                basisFunctionsDerivativeValues[i] = vanderInternal.leftCols(Nkm1) * pi0km1Der[i];

            return basisFunctionsDerivativeValues;
        }
        default:
            throw std::runtime_error("Unknown projector type");
        }
    }

    Eigen::MatrixXd ComputeBasisFunctionsLaplacianValues(const ProjectionTypes &projectionType,
                                                         const unsigned int &Nkm1,
                                                         const std::vector<Eigen::MatrixXd> &vanderInternalDerivatives,
                                                         const std::vector<Eigen::MatrixXd> &pi0km1Der) const
    {
        switch (projectionType)
        {
        case ProjectionTypes::Pi0km1Der: {
            Eigen::MatrixXd basisFunctionsLaplacianValues = vanderInternalDerivatives[0].leftCols(Nkm1) * pi0km1Der[0];
            for (unsigned int d = 1; d < dimension; ++d)
                basisFunctionsLaplacianValues += vanderInternalDerivatives[d].leftCols(Nkm1) * pi0km1Der[d];

            return basisFunctionsLaplacianValues;
        }
        default:
            throw std::runtime_error("Unknown projector type");
        }
    }

    inline Eigen::MatrixXd ComputeBasisFunctionsValues(const ProjectionTypes &projectionType,
                                                       const unsigned int &Nkm1,
                                                       const Eigen::MatrixXd &pi0km1,
                                                       const Eigen::MatrixXd &pi0k,
                                                       const Eigen::MatrixXd &vanderInternal) const
    {
        switch (projectionType)
        {
        case ProjectionTypes::Pi0km1:
            return vanderInternal.leftCols(Nkm1) * pi0km1;
        case ProjectionTypes::Pi0k:
            return vanderInternal * pi0k;
        default:
            throw std::runtime_error("Unknown projector type");
        }
    }

    inline Eigen::MatrixXd ComputePolynomialsValues(const Eigen::MatrixXd &vanderInternal) const
    {
        return vanderInternal;
    }

    template <typename VEM_MonomialType>
    inline Eigen::MatrixXd ComputePolynomialsValues(const Utilities::VEM_Monomials_Data &data,
                                                    const VEM_MonomialType &monomials,
                                                    const Eigen::Vector3d &centroid,
                                                    const double &diameter,
                                                    const Eigen::MatrixXd &points) const
    {
        return monomials.Vander(data, points, centroid, diameter);
    }

    inline std::vector<Eigen::MatrixXd> ComputePolynomialsDerivativeValues(const std::vector<Eigen::MatrixXd> &vanderInternalDerivatives) const
    {
        return vanderInternalDerivatives;
    }

    template <typename VEM_MonomialType>
    inline std::vector<Eigen::MatrixXd> ComputePolynomialsDerivativeValues(const Utilities::VEM_Monomials_Data &data,
                                                                           const VEM_MonomialType &monomials,
                                                                           const double &diameter,
                                                                           const Eigen::MatrixXd &vander) const
    {
        return monomials.VanderDerivatives(data, vander, diameter);
    }

    template <typename VEM_MonomialType>
    inline Eigen::MatrixXd ComputePolynomialsLaplacianValues(const Utilities::VEM_Monomials_Data &data,
                                                             const VEM_MonomialType &monomials,
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
        return Interpolation::Lagrange::Lagrange_1D_values(interpolation_points_x, edgeBasisCoefficients, pointsCurvilinearCoordinates);
    }

    Eigen::MatrixXd ComputeDofiDofiStabilizationMatrix(const Eigen::MatrixXd &projector,
                                                       const double &coefficient,
                                                       const Eigen::MatrixXd &Dmatrix) const
    {
        Eigen::MatrixXd stabMatrix = Dmatrix * projector;
        stabMatrix.diagonal().array() -= 1;
        // stabMatrix = (\Pi^{\nabla,dofs}_order - I)^T * (\Pi^{\nabla,dofs}_order - I).
        stabMatrix = coefficient * stabMatrix.transpose() * stabMatrix;

        return stabMatrix;
    }
};
} // namespace PCC
} // namespace VEM
} // namespace Polydim

#endif
