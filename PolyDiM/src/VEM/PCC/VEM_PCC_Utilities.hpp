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
#include "Gedim_Macro.hpp"
#include "Monomials_Data.hpp"
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

struct VEM_PCC_Utilities final
{
    Eigen::VectorXd ComputeEdgeBasisCoefficients(const unsigned int &order, const Eigen::VectorXd &edgeInternalPoints) const
    {
        // Compute basis function coefficients on the generic edge.
        Eigen::VectorXd interpolation_points_x(order + 1);
        interpolation_points_x << 0.0, 1.0, edgeInternalPoints;
        return Polydim::Interpolation::Lagrange::Lagrange_1D_coefficients(interpolation_points_x);
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

    inline Eigen::MatrixXd EdgeDOFsCoordinates(const Eigen::RowVectorXd &referenceEdgeDOFsPoint,
                                               const Eigen::MatrixXd &vertices,
                                               const Eigen::MatrixXi &edges,
                                               const std::vector<bool> &edgesDirection,
                                               const Eigen::MatrixXd &edgesTangent,
                                               const unsigned int &edge_local_index) const
    {

        const unsigned int num_edge_dofs = referenceEdgeDOFsPoint.cols();

        if (num_edge_dofs == 0)
            return Eigen::MatrixXd(0, 0);

        const Eigen::Vector3d edge_origin = edgesDirection.at(edge_local_index) ? vertices.col(edges(0, edge_local_index))
                                                                                : vertices.col(edges(1, edge_local_index));

        const Eigen::Vector3d edge_tangent = edgesTangent.col(edge_local_index);
        const double edge_direction = edgesDirection[edge_local_index] ? 1.0 : -1.0;

        Eigen::MatrixXd edge_dofs_coordinates = Eigen::MatrixXd::Zero(3, num_edge_dofs);
        for (unsigned int r = 0; r < num_edge_dofs; r++)
        {
            edge_dofs_coordinates.col(r) << edge_origin + edge_direction * referenceEdgeDOFsPoint(0, r) * edge_tangent;
        }
        return edge_dofs_coordinates;
    }

    std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const unsigned int dimension,
                                                                       const Polydim::VEM::PCC::ProjectionTypes &projectionType,
                                                                       const unsigned int &Nkm1,
                                                                       const Eigen::MatrixXd &vanderInternal,
                                                                       const std::vector<Eigen::MatrixXd> &vanderInternalDerivatives,
                                                                       const Eigen::MatrixXd &piNabla,
                                                                       const std::vector<Eigen::MatrixXd> &pi0km1Der) const
    {
        switch (projectionType)
        {
        case Polydim::VEM::PCC::ProjectionTypes::PiNabla: {
            std::vector<Eigen::MatrixXd> basisFunctionsDerivativeValues(dimension);

            for (unsigned short i = 0; i < dimension; ++i)
                basisFunctionsDerivativeValues[i] = vanderInternalDerivatives[i] * piNabla;

            return basisFunctionsDerivativeValues;
        }
        case Polydim::VEM::PCC::ProjectionTypes::Pi0km1Der: {
            std::vector<Eigen::MatrixXd> basisFunctionsDerivativeValues(dimension);

            for (unsigned short i = 0; i < dimension; ++i)
                basisFunctionsDerivativeValues[i] = vanderInternal.leftCols(Nkm1) * pi0km1Der[i];

            return basisFunctionsDerivativeValues;
        }
        default:
            throw std::runtime_error("Unknown projector type");
        }
    }

    Eigen::MatrixXd ComputeBasisFunctionsLaplacianValues(const unsigned int dimension,
                                                         const Polydim::VEM::PCC::ProjectionTypes &projectionType,
                                                         const unsigned int &Nkm1,
                                                         const std::vector<Eigen::MatrixXd> &vanderInternalDerivatives,
                                                         const std::vector<Eigen::MatrixXd> &pi0km1Der) const
    {
        switch (projectionType)
        {
        case Polydim::VEM::PCC::ProjectionTypes::Pi0km1Der: {
            Eigen::MatrixXd basisFunctionsLaplacianValues = vanderInternalDerivatives[0].leftCols(Nkm1) * pi0km1Der[0];
            for (unsigned int d = 1; d < dimension; ++d)
                basisFunctionsLaplacianValues += vanderInternalDerivatives[d].leftCols(Nkm1) * pi0km1Der[d];

            return basisFunctionsLaplacianValues;
        }
        default:
            throw std::runtime_error("Unknown projector type");
        }
    }

    inline Eigen::MatrixXd ComputeBasisFunctionsValues(const Polydim::VEM::PCC::ProjectionTypes &projectionType,
                                                       const unsigned int &Nkm1,
                                                       const Eigen::MatrixXd &pi0km1,
                                                       const Eigen::MatrixXd &pi0k,
                                                       const Eigen::MatrixXd &vanderInternal) const
    {
        switch (projectionType)
        {
        case Polydim::VEM::PCC::ProjectionTypes::Pi0km1:
            return vanderInternal.leftCols(Nkm1) * pi0km1;
        case Polydim::VEM::PCC::ProjectionTypes::Pi0k:
            return vanderInternal * pi0k;
        default:
            throw std::runtime_error("Unknown projector type");
        }
    }

#if PYBIND == 1
    template <typename MonomialType>
    inline Eigen::MatrixXd ComputePolynomialsValues(const Eigen::MatrixXd &vanderInternal, const MonomialType &) const
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
                                                                           const MonomialType &) const
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

    Eigen::MatrixXd ComputeDRecipeStabilizationMatrix(const Eigen::MatrixXd &projector,
                                                      const Eigen::MatrixXd &coercivity_matrix,
                                                      const Eigen::VectorXd &vector_coefficients,
                                                      const Eigen::MatrixXd &Dmatrix) const
    {
        Eigen::MatrixXd stabMatrix = Dmatrix * projector;
        stabMatrix.diagonal().array() -= 1;

        const Eigen::VectorXd diagonal_coercivity = coercivity_matrix.diagonal();
        Eigen::MatrixXd max_matrix = Eigen::MatrixXd::Zero(coercivity_matrix.cols(), 2);
        max_matrix << diagonal_coercivity, vector_coefficients;

        const Eigen::VectorXd weights = max_matrix.rowwise().maxCoeff();

        // stabMatrix = (\Pi^{\nabla,dofs}_order - I)^T * (\Pi^{\nabla,dofs}_order - I).
        stabMatrix = stabMatrix.transpose() * weights.asDiagonal() * stabMatrix;

        return stabMatrix;
    }
};
} // namespace PCC
} // namespace VEM
} // namespace Polydim

#endif
