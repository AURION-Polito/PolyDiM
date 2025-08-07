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

#ifndef __I_VEM_PCC_3D_LocalSpace_HPP
#define __I_VEM_PCC_3D_LocalSpace_HPP

#include "Eigen/Eigen"
#include "I_VEM_PCC_2D_ReferenceElement.hpp"
#include "I_VEM_PCC_3D_ReferenceElement.hpp"
#include "VEM_PCC_3D_LocalSpace_Data.hpp"
#include "VEM_PCC_Utilities.hpp"
#include <vector>

namespace Polydim
{
namespace VEM
{
namespace PCC
{
/// \brief Interface class for Primal Conforming Constant degree 3D Virtual Element Methods \cite DassiMascotto2018
/// \cite Teora2024.
class I_VEM_PCC_3D_LocalSpace
{
public:
    virtual ~I_VEM_PCC_3D_LocalSpace()
    {
    }

    /// \brief Compute data of VEM space on a polygon.
    /// \param reference_element_data Data of the reference element of the VEM space.
    /// \param polygon The geometry of the polygon
    /// \return A data structure containing the information about the local VEM space.
    virtual VEM_PCC_3D_LocalSpace_Data CreateLocalSpace(const VEM_PCC_2D_ReferenceElement_Data &reference_element_data_2D,
                                                        const VEM_PCC_3D_ReferenceElement_Data &reference_element_data_3D,
                                                        const std::vector<VEM_PCC_2D_Polygon_Geometry> &polygonalFaces,
                                                        const VEM_PCC_3D_Polyhedron_Geometry &polyhedron) const = 0;
    /// \brief Compute matrix representation of dofi-dofi stabilization.
    /// \param localSpace Data of the local space.
    /// \param projectionType Type of projection operator to be used in the stabilization.
    /// \return The matrix.
    virtual Eigen::MatrixXd ComputeDofiDofiStabilizationMatrix(const VEM_PCC_3D_LocalSpace_Data &localSpace,
                                                               const ProjectionTypes &projectionType) const = 0;

    virtual Eigen::MatrixXd ComputeDRecipeStabilizationMatrix(const VEM_PCC_3D_LocalSpace_Data &localSpace,
                                                              const ProjectionTypes &projectionType,
                                                              const Eigen::MatrixXd &coercivity_matrix,
                                                              const Eigen::VectorXd &vector_coefficients) const = 0;

    /// \brief Compute values of a suitable projection of basis functions at default quadrature points.
    /// \param localSpace Data of the local space.
    /// \param projectionType Type of projection operator to be applied.
    /// \return A matrix containing the values of a basis function's projection on each column.
    virtual Eigen::MatrixXd ComputeBasisFunctionsValues(const VEM_PCC_3D_LocalSpace_Data &localSpace,
                                                        const ProjectionTypes &projectionType) const = 0;

    /// \brief Compute values of a suitable projection of derivatives of basis functions at default quadrature points.
    /// \param localSpace Data of the local space.
    /// \param projectionType Type of projection operator to be applied.
    /// \return A vector of three matrices, each one containing values of projections of derivatives on columns.
    virtual std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const VEM_PCC_3D_LocalSpace_Data &localSpace,
                                                                               const ProjectionTypes &projectionType) const = 0;

    /// \brief Compute values of a suitable projection of laplacian of basis functions at default quadrature points.
    /// \param localSpace Data of the local space.
    /// \param projectionType Type of projection operator to be applied.
    /// \return A matrix containing the values of a basis function's laplacian projection on each column.
    virtual Eigen::MatrixXd ComputeBasisFunctionsLaplacianValues(const VEM_PCC_3D_LocalSpace_Data &localSpace,
                                                                 const ProjectionTypes &projectionType) const = 0;

    /// \brief Compute values of a suitable projection of basis functions at given points.
    /// \param reference_element_data Data of the reference element of the VEM space.
    /// \param localSpace Data of the local space.
    /// \param projectionType Type of projection operator to be applied.
    /// \param points Evaluation points.
    /// \return A matrix containing the values of a basis function's projection on each column.
    virtual Eigen::MatrixXd ComputeBasisFunctionsValues(const VEM_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                        const VEM_PCC_3D_LocalSpace_Data &localSpace,
                                                        const ProjectionTypes &projectionType,
                                                        const Eigen::MatrixXd &points) const = 0;

    /// \brief Compute values of a suitable projection of derivatives of basis functions at given points.
    /// \param reference_element_data Data of the reference element of the VEM space.
    /// \param localSpace Data of the local space.
    /// \param projectionType Type of projection operator to be applied.
    /// \param points Evaluation points.
    /// \return A vector of three matrices, each one containing values of projections of derivatives on columns.
    virtual std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const VEM_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                                               const VEM_PCC_3D_LocalSpace_Data &localSpace,
                                                                               const ProjectionTypes &projectionType,
                                                                               const Eigen::MatrixXd &points) const = 0;
    /// \brief Compute values of a suitable projection of laplacian of basis functions at default quadrature points.
    /// \param reference_element_data Data of the reference element of the VEM space.
    /// \param localSpace  Data of the local space.
    /// \param projectionType Type of projection operator to be applied.
    /// \param points Evaluation points.
    /// \return A matrix containing the values of a basis function's laplacian projection on each column.
    virtual Eigen::MatrixXd ComputeBasisFunctionsLaplacianValues(const VEM_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                                 const VEM_PCC_3D_LocalSpace_Data &localSpace,
                                                                 const ProjectionTypes &projectionType,
                                                                 const Eigen::MatrixXd &points) const = 0;

    /// \brief Compute values of the polynomial basis at default quadrature points.
    /// \param localSpace Data of the local space.
    /// \return A matrix containing the values of a basis function on each column.
    virtual Eigen::MatrixXd ComputePolynomialsValues(const VEM_PCC_3D_LocalSpace_Data &localSpace) const = 0;

    /// \brief Compute values of the polynomial basis at given points.
    /// \param reference_element_data Data of the reference element of the VEM space.
    /// \param localSpace Data of the local space.
    /// \param points Evaluation points.
    /// \return A matrix containing the values of a basis function on each column.
    virtual Eigen::MatrixXd ComputePolynomialsValues(const VEM_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                     const VEM_PCC_3D_LocalSpace_Data &localSpace,
                                                     const Eigen::MatrixXd &points) const = 0;

    /// \brief Compute values of the derivatives of the polynomial basis at default quadrature points.
    /// \param localSpace Data of the local space.
    /// \return A vector of three matrices, each one containing values of derivatives on columns.
    virtual std::vector<Eigen::MatrixXd> ComputePolynomialsDerivativeValues(const VEM_PCC_3D_LocalSpace_Data &localSpace) const = 0;

    /// \brief Compute values of the derivatives of the polynomial basis at given points.
    /// \param reference_element_data Data of the reference element of the VEM space.
    /// \param localSpace Data of the local space.
    /// \param points Evaluation points.
    /// \return A vector of three matrices, each one containing values of derivatives on columns.
    virtual std::vector<Eigen::MatrixXd> ComputePolynomialsDerivativeValues(const VEM_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                                            const VEM_PCC_3D_LocalSpace_Data &localSpace,
                                                                            const Eigen::MatrixXd &points) const = 0;

    /// \brief Compute values of the laplacian of the polynomial basis at given points.
    /// \param reference_element_data Data of the reference element of the VEM space.
    /// \param localSpace Data of the local space.
    /// \param points Evaluation points.
    /// \return A matrix containing the values of the laplacian of a basis function on each column.
    virtual Eigen::MatrixXd ComputePolynomialsLaplacianValues(const VEM_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                              const VEM_PCC_3D_LocalSpace_Data &localSpace,
                                                              const Eigen::MatrixXd &points) const = 0;

    /// \brief Compute values of the trace of VEM basis functions at given points on an edge.
    /// \param localSpace Data of the local space.
    /// \param edgeInternalPoints Points where edge internal degrees of freedom are located.
    /// \param pointsCurvilinearCoordinates Curvilinear coordinates of evaluation points in [0,1].
    /// \return A matrix containing the values of a basis function on each column.    
    virtual Eigen::MatrixXd ComputeValuesOnEdge(const VEM_PCC_3D_ReferenceElement_Data &reference_element_data,
                                                const VEM_PCC_3D_LocalSpace_Data &localSpace,
                                                const Eigen::VectorXd &pointsCurvilinearCoordinates) const = 0;
};
} // namespace PCC
} // namespace VEM
} // namespace Polydim

#endif
