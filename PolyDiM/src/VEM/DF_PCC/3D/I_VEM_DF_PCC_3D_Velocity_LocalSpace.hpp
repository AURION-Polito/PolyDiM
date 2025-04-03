#ifndef __I_VEM_DF_PCC_3D_Velocity_LocalSpace_HPP
#define __I_VEM_DF_PCC_3D_Velocity_LocalSpace_HPP

#include "Eigen/Eigen"
#include "I_VEM_DF_PCC_3D_ReferenceElement.hpp"
#include "VEM_DF_PCC_3D_LocalSpace_Data.hpp"
#include "VEM_DF_PCC_Utilities.hpp"
#include "VEM_Monomials_3D.hpp"
#include "VEM_PCC_2D_LocalSpace.hpp"
#include "VEM_PCC_2D_ReferenceElement.hpp"
#include <vector>

namespace Polydim
{
namespace VEM
{
namespace DF_PCC
{

/// \brief Class used for computing values of basis functions related to velocitye space of 3D
/// Divergence-Free Primal Conforming Constant degree Virtual Element Methods \cite DaveigaDassi2020.
class I_VEM_DF_PCC_3D_Velocity_LocalSpace
{
  public:
    /// \brief Compute data of VEM space on a polyhedron.
    /// \param reference_element_data_2D Data of the reference element of the VEM space on faces.
    /// \param reference_element_data_3D Data of the reference element of the VEM space.
    /// \param polygonalFaces Vector of geometries of polyhedron faces.
    /// \param polyhedron The geometry of the polyhedron
    /// \return A data structure containing the information about the local space.
    virtual VEM_DF_PCC_3D_Velocity_LocalSpace_Data CreateLocalSpace(const PCC::VEM_PCC_2D_ReferenceElement_Data &reference_element_data_2D,
                                                                    const VEM_DF_PCC_3D_Velocity_ReferenceElement_Data &reference_element_data_3D,
                                                                    const std::vector<PCC::VEM_PCC_2D_Polygon_Geometry> &polygonalFaces,
                                                                    const VEM_DF_PCC_3D_Polyhedron_Geometry &polyhedron) const = 0;

    /// \brief Compute matrix representation of dofi-dofi stabilization.
    /// \param localSpace Data of the local space.
    /// \param projectionType Type of projection operator to be used in the stabilization.
    /// \return The matrix.
    virtual Eigen::MatrixXd ComputeDofiDofiStabilizationMatrix(const VEM_DF_PCC_3D_Velocity_LocalSpace_Data &localSpace,
                                                               const ProjectionTypes &projectionType) const = 0;

    /// \brief Compute values of a suitable projection of basis functions at default quadrature points.
    /// \param localSpace Data of the local space.
    /// \param projectionType Type of projection operator to be applied.
    /// \return A vector of three matrices containing the values of the three components of a basis function's projection on
    /// each column.
    virtual std::vector<Eigen::MatrixXd> ComputeBasisFunctionsValues(const VEM_DF_PCC_3D_Velocity_LocalSpace_Data &localSpace,
                                                                     const ProjectionTypes &projectionType) const = 0;

    /// \brief Compute values of projected derivatives of each component of a basis function at default quadrature points.
    /// \param localSpace  Data of the local space.
    /// \param projectionType Type of projection operator to be applied.
    /// \return A vector with nine components such that the 2 * j + i component is the i-th derivative of the j-th component.
    virtual std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const VEM_DF_PCC_3D_Velocity_LocalSpace_Data &localSpace,
                                                                               const ProjectionTypes &projectionType) const = 0;

    /// \brief Compute values of a suitable projection of basis functions at default quadrature points.
    /// \param reference_element_data Data of the reference element of the VEM space.
    /// \param polyhedron The geometry of the polyhedron.
    /// \param localSpace Data of the local space.
    /// \param projectionType Type of projection operator to be applied.
    /// \param points Evaluation points.
    /// \return A vector of three matrices containing the values of the three components of a basis function's projection on
    /// each column.
    virtual std::vector<Eigen::MatrixXd> ComputeBasisFunctionsValues(const VEM_DF_PCC_3D_Velocity_ReferenceElement_Data &reference_element_data,
                                                                     const VEM_DF_PCC_3D_Polyhedron_Geometry &polyhedron,
                                                                     const VEM_DF_PCC_3D_Velocity_LocalSpace_Data &localSpace,
                                                                     const ProjectionTypes &projectionType,
                                                                     const Eigen::MatrixXd &points) const = 0;

    /// \brief Compute values of projected derivatives of each component of a basis function at default quadrature points.
    /// \param reference_element_data Data of the reference element of the VEM space.
    /// \param polyhedron The geometry of the polyhedron.
    /// \param localSpace Data of the local space.
    /// \param projectionType Type of projection operator to be applied.
    /// \param points Evaluation points.
    /// \return A vector of nine components such that the 2 * j + i component is the i-th derivative of the j-th component.
    virtual std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const VEM_DF_PCC_3D_Velocity_ReferenceElement_Data &reference_element_data,
                                                                               const VEM_DF_PCC_3D_Polyhedron_Geometry &polyhedron,
                                                                               const VEM_DF_PCC_3D_Velocity_LocalSpace_Data &localSpace,
                                                                               const ProjectionTypes &projectionType,
                                                                               const Eigen::MatrixXd &points) const = 0;

    /// \brief Compute values of the polynomial basis at default quadrature points.
    /// \param localSpace Data of the local space.
    /// \return A matrix containing the values of a basis function on each column.
    virtual Eigen::MatrixXd ComputePolynomialsValues(const VEM_DF_PCC_3D_Velocity_LocalSpace_Data &localSpace) const = 0;

    /// \brief Compute values of the polynomial basis at given points.
    /// \param reference_element_data Data of the reference element of the VEM space.
    /// \param polyhedron The geometry of the polyhedron.
    /// \param points Evaluation points.
    /// \return A matrix containing the values of a basis function on each column.
    virtual Eigen::MatrixXd ComputePolynomialsValues(const VEM_DF_PCC_3D_Velocity_ReferenceElement_Data &reference_element_data,
                                                     const VEM_DF_PCC_3D_Polyhedron_Geometry &polyhedron,
                                                     const Eigen::MatrixXd &points) const = 0;

    /// \brief Compute values of the derivatives of the polynomial basis at default quadrature points.
    /// \param localSpace Data of the local space.
    /// \return A vector of three matrices, each one containing values of derivatives on columns.
    virtual std::vector<Eigen::MatrixXd> ComputePolynomialsDerivativeValues(const VEM_DF_PCC_3D_Velocity_LocalSpace_Data &localSpace) const = 0;

    /// \brief Compute values of the derivatives of the polynomial basis at given points.
    /// \param reference_element_data Data of the reference element of the VEM space.
    /// \param polyhedron The geometry of the polyhedron.
    /// \param points Evaluation points.
    /// \return A vector of three matrices, each one containing values of derivatives on columns.
    virtual std::vector<Eigen::MatrixXd> ComputePolynomialsDerivativeValues(const VEM_DF_PCC_3D_Velocity_ReferenceElement_Data &reference_element_data,
                                                                            const VEM_DF_PCC_3D_Polyhedron_Geometry &polyhedron,
                                                                            const Eigen::MatrixXd &points) const = 0;

    /// \brief Compute values of the laplacian of the polynomial basis at given points.
    /// \param reference_element_data Data of the reference element of the VEM space.
    /// \param polyhedron The geometry of the polyhedron.
    /// \param points Evaluation points.
    /// \return A matrix containing the values of the laplacian of a basis function on each column.
    virtual Eigen::MatrixXd ComputePolynomialsLaplacianValues(const VEM_DF_PCC_3D_Velocity_ReferenceElement_Data &reference_element_data,
                                                              const VEM_DF_PCC_3D_Polyhedron_Geometry &polyhedron,
                                                              const Eigen::MatrixXd &points) const = 0;

    /// \brief Compute values of the divergence of basis functions at default quadrature points.
    /// \param localSpace Data of the local space.
    /// \return A matrix containing the values of the divergence of basis function on each column.
    virtual Eigen::MatrixXd ComputeBasisFunctionsDivergenceValues(const VEM_DF_PCC_3D_Velocity_LocalSpace_Data &localSpace) const = 0;
};
} // namespace DF_PCC
} // namespace VEM
} // namespace Polydim

#endif
