#ifndef __I_VEM_DF_PCC_2D_Velocity_LocalSpace_HPP
#define __I_VEM_DF_PCC_2D_Velocity_LocalSpace_HPP

#include "Eigen/Eigen"
#include "I_VEM_DF_PCC_2D_ReferenceElement.hpp"
#include "VEM_DF_PCC_2D_LocalSpace_Data.hpp"
#include "VEM_DF_PCC_Utilities.hpp"
#include <vector>

namespace Polydim
{
namespace VEM
{
namespace DF_PCC
{

/// \brief Class used for computing values of basis functions of 2D
/// Divergence Free Primal Conforming Constant degree Virtual Element Methods.
class I_VEM_DF_PCC_2D_Velocity_LocalSpace
{
  public:
    /// \brief Class destructor.
    virtual ~I_VEM_DF_PCC_2D_Velocity_LocalSpace()
    {
    }

    /// \brief Compute data of VEM space on a polygon.
    /// \param reference_element_data Data of the reference element of the VEM space.
    /// \param polygon The geometry of the polygon.
    /// \return A data structure containing the information about the local space.
    virtual VEM_DF_PCC_2D_Velocity_LocalSpace_Data CreateLocalSpace(const VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                                                    const VEM_DF_PCC_2D_Polygon_Geometry &polygon) const = 0;

    /// \brief Compute matrix representation of dofi-dofi stabilization.
    /// \param localSpace Data of the local space.
    /// \param projectionType Type of projection operator to be used in the stabilization.
    /// \return The matrix.
    virtual Eigen::MatrixXd ComputeDofiDofiStabilizationMatrix(const VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace,
                                                               const ProjectionTypes &projectionType) const = 0;

    /// \brief Compute values of the basis functions of each component at given points on an edge.
    /// \param reference_element_data Data of the reference element of the VEM space.
    /// \param pointsCurvilinearCoordinates Curvilinear coordinates of evaluation points in [0,1].
    /// \return A matrix containing the values of a basis function on each column.
    virtual Eigen::MatrixXd ComputeValuesOnEdge(const VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                                const Eigen::VectorXd &pointsCurvilinearCoordinates) const = 0;

    /// \brief Compute values of a suitable projection of basis functions at default quadrature points.
    /// \param localSpace Data of the local space.
    /// \param projectionType Type of projection operator to be applied.
    /// \return A vector of two matrices containing the values of the two components of a basis function's projection on
    /// each column.
    virtual std::vector<Eigen::MatrixXd> ComputeBasisFunctionsValues(const VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace,
                                                                     const ProjectionTypes &projectionType) const = 0;

    /// \brief Compute values of projected derivatives of each component of a basis function at default quadrature points.
    /// \param localSpace  Data of the local space.
    /// \param projectionType Type of projection operator to be applied.
    /// \return A vector with four components such that the 2 * j + i component is the i-th derivative of the j-th component.
    virtual std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace,
                                                                               const ProjectionTypes &projectionType) const = 0;

    /// \brief Compute values of a suitable projection of basis functions at default quadrature points.
    /// \param reference_element_data Data of the reference element of the VEM space.
    /// \param polygon The geometry of the polygon.
    /// \param localSpace Data of the local space.
    /// \param projectionType Type of projection operator to be applied.
    /// \param points Evaluation points.
    /// \return A vector of two matrices containing the values of the two components of a basis function's projection on
    /// each column.
    virtual std::vector<Eigen::MatrixXd> ComputeBasisFunctionsValues(const VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                                                     const VEM_DF_PCC_2D_Polygon_Geometry &polygon,
                                                                     const VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace,
                                                                     const ProjectionTypes &projectionType,
                                                                     const Eigen::MatrixXd &points) const = 0;

    /// \brief Compute values of projected derivatives of each component of a basis function at default quadrature points.
    /// \param reference_element_data Data of the reference element of the VEM space.
    /// \param polygon The geometry of the polygon.
    /// \param localSpace Data of the local space.
    /// \param projectionType Type of projection operator to be applied.
    /// \param points Evaluation points.
    /// \return A vector of four components such that the 2 * j + i component is the i-th derivative of the j-th component.
    virtual std::vector<Eigen::MatrixXd> ComputeBasisFunctionsDerivativeValues(const VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                                                               const VEM_DF_PCC_2D_Polygon_Geometry &polygon,
                                                                               const VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace,
                                                                               const ProjectionTypes &projectionType,
                                                                               const Eigen::MatrixXd &points) const = 0;

    /// \brief Compute values of the polynomial basis at default quadrature points.
    /// \param localSpace Data of the local space.
    /// \return A matrix containing the values of a basis function on each column.
    virtual Eigen::MatrixXd ComputePolynomialsValues(const VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace) const = 0;

    /// \brief Compute values of the polynomial basis at given points.
    /// \param reference_element_data Data of the reference element of the VEM space.
    /// \param polygon The geometry of the polygon.
    /// \param points Evaluation points.
    /// \return A matrix containing the values of a basis function on each column.
    virtual Eigen::MatrixXd ComputePolynomialsValues(const VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                                     const VEM_DF_PCC_2D_Polygon_Geometry &polygon,
                                                     const Eigen::MatrixXd &points) const = 0;

    /// \brief Compute values of the derivatives of the polynomial basis at default quadrature points.
    /// \param localSpace Data of the local space.
    /// \return A vector of two matrices, each one containing values of derivatives on columns.
    virtual std::vector<Eigen::MatrixXd> ComputePolynomialsDerivativeValues(const VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace) const = 0;

    /// \brief Compute values of the derivatives of the polynomial basis at given points.
    /// \param reference_element_data Data of the reference element of the VEM space.
    /// \param polygon The geometry of the polygon.
    /// \param points Evaluation points.
    /// \return A vector of two matrices, each one containing values of derivatives on columns.
    virtual std::vector<Eigen::MatrixXd> ComputePolynomialsDerivativeValues(const VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                                                            const VEM_DF_PCC_2D_Polygon_Geometry &polygon,
                                                                            const Eigen::MatrixXd &points) const = 0;

    /// \brief Compute values of the laplacian of the polynomial basis at given points.
    /// \param reference_element_data Data of the reference element of the VEM space.
    /// \param polygon The geometry of the polygon.
    /// \param points Evaluation points.
    /// \return A matrix containing the values of the laplacian of a basis function on each column.
    virtual Eigen::MatrixXd ComputePolynomialsLaplacianValues(const VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &reference_element_data,
                                                              const VEM_DF_PCC_2D_Polygon_Geometry &polygon,
                                                              const Eigen::MatrixXd &points) const = 0;

    /// \brief Compute values of the divergence of basis functions at default quadrature points.
    /// \param localSpace Data of the local space.
    /// \return A matrix containing the values of the divergence of basis function on each column.
    virtual Eigen::MatrixXd ComputeBasisFunctionsDivergenceValues(const VEM_DF_PCC_2D_Velocity_LocalSpace_Data &localSpace) const = 0;
};
} // namespace DF_PCC
} // namespace VEM
} // namespace Polydim

#endif
