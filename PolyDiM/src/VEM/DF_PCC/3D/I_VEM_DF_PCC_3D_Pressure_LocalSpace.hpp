#ifndef __I_VEM_DF_PCC_3D_Pressure_LocalSpace_HPP
#define __I_VEM_DF_PCC_3D_Pressure_LocalSpace_HPP

#include "Eigen/Eigen"
#include "I_VEM_DF_PCC_3D_ReferenceElement.hpp"
#include "VEM_DF_PCC_3D_LocalSpace_Data.hpp"
#include <vector>

namespace Polydim
{
namespace VEM
{
namespace DF_PCC
{

/// \brief Class used for computing values of basis functions of 2D
/// Divergence Free Primal Conforming Constant degree Virtual Element Methods.
///
/// Please cite the following article:
///     - <a href="https://doi.org/10.1016/j.matcom.2023.10.003">"Improving high-order VEM stability on badly-shaped
///     elements. Stefano Berrone, Gioana Teora and Fabio Vicini. (2024)"</a>

class I_VEM_DF_PCC_3D_Pressure_LocalSpace
{
  public:
    virtual ~I_VEM_DF_PCC_3D_Pressure_LocalSpace()
    {
    }

    virtual VEM_DF_PCC_3D_Pressure_LocalSpace_Data CreateLocalSpace(const VEM_DF_PCC_3D_Pressure_ReferenceElement_Data &reference_element_data,
                                                                    const VEM_DF_PCC_3D_Polyhedron_Geometry &polygon) const = 0;

    /// \brief Compute the values of projections of VEM basis functions at the internal quadrature points.
    /// \param localSpace: an object of type \ref VEM::PCC::VEM_DF_PCC_3D_Pressure_LocalSpace_Data which contains local
    /// matrices. \param projectionType: the \ref VEM::PCC::ProjectionTypes reporting the kind of projector used to
    /// access to the point-wise evalution of VE basis functions. \return A matrix of size numQuadrature \f$\times\f$
    /// numDOFs whose columns contain the evaluation of the projection of each basis function at the internal quadrature
    /// points.
    virtual Eigen::MatrixXd ComputeBasisFunctionsValues(const VEM_DF_PCC_3D_Pressure_LocalSpace_Data &localSpace) const = 0;

    /// \brief Compute the values of projections of VEM basis functions at points.
    /// \param reference_element_data: an object of type \ref VEM::PCC::VEM_DF_PCC_3D_Pressure_ReferenceElement_Data
    /// which contains monomials stuff. \param polygon: an object of type \ref VEM::PCC::VEM_DF_PCC_3D_Polygon_Geometry
    /// which contains the geoemtric properties of the elements. \param localSpace: an object of type \ref
    /// VEM::PCC::VEM_DF_PCC_3D_Pressure_LocalSpace_Data which contains local matrices. \param projectionType: the \ref
    /// VEM::PCC::ProjectionTypes reporting the kind of projector used to access to the point-wise evalution of VE basis
    /// functions. \param points: a matrix 3 \f$\times\f$ numPoints reporting the coordinates of points. \return A
    /// matrix of size numPoints \f$\times\f$ numDOFs whose columns contain the evaluation of the projection of each
    /// basis function at points.
    virtual Eigen::MatrixXd ComputeBasisFunctionsValues(const VEM_DF_PCC_3D_Pressure_ReferenceElement_Data &reference_element_data,
                                                        const VEM_DF_PCC_3D_Pressure_LocalSpace_Data &localSpace,
                                                        const Eigen::MatrixXd &points) const = 0;

    /// \brief Compute the values of monomial basis functions at the internal quadrature points.
    /// \param localSpace: an object of type \ref VEM::PCC::VEM_DF_PCC_3D_Pressure_LocalSpace_Data which contains local
    /// matrices. \return A matrix of size numQuadrature \f$\times\f$ numMonomials whose columns contain the evaluation
    /// of monomials at the internal quadrature points.
    virtual Eigen::MatrixXd ComputePolynomialsValues(const VEM_DF_PCC_3D_Pressure_LocalSpace_Data &localSpace) const = 0;

    /// \brief Compute the values of monomial basis functions at points.
    /// \param reference_element_data: an object of type \ref VEM::PCC::VEM_DF_PCC_3D_Pressure_ReferenceElement_Data
    /// which contains monomials stuff. \param polygon: an object of type \ref VEM::PCC::VEM_DF_PCC_3D_Polygon_Geometry
    /// which contains the geoemtric properties of the elements. \param points: a matrix 3 \f$\times\f$ numPoints
    /// reporting the coordinates of points. \return A matrix of size numPoints \f$\times\f$ numMonomials whose columns
    /// contain the evaluation of monomials at points.
    virtual Eigen::MatrixXd ComputePolynomialsValues(const VEM_DF_PCC_3D_Pressure_ReferenceElement_Data &reference_element_data,
                                                     const VEM_DF_PCC_3D_Pressure_LocalSpace_Data &localSpace,
                                                     const Eigen::MatrixXd &points) const = 0;
};
} // namespace DF_PCC
} // namespace VEM
} // namespace Polydim

#endif
