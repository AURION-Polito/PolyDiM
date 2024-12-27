#ifndef __VEM_DF_PCC_2D_Presure_LocalSpace_HPP
#define __VEM_DF_PCC_2D_Presure_LocalSpace_HPP

#include "Eigen/Eigen"
#include "I_VEM_DF_PCC_2D_Pressure_LocalSpace.hpp"
#include "I_VEM_DF_PCC_2D_ReferenceElement.hpp"
#include "VEM_DF_PCC_2D_LocalSpace_Data.hpp"
#include "VEM_Monomials_2D.hpp"
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

class VEM_DF_PCC_2D_Pressure_LocalSpace final : public I_VEM_DF_PCC_2D_Pressure_LocalSpace
{
  private:
    Monomials::VEM_Monomials_2D monomials;

    void InitializeProjectorsComputation(const VEM_DF_PCC_2D_Pressure_ReferenceElement_Data &reference_element_data,
                                         const Eigen::Vector3d &polygonCentroid,
                                         const double &polygonDiameter,
                                         const Eigen::MatrixXd &internalQuadraturePoints,
                                         const Eigen::VectorXd &internalQuadratureWeights,
                                         VEM_DF_PCC_2D_Pressure_LocalSpace_Data &localSpace) const;

  public:
    virtual ~VEM_DF_PCC_2D_Pressure_LocalSpace()
    {
    }

    VEM_DF_PCC_2D_Pressure_LocalSpace_Data CreateLocalSpace(const VEM_DF_PCC_2D_Pressure_ReferenceElement_Data &reference_element_data,
                                                            const VEM_DF_PCC_2D_Polygon_Geometry &polygon) const;

    /// \brief Compute the values of projections of VEM basis functions at the internal quadrature points.
    /// \param localSpace: an object of type \ref VEM::PCC::VEM_DF_PCC_2D_Pressure_LocalSpace_Data which contains local
    /// matrices. \param projectionType: the \ref VEM::PCC::ProjectionTypes reporting the kind of projector used to
    /// access to the point-wise evalution of VE basis functions. \return A matrix of size numQuadrature \f$\times\f$
    /// numDOFs whose columns contain the evaluation of the projection of each basis function at the internal quadrature
    /// points.
    inline Eigen::MatrixXd ComputeBasisFunctionsValues(const VEM_DF_PCC_2D_Pressure_LocalSpace_Data &localSpace) const
    {
        return localSpace.VanderInternal;
    }

    /// \brief Compute the values of projections of VEM basis functions at points.
    /// \param reference_element_data: an object of type \ref VEM::PCC::VEM_DF_PCC_2D_Pressure_ReferenceElement_Data
    /// which contains monomials stuff. \param polygon: an object of type \ref VEM::PCC::VEM_DF_PCC_2D_Polygon_Geometry
    /// which contains the geoemtric properties of the elements. \param localSpace: an object of type \ref
    /// VEM::PCC::VEM_DF_PCC_2D_Pressure_LocalSpace_Data which contains local matrices. \param projectionType: the \ref
    /// VEM::PCC::ProjectionTypes reporting the kind of projector used to access to the point-wise evalution of VE basis
    /// functions. \param points: a matrix 3 \f$\times\f$ numPoints reporting the coordinates of points. \return A
    /// matrix of size numPoints \f$\times\f$ numDOFs whose columns contain the evaluation of the projection of each
    /// basis function at points.
    inline Eigen::MatrixXd ComputeBasisFunctionsValues(const VEM_DF_PCC_2D_Pressure_ReferenceElement_Data &reference_element_data,
                                                       const VEM_DF_PCC_2D_Pressure_LocalSpace_Data &localSpace,
                                                       const Eigen::MatrixXd &points) const
    {
        return monomials.Vander(reference_element_data.Monomials, points, localSpace.Centroid, localSpace.Diameter);
    }

    /// \brief Compute the values of monomial basis functions at the internal quadrature points.
    /// \param localSpace: an object of type \ref VEM::PCC::VEM_DF_PCC_2D_Pressure_LocalSpace_Data which contains local
    /// matrices. \return A matrix of size numQuadrature \f$\times\f$ numMonomials whose columns contain the evaluation
    /// of monomials at the internal quadrature points.
    inline Eigen::MatrixXd ComputePolynomialsValues(const VEM_DF_PCC_2D_Pressure_LocalSpace_Data &localSpace) const
    {
        return localSpace.VanderInternal;
    }

    /// \brief Compute the values of monomial basis functions at points.
    /// \param reference_element_data: an object of type \ref VEM::PCC::VEM_DF_PCC_2D_Pressure_ReferenceElement_Data
    /// which contains monomials stuff. \param polygon: an object of type \ref VEM::PCC::VEM_DF_PCC_2D_Polygon_Geometry
    /// which contains the geoemtric properties of the elements. \param points: a matrix 3 \f$\times\f$ numPoints
    /// reporting the coordinates of points. \return A matrix of size numPoints \f$\times\f$ numMonomials whose columns
    /// contain the evaluation of monomials at points.
    inline Eigen::MatrixXd ComputePolynomialsValues(const VEM_DF_PCC_2D_Pressure_ReferenceElement_Data &reference_element_data,
                                                    const VEM_DF_PCC_2D_Pressure_LocalSpace_Data &localSpace,
                                                    const Eigen::MatrixXd &points) const
    {
        return monomials.Vander(reference_element_data.Monomials, points, localSpace.Centroid, localSpace.Diameter);
    }
};
} // namespace DF_PCC
} // namespace VEM
} // namespace Polydim

#endif
