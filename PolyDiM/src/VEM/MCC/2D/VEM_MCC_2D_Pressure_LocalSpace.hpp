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

#ifndef __VEM_MCC_2D_Presure_LocalSpace_HPP
#define __VEM_MCC_2D_Presure_LocalSpace_HPP

#include "Eigen/Eigen"
#include "I_VEM_MCC_2D_Pressure_LocalSpace.hpp"
#include "I_VEM_MCC_2D_ReferenceElement.hpp"
#include "Monomials_2D.hpp"
#include "VEM_MCC_2D_LocalSpace_Data.hpp"

namespace Polydim
{
namespace VEM
{
namespace MCC
{

class VEM_MCC_2D_Pressure_LocalSpace final : public I_VEM_MCC_2D_Pressure_LocalSpace
{
  private:
    Polydim::Utilities::Monomials_2D monomials;

    void InitializeProjectorsComputation(const Polydim::VEM::MCC::VEM_MCC_2D_Pressure_ReferenceElement_Data &reference_element_data,
                                         const Eigen::Vector3d &polygonCentroid,
                                         const double &polygonDiameter,
                                         const Eigen::MatrixXd &internalQuadraturePoints,
                                         const Eigen::VectorXd &internalQuadratureWeights,
                                         Polydim::VEM::MCC::VEM_MCC_2D_Pressure_LocalSpace_Data &localSpace) const;

  public:
    virtual ~VEM_MCC_2D_Pressure_LocalSpace()
    {
    }

    Polydim::VEM::MCC::VEM_MCC_2D_Pressure_LocalSpace_Data CreateLocalSpace(
        const Polydim::VEM::MCC::VEM_MCC_2D_Pressure_ReferenceElement_Data &reference_element_data,
        const Polydim::VEM::MCC::VEM_MCC_2D_Polygon_Geometry &polygon) const;

    inline Eigen::MatrixXd ComputeBasisFunctionsValues(const Polydim::VEM::MCC::VEM_MCC_2D_Pressure_LocalSpace_Data &localSpace) const
    {
        return localSpace.VanderInternal * localSpace.Qmatrix.transpose();
    }

    inline Eigen::MatrixXd ComputeBasisFunctionsValues(const Polydim::VEM::MCC::VEM_MCC_2D_Pressure_ReferenceElement_Data &reference_element_data,
                                                       const Polydim::VEM::MCC::VEM_MCC_2D_Pressure_LocalSpace_Data &localSpace,
                                                       const Eigen::MatrixXd &points) const
    {
        return monomials.Vander(reference_element_data.Monomials, points, localSpace.Centroid, localSpace.Diameter) *
               localSpace.Qmatrix.transpose();
    }

    inline Eigen::MatrixXd ComputePolynomialsValues(const Polydim::VEM::MCC::VEM_MCC_2D_Pressure_LocalSpace_Data &localSpace) const
    {
        return localSpace.VanderInternal * localSpace.Qmatrix.transpose();
    }

    inline Eigen::MatrixXd ComputePolynomialsValues(const Polydim::VEM::MCC::VEM_MCC_2D_Pressure_ReferenceElement_Data &reference_element_data,
                                                    const Polydim::VEM::MCC::VEM_MCC_2D_Pressure_LocalSpace_Data &localSpace,
                                                    const Eigen::MatrixXd &points) const
    {
        return monomials.Vander(reference_element_data.Monomials, points, localSpace.Centroid, localSpace.Diameter) *
               localSpace.Qmatrix.transpose();
    }
};
} // namespace MCC
} // namespace VEM
} // namespace Polydim

#endif
