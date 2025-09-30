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

#ifndef __VEM_DF_PCC_3D_Reduced_Pressure_LocalSpace_HPP
#define __VEM_DF_PCC_3D_Reduced_Pressure_LocalSpace_HPP

#include "Eigen/Eigen"
#include "I_VEM_DF_PCC_3D_Pressure_LocalSpace.hpp"
#include "I_VEM_DF_PCC_3D_ReferenceElement.hpp"
#include "Monomials_3D.hpp"
#include "VEM_DF_PCC_3D_LocalSpace_Data.hpp"
#include "VEM_DF_PCC_Utilities.hpp"

namespace Polydim
{
namespace VEM
{
namespace DF_PCC
{

class VEM_DF_PCC_3D_Reduced_Pressure_LocalSpace final : public Polydim::VEM::DF_PCC::I_VEM_DF_PCC_3D_Pressure_LocalSpace
{
  private:
    Polydim::VEM::DF_PCC::VEM_DF_PCC_Utilities<3> utilities;
    Polydim::Utilities::Monomials_3D monomials;

    void InitializeProjectorsComputation(const VEM_DF_PCC_3D_Pressure_ReferenceElement_Data &reference_element_data,
                                         const Eigen::Vector3d &polyhedronCentroid,
                                         const double &polyhedronDiameter,
                                         const Eigen::MatrixXd &internalQuadraturePoints,
                                         const Eigen::VectorXd &internalQuadratureWeights,
                                         Polydim::VEM::DF_PCC::VEM_DF_PCC_3D_Pressure_LocalSpace_Data &localSpace) const;

  public:
    Polydim::VEM::DF_PCC::VEM_DF_PCC_3D_Pressure_LocalSpace_Data CreateLocalSpace(
        const Polydim::VEM::DF_PCC::VEM_DF_PCC_3D_Pressure_ReferenceElement_Data &reference_element_data_3D,
        const Polydim::VEM::DF_PCC::VEM_DF_PCC_3D_Polyhedron_Geometry &polyhedron) const;

    inline Eigen::MatrixXd ComputeBasisFunctionsValues(const Polydim::VEM::DF_PCC::VEM_DF_PCC_3D_Pressure_LocalSpace_Data &localSpace) const
    {
        return localSpace.VanderInternal;
    }

    inline Eigen::MatrixXd ComputeBasisFunctionsValues(const Polydim::VEM::DF_PCC::VEM_DF_PCC_3D_Pressure_ReferenceElement_Data &reference_element_data,
                                                       const Polydim::VEM::DF_PCC::VEM_DF_PCC_3D_Pressure_LocalSpace_Data &localSpace,
                                                       const Eigen::MatrixXd &points) const
    {
        return monomials.Vander(reference_element_data.Monomials, points, localSpace.Centroid, localSpace.Diameter);
    }

    inline Eigen::MatrixXd ComputePolynomialsValues(const Polydim::VEM::DF_PCC::VEM_DF_PCC_3D_Pressure_LocalSpace_Data &localSpace) const
    {
        return localSpace.VanderInternal;
    }

    inline Eigen::MatrixXd ComputePolynomialsValues(const Polydim::VEM::DF_PCC::VEM_DF_PCC_3D_Pressure_ReferenceElement_Data &reference_element_data,
                                                    const Polydim::VEM::DF_PCC::VEM_DF_PCC_3D_Pressure_LocalSpace_Data &localSpace,
                                                    const Eigen::MatrixXd &points) const
    {
        return monomials.Vander(reference_element_data.Monomials, points, localSpace.Centroid, localSpace.Diameter);
    }
};
} // namespace DF_PCC
} // namespace VEM
} // namespace Polydim

#endif
