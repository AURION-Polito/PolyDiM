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

#ifndef __FEM_Triangle_RT_MCC_2D_LocalSpace_HPP
#define __FEM_Triangle_RT_MCC_2D_LocalSpace_HPP

#include "FEM_MCC_2D_LocalSpace_Data.hpp"
#include "FEM_Triangle_RT_MCC_2D_ReferenceElement.hpp"
#include "MapTriangle.hpp"

namespace Polydim
{
namespace FEM
{
namespace MCC
{

class FEM_Triangle_RT_MCC_2D_LocalSpace final
{
  private:
    std::vector<Eigen::MatrixXd> MapVelocityValues(const Polydim::FEM::MCC::FEM_Triangle_RT_MCC_2D_LocalSpace_Data &local_space,
                                                   const std::vector<Eigen::MatrixXd> &referenceValues) const;

    Eigen::MatrixXd MapPressureValues(const Polydim::FEM::MCC::FEM_Triangle_RT_MCC_2D_LocalSpace_Data &local_space,
                                      const Eigen::MatrixXd &referenceValues) const;

    Eigen::MatrixXd MapVelocityDivergenceValues(const Polydim::FEM::MCC::FEM_Triangle_RT_MCC_2D_LocalSpace_Data &local_space,
                                                const Eigen::MatrixXd &referenceDerivateValues) const;

    Gedim::Quadrature::QuadratureData InternalQuadrature(const Gedim::Quadrature::QuadratureData &reference_quadrature,
                                                         const Gedim::MapTriangle::MapTriangleData &mapData) const;

    std::vector<Gedim::Quadrature::QuadratureData> BoundaryQuadrature(const Gedim::Quadrature::QuadratureData &reference_quadrature,
                                                                      const Polydim::FEM::MCC::FEM_MCC_2D_Polygon_Geometry &polygon) const;

  public:
    Polydim::FEM::MCC::FEM_Triangle_RT_MCC_2D_LocalSpace_Data CreateLocalSpace(
        const Polydim::FEM::MCC::FEM_Triangle_RT_MCC_2D_ReferenceElement_Data &reference_element_data,
        const Polydim::FEM::MCC::FEM_MCC_2D_Polygon_Geometry &polygon) const;

    Eigen::MatrixXd ComputePressureBasisFunctionsValues(const Polydim::FEM::MCC::FEM_Triangle_RT_MCC_2D_ReferenceElement_Data &reference_element_data,
                                                        const Polydim::FEM::MCC::FEM_Triangle_RT_MCC_2D_LocalSpace_Data &local_space) const
    {
        return MapPressureValues(local_space, reference_element_data.reference_element_data_pressure.ReferenceBasisFunctionValues);
    }

    std::vector<Eigen::MatrixXd> ComputeVelocityBasisFunctionsValues(
        const Polydim::FEM::MCC::FEM_Triangle_RT_MCC_2D_ReferenceElement_Data &reference_element_data,
        const Polydim::FEM::MCC::FEM_Triangle_RT_MCC_2D_LocalSpace_Data &local_space) const
    {
        return MapVelocityValues(local_space, reference_element_data.reference_element_data_velocity.ReferenceBasisFunctionValues);
    }

    Eigen::MatrixXd ComputeVelocityBasisFunctionsDivergenceValues(
        const Polydim::FEM::MCC::FEM_Triangle_RT_MCC_2D_ReferenceElement_Data &reference_element_data,
        const Polydim::FEM::MCC::FEM_Triangle_RT_MCC_2D_LocalSpace_Data &local_space) const
    {
        return MapVelocityDivergenceValues(local_space, reference_element_data.reference_element_data_velocity.ReferenceBasisFunctionDivergenceValues);
    }

    Eigen::MatrixXd ComputePressureBasisFunctionsValues(const Polydim::FEM::MCC::FEM_Triangle_RT_MCC_2D_ReferenceElement_Data &reference_element_data,
                                                        const Polydim::FEM::MCC::FEM_Triangle_RT_MCC_2D_LocalSpace_Data &local_space,
                                                        const Eigen::MatrixXd &points) const
    {
        Gedim::MapTriangle mapTriangle;
        const Eigen::MatrixXd referencePoints = mapTriangle.FInv(local_space.MapData, points);

        Polydim::FEM::MCC::FEM_Triangle_RT_MCC_2D_ReferenceElement reference_element;

        return MapPressureValues(local_space, reference_element.EvaluatePressureBasisFunctions(referencePoints, reference_element_data));
    }

    std::vector<Eigen::MatrixXd> ComputeVelocityBasisFunctionsValues(
        const Polydim::FEM::MCC::FEM_Triangle_RT_MCC_2D_ReferenceElement_Data &reference_element_data,
        const Polydim::FEM::MCC::FEM_Triangle_RT_MCC_2D_LocalSpace_Data &local_space,
        const Eigen::MatrixXd &points) const
    {
        Gedim::MapTriangle mapTriangle;
        const Eigen::MatrixXd referencePoints = mapTriangle.FInv(local_space.MapData, points);

        Polydim::FEM::MCC::FEM_Triangle_RT_MCC_2D_ReferenceElement reference_element;

        return MapVelocityValues(local_space, reference_element.EvaluateVelocityBasisFunctions(referencePoints, reference_element_data));
    }

    Eigen::MatrixXd ComputeVelocityBasisFunctionsDivergenceValues(
        const Polydim::FEM::MCC::FEM_Triangle_RT_MCC_2D_ReferenceElement_Data &reference_element_data,
        const Polydim::FEM::MCC::FEM_Triangle_RT_MCC_2D_LocalSpace_Data &local_space,
        const Eigen::MatrixXd &points) const
    {
        Gedim::MapTriangle mapTriangle;
        const Eigen::MatrixXd referencePoints = mapTriangle.FInv(local_space.MapData, points);

        Polydim::FEM::MCC::FEM_Triangle_RT_MCC_2D_ReferenceElement reference_element;

        return MapVelocityDivergenceValues(
            local_space,
            reference_element.EvaluateVelociytBasisFunctionsDivergence(referencePoints, reference_element_data));
    }
};
} // namespace MCC
} // namespace FEM
} // namespace Polydim

#endif
