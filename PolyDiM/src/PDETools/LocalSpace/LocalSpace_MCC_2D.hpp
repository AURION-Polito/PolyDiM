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

#ifndef __LocalSpace_MCC_2D_H
#define __LocalSpace_MCC_2D_H

#include "I_VEM_MCC_2D_ReferenceElement.hpp"
#include "MeshUtilities.hpp"
#include "QuadratureData.hpp"
#include "VEM_MCC_2D_Creator.hpp"
#include "VEM_MCC_2D_LocalSpace_Data.hpp"
#include "VEM_MCC_PerformanceAnalysis.hpp"

namespace Polydim
{
namespace PDETools
{
namespace LocalSpace_MCC_2D
{
enum class MethodTypes
{
    VEM_MCC = 1,
    VEM_MCC_Partial = 2,
    VEM_MCC_Ortho = 3,
    VEM_MCC_EdgeOrtho = 4,
    VEM_MCC_Ortho_EdgeOrtho = 5
};

class ReferenceElement_Data final
{
public:
    Polydim::PDETools::LocalSpace_MCC_2D::MethodTypes Method_Type;
    unsigned int Order;

    std::unique_ptr<Polydim::VEM::MCC::I_VEM_MCC_2D_Velocity_ReferenceElement> VEM_ReferenceElement_Velocity;
    std::unique_ptr<Polydim::VEM::MCC::I_VEM_MCC_2D_Pressure_ReferenceElement> VEM_ReferenceElement_Pressure;
    Polydim::VEM::MCC::VEM_MCC_2D_Velocity_ReferenceElement_Data VEM_ReferenceElement_Data_Velocity;
    Polydim::VEM::MCC::VEM_MCC_2D_Pressure_ReferenceElement_Data VEM_ReferenceElement_Data_Pressure;
    Polydim::VEM::MCC::VEM_MCC_2D_LocalSpace_Types VEM_Type;
    std::unique_ptr<VEM::MCC::I_VEM_MCC_2D_Velocity_LocalSpace> VEM_LocalSpace_Velocity;
    std::unique_ptr<VEM::MCC::I_VEM_MCC_2D_Pressure_LocalSpace> VEM_LocalSpace_Pressure;
};

class LocalSpace_Data final
{
public:
    Polydim::VEM::MCC::VEM_MCC_2D_Polygon_Geometry VEM_Geometry;
    Polydim::VEM::MCC::VEM_MCC_2D_Velocity_LocalSpace_Data VEM_LocalSpace_Data_Velocity;
    Polydim::VEM::MCC::VEM_MCC_2D_Pressure_LocalSpace_Data VEM_LocalSpace_Data_Pressure;
};

class Performance_Data final
{
public:
    class Cell2D_Performance final
    {
    public:
        unsigned int NumBoundaryQuadraturePoints = 0;
        unsigned int NumInternalQuadraturePoints = 0;
        Polydim::VEM::MCC::VEM_MCC_PerformanceAnalysis_Data Analysis;
    };

    Polydim::PDETools::LocalSpace_MCC_2D::Performance_Data::Cell2D_Performance VEM_Performance_Data;
};

Polydim::PDETools::LocalSpace_MCC_2D::ReferenceElement_Data CreateReferenceElement(const Polydim::PDETools::LocalSpace_MCC_2D::MethodTypes &method_type,
                                                                                   const unsigned int method_order);

std::array<std::array<unsigned int, 4>, 2> ReferenceElementNumDOFs(const Polydim::PDETools::LocalSpace_MCC_2D::ReferenceElement_Data &reference_element_data);

Polydim::PDETools::LocalSpace_MCC_2D::LocalSpace_Data CreateLocalSpace(
    const double &geometric_tolerance_1D,
    const double &geometric_tolerance_2D,
    const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
    const unsigned int cell2D_index,
    const Polydim::PDETools::LocalSpace_MCC_2D::ReferenceElement_Data &reference_element_data);

std::vector<Eigen::MatrixXd> VelocityBasisFunctionsValues(
    const ReferenceElement_Data &reference_element_data,
    const LocalSpace_Data &local_space_data,
    const Polydim::VEM::MCC::ProjectionTypes &projectionType = Polydim::VEM::MCC::ProjectionTypes::Pi0k);

Eigen::MatrixXd PressureBasisFunctionsValues(const Polydim::PDETools::LocalSpace_MCC_2D::ReferenceElement_Data &reference_element_data,
                                             const Polydim::PDETools::LocalSpace_MCC_2D::LocalSpace_Data &local_space_data);

Eigen::MatrixXd VelocityBasisFunctionsDivergenceValues(const Polydim::PDETools::LocalSpace_MCC_2D::ReferenceElement_Data &reference_element_data,
                                                       const Polydim::PDETools::LocalSpace_MCC_2D::LocalSpace_Data &local_space_data);

Eigen::MatrixXd StabilizationMatrix(const Polydim::PDETools::LocalSpace_MCC_2D::ReferenceElement_Data &reference_element_data,
                                    const Polydim::PDETools::LocalSpace_MCC_2D::LocalSpace_Data &local_space_data,
                                    const Polydim::VEM::MCC::ProjectionTypes &projectionType = Polydim::VEM::MCC::ProjectionTypes::Pi0k);

Gedim::Quadrature::QuadratureData EdgeDofsCoordinates(const Polydim::PDETools::LocalSpace_MCC_2D::ReferenceElement_Data &reference_element_data,
                                                      const Polydim::PDETools::LocalSpace_MCC_2D::LocalSpace_Data &local_space_data,
                                                      const unsigned int edge_local_index);

Eigen::VectorXd EdgeDofs(const Polydim::PDETools::LocalSpace_MCC_2D::ReferenceElement_Data &reference_element_data,
                         const Polydim::PDETools::LocalSpace_MCC_2D::LocalSpace_Data &local_space_data,
                         const unsigned int edge_local_index,
                         const Gedim::Quadrature::QuadratureData &edge_dofs_coordinates,
                         const Eigen::VectorXd &strong_values);

Eigen::MatrixXd VelocityBasisFunctionsValuesOnEdges(const unsigned int &edge_local_index,
                                                    const Polydim::PDETools::LocalSpace_MCC_2D::ReferenceElement_Data &reference_element_data,
                                                    const Polydim::PDETools::LocalSpace_MCC_2D::LocalSpace_Data &local_space_data,
                                                    const Eigen::MatrixXd &edge_quadrature_points);

Gedim::Quadrature::QuadratureData EdgeQuadrature(const Polydim::PDETools::LocalSpace_MCC_2D::ReferenceElement_Data &reference_element_data,
                                                 const Polydim::PDETools::LocalSpace_MCC_2D::LocalSpace_Data &local_space_data,
                                                 const unsigned int edge_local_index);

Gedim::Quadrature::QuadratureData InternalQuadrature(const Polydim::PDETools::LocalSpace_MCC_2D::ReferenceElement_Data &reference_element_data,
                                                     const Polydim::PDETools::LocalSpace_MCC_2D::LocalSpace_Data &local_space_data);

unsigned int VelocitySize(const Polydim::PDETools::LocalSpace_MCC_2D::ReferenceElement_Data &reference_element_data,
                          const Polydim::PDETools::LocalSpace_MCC_2D::LocalSpace_Data &local_space_data);

Polydim::PDETools::LocalSpace_MCC_2D::Performance_Data ComputePerformance(
    const Polydim::PDETools::LocalSpace_MCC_2D::ReferenceElement_Data &reference_element_data,
    const Polydim::PDETools::LocalSpace_MCC_2D::LocalSpace_Data &local_space_data);

} // namespace LocalSpace_MCC_2D
} // namespace PDETools
} // namespace Polydim

#endif
