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

#ifndef __LocalSpace_DF_PCC_2D_H
#define __LocalSpace_DF_PCC_2D_H

#include "Assembler_Utilities.hpp"
#include "DOFsManager.hpp"
#include "FEM_PCC_2D_LocalSpace.hpp"
#include "IArray.hpp"
#include "I_VEM_DF_PCC_2D_ReferenceElement.hpp"
#include "VEM_DF_PCC_2D_Creator.hpp"
#include "VEM_DF_PCC_2D_LocalSpace_Data.hpp"
#include "VEM_DF_PCC_PerformanceAnalysis.hpp"
#include <memory>

namespace Polydim
{
namespace PDETools
{
namespace LocalSpace_DF_PCC_2D
{
enum struct MethodTypes
{
    TAYLOR_HOOD = 0,
    VEM_DF_PCC_FULL = 1,
    VEM_DF_PCC_REDUCED = 2
};

class ReferenceElement_Data final
{
  public:
    Polydim::PDETools::LocalSpace_DF_PCC_2D::MethodTypes Method_Type;
    unsigned int Order;
    unsigned int Dimension = 2;

    Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_LocalSpace_Types VEM_Type;
    std::unique_ptr<Polydim::VEM::DF_PCC::I_VEM_DF_PCC_2D_Velocity_ReferenceElement> VEM_Velocity_ReferenceElement;
    Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_ReferenceElement_Data VEM_Velocity_ReferenceElement_Data;
    std::unique_ptr<Polydim::VEM::DF_PCC::I_VEM_DF_PCC_2D_Pressure_ReferenceElement> VEM_Pressure_ReferenceElement;
    Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Pressure_ReferenceElement_Data VEM_Pressure_ReferenceElement_Data;
    std::unique_ptr<Polydim::VEM::DF_PCC::I_VEM_DF_PCC_2D_Velocity_LocalSpace> VEM_Velocity_LocalSpace;
    std::unique_ptr<Polydim::VEM::DF_PCC::I_VEM_DF_PCC_2D_Pressure_LocalSpace> VEM_Pressure_LocalSpace;

    std::unique_ptr<Polydim::FEM::PCC::FEM_PCC_2D_ReferenceElement> FEM_ReferenceElement;
    Polydim::FEM::PCC::FEM_PCC_2D_ReferenceElement_Data FEM_ReferenceElement_Data;
    std::unique_ptr<Polydim::FEM::PCC::FEM_PCC_2D_LocalSpace> FEM_LocalSpace;
};

class LocalSpace_Data final
{
  public:
    Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Polygon_Geometry VEM_Geometry;
    Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_LocalSpace_Data VEM_Velocity_LocalSpace_Data;
    Polydim::VEM::DF_PCC::VEM_DF_PCC_2D_Pressure_LocalSpace_Data VEM_Pressure_LocalSpace_Data;

    Polydim::FEM::PCC::FEM_PCC_2D_Polygon_Geometry FEM_Geometry;
    Polydim::FEM::PCC::FEM_PCC_2D_LocalSpace_Data FEM_LocalSpace_Data;
};

class Performance_Data final
{
  public:
    class Cell2D_Performance final
    {
      public:
        unsigned int NumBoundaryQuadraturePoints = 0;
        unsigned int NumInternalQuadraturePoints = 0;
        Polydim::VEM::DF_PCC::VEM_DF_PCC_PerformanceAnalysis_Data vem_analysis_data;
    };

    Polydim::PDETools::LocalSpace_DF_PCC_2D::Performance_Data::Cell2D_Performance performance_data;
};

Polydim::PDETools::LocalSpace_DF_PCC_2D::ReferenceElement_Data CreateReferenceElement(const Polydim::PDETools::LocalSpace_DF_PCC_2D::MethodTypes &method_type,
                                                                                      const unsigned int method_order);

std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> SetMeshDOFsInfo(
    const Polydim::PDETools::LocalSpace_DF_PCC_2D::ReferenceElement_Data &reference_element_data,
    const Gedim::MeshMatricesDAO &mesh,
    const std::array<std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo>, 2> &boundary_info);

LocalSpace_Data CreateLocalSpace(const double &geometric_tolerance_1D,
                                 const double &geometric_tolerance_2D,
                                 const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                 const unsigned int cell2D_index,
                                 const Polydim::PDETools::LocalSpace_DF_PCC_2D::ReferenceElement_Data &reference_element_data);

std::vector<Eigen::MatrixXd> VelocityBasisFunctionsValues(
    const Polydim::PDETools::LocalSpace_DF_PCC_2D::ReferenceElement_Data &reference_element_data,
    const Polydim::PDETools::LocalSpace_DF_PCC_2D::LocalSpace_Data &local_space_data,
    const Polydim::VEM::DF_PCC::ProjectionTypes &projectionType = Polydim::VEM::DF_PCC::ProjectionTypes::Pi0k);

std::vector<Eigen::MatrixXd> VelocityBasisFunctionsValues(
    const ReferenceElement_Data &reference_element_data,
    const LocalSpace_Data &local_space_data,
    const Eigen::MatrixXd &points,
    const Polydim::VEM::DF_PCC::ProjectionTypes &projectionType = Polydim::VEM::DF_PCC::ProjectionTypes::Pi0k);

Eigen::MatrixXd PressureBasisFunctionsValues(const Polydim::PDETools::LocalSpace_DF_PCC_2D::ReferenceElement_Data &reference_element_data,
                                             const Polydim::PDETools::LocalSpace_DF_PCC_2D::LocalSpace_Data &local_space_data);

Eigen::MatrixXd PressureBasisFunctionsValues(const ReferenceElement_Data &reference_element_data,
                                             const LocalSpace_Data &local_space_data,
                                             const Eigen::MatrixXd &points);

std::vector<Eigen::MatrixXd> VelocityBasisFunctionsDerivativeValues(
    const Polydim::PDETools::LocalSpace_DF_PCC_2D::ReferenceElement_Data &reference_element_data,
    const Polydim::PDETools::LocalSpace_DF_PCC_2D::LocalSpace_Data &local_space_data,
    const Polydim::VEM::DF_PCC::ProjectionTypes &projectionType = Polydim::VEM::DF_PCC::ProjectionTypes::PiNabla);

std::vector<Eigen::MatrixXd> VelocityBasisFunctionsDerivativeValues(
    const Polydim::PDETools::LocalSpace_DF_PCC_2D::ReferenceElement_Data &reference_element_data,
    const Polydim::PDETools::LocalSpace_DF_PCC_2D::LocalSpace_Data &local_space_data,
    const Eigen::MatrixXd &points,
    const Polydim::VEM::DF_PCC::ProjectionTypes &projectionType = Polydim::VEM::DF_PCC::ProjectionTypes::PiNabla);

Eigen::MatrixXd VelocityBasisFunctionsDivergenceValues(const Polydim::PDETools::LocalSpace_DF_PCC_2D::ReferenceElement_Data &reference_element_data,
                                                       const Polydim::PDETools::LocalSpace_DF_PCC_2D::LocalSpace_Data &local_space_data);

Gedim::Quadrature::QuadratureData InternalQuadrature(const Polydim::PDETools::LocalSpace_DF_PCC_2D::ReferenceElement_Data &reference_element_data,
                                                     const Polydim::PDETools::LocalSpace_DF_PCC_2D::LocalSpace_Data &local_space_data);

Eigen::MatrixXd VelocityBasisFunctionsValuesOnEdge(const unsigned int &edge_local_index,
                                                   const Polydim::PDETools::LocalSpace_DF_PCC_2D::ReferenceElement_Data &reference_element_data,
                                                   const Polydim::PDETools::LocalSpace_DF_PCC_2D::LocalSpace_Data &local_space_data,
                                                   const Eigen::MatrixXd &pointsCurvilinearCoordinates);

Eigen::MatrixXd VelocityStabilizationMatrix(
    const Polydim::PDETools::LocalSpace_DF_PCC_2D::ReferenceElement_Data &reference_element_data,
    const Polydim::PDETools::LocalSpace_DF_PCC_2D::LocalSpace_Data &local_space_data,
    const Polydim::VEM::DF_PCC::ProjectionTypes &projectionType = Polydim::VEM::DF_PCC::ProjectionTypes::PiNabla);

Eigen::MatrixXd VelocityEdgeDofsCoordinates(const Polydim::PDETools::LocalSpace_DF_PCC_2D::ReferenceElement_Data &reference_element_data,
                                            const Polydim::PDETools::LocalSpace_DF_PCC_2D::LocalSpace_Data &local_space_data,
                                            const unsigned int edge_local_index);

unsigned int VelocitySize(const Polydim::PDETools::LocalSpace_DF_PCC_2D::ReferenceElement_Data &reference_element_data,
                          const Polydim::PDETools::LocalSpace_DF_PCC_2D::LocalSpace_Data &local_space_data);

Polydim::PDETools::LocalSpace_DF_PCC_2D::Performance_Data ComputePerformance(
    const Polydim::PDETools::LocalSpace_DF_PCC_2D::ReferenceElement_Data &reference_element_data,
    const Polydim::PDETools::LocalSpace_DF_PCC_2D::LocalSpace_Data &local_space_data);

void export_velocity_dofs(const Gedim::GeometryUtilities &geometry_utilities,
                          const Gedim::MeshMatricesDAO &mesh,
                          const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                          const std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> &mesh_dofs_info,
                          const std::vector<DOFs::DOFsManager::DOFsData> &dofs_data,
                          const Polydim::PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                          const Gedim::IArray &right_hand_side,
                          const Gedim::IArray &solution,
                          const Gedim::IArray &solution_strongs,
                          const std::string &file_path);

} // namespace LocalSpace_DF_PCC_2D
} // namespace PDETools
} // namespace Polydim

#endif
