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

#ifndef __LocalSpace_PCC_2D_H
#define __LocalSpace_PCC_2D_H

#include "DOFsManager.hpp"
#include "FEM_PCC_2D_LocalSpace.hpp"
#include "I_VEM_PCC_2D_ReferenceElement.hpp"
#include "MeshMatricesDAO.hpp"
#include "MeshUtilities.hpp"
#include "QuadratureData.hpp"
#include "VEM_PCC_2D_Creator.hpp"
#include "VEM_PCC_2D_LocalSpace_Data.hpp"
#include "VEM_PCC_PerformanceAnalysis.hpp"

namespace Polydim
{
namespace PDETools
{
namespace LocalSpace_PCC_2D
{
enum struct MethodTypes
{
    FEM_PCC = 0,
    VEM_PCC = 1,
    VEM_PCC_Inertia = 2,
    VEM_PCC_Ortho = 3
};

class ReferenceElement_Data final
{
  public:
    Polydim::PDETools::LocalSpace_PCC_2D::MethodTypes Method_Type;
    unsigned int Order;

    std::unique_ptr<Polydim::VEM::PCC::I_VEM_PCC_2D_ReferenceElement> VEM_ReferenceElement;
    Polydim::VEM::PCC::VEM_PCC_2D_ReferenceElement_Data VEM_ReferenceElement_Data;
    Polydim::VEM::PCC::VEM_PCC_2D_LocalSpace_Types VEM_Type;
    std::unique_ptr<Polydim::VEM::PCC::I_VEM_PCC_2D_LocalSpace> VEM_LocalSpace;

    std::unique_ptr<Polydim::FEM::PCC::FEM_PCC_2D_ReferenceElement> FEM_ReferenceElement;
    Polydim::FEM::PCC::FEM_PCC_2D_ReferenceElement_Data FEM_ReferenceElement_Data;
    std::unique_ptr<Polydim::FEM::PCC::FEM_PCC_2D_LocalSpace> FEM_LocalSpace;
};

class LocalSpace_Data final
{
  public:
    Polydim::VEM::PCC::VEM_PCC_2D_Polygon_Geometry VEM_Geometry;
    Polydim::VEM::PCC::VEM_PCC_2D_LocalSpace_Data VEM_LocalSpace_Data;

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
        Polydim::VEM::PCC::VEM_PCC_PerformanceAnalysis_Data Analysis;
    };

    Polydim::PDETools::LocalSpace_PCC_2D::Performance_Data::Cell2D_Performance VEM_Performance_Data;
};

Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data CreateReferenceElement(const Polydim::PDETools::LocalSpace_PCC_2D::MethodTypes &method_type,
                                                                                   const unsigned int method_order);

Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo SetMeshDOFsInfo(
    const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
    const Gedim::MeshMatricesDAO &mesh,
    const std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> &boundary_info);

LocalSpace_Data CreateLocalSpace(const double &geometric_tolerance_1D,
                                 const double &geometric_tolerance_2D,
                                 const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                 const unsigned int cell2D_index,
                                 const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data);

Eigen::MatrixXd BasisFunctionsValues(const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                                     const Polydim::PDETools::LocalSpace_PCC_2D::LocalSpace_Data &local_space_data,
                                     const Polydim::VEM::PCC::ProjectionTypes &projectionType = Polydim::VEM::PCC::ProjectionTypes::Pi0km1);

Eigen::MatrixXd BasisFunctionsValuesOnEdge(const unsigned int &edge_local_index,
                                           const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                                           const Polydim::PDETools::LocalSpace_PCC_2D::LocalSpace_Data &local_space_data,
                                           const Eigen::MatrixXd &pointsCurvilinearCoordinates);

std::vector<Eigen::MatrixXd> BasisFunctionsDerivativeValues(
    const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
    const Polydim::PDETools::LocalSpace_PCC_2D::LocalSpace_Data &local_space_data,
    const Polydim::VEM::PCC::ProjectionTypes &projectionType = Polydim::VEM::PCC::ProjectionTypes::Pi0km1Der);

Eigen::MatrixXd BasisFunctionsLaplacianValues(
    const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
    const Polydim::PDETools::LocalSpace_PCC_2D::LocalSpace_Data &local_space_data,
    const Polydim::VEM::PCC::ProjectionTypes &projectionType = Polydim::VEM::PCC::ProjectionTypes::Pi0km1Der);

Eigen::MatrixXd StabilizationMatrix(const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                                    const Polydim::PDETools::LocalSpace_PCC_2D::LocalSpace_Data &local_space_data,
                                    const Polydim::VEM::PCC::ProjectionTypes &projectionType = Polydim::VEM::PCC::ProjectionTypes::PiNabla);

Eigen::MatrixXd EdgeDofsCoordinates(const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                                    const Polydim::PDETools::LocalSpace_PCC_2D::LocalSpace_Data &local_space_data,
                                    const unsigned int edge_local_index);

Gedim::Quadrature::QuadratureData InternalDofsCoordinates(const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                                                          const Polydim::PDETools::LocalSpace_PCC_2D::LocalSpace_Data &local_space_data);

Eigen::VectorXd InternalDofs(const ReferenceElement_Data &reference_element_data,
                             const Polydim::PDETools::LocalSpace_PCC_2D::LocalSpace_Data &local_space_data,
                             const Eigen::VectorXd &values_at_dofs,
                             const Gedim::Quadrature::QuadratureData &internal_dofs_coordinates);

Gedim::Quadrature::QuadratureData InternalQuadrature(const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                                                     const Polydim::PDETools::LocalSpace_PCC_2D::LocalSpace_Data &local_space_data);

unsigned int Size(const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                  const Polydim::PDETools::LocalSpace_PCC_2D::LocalSpace_Data &local_space_data);

Polydim::PDETools::LocalSpace_PCC_2D::Performance_Data ComputePerformance(
    const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
    const Polydim::PDETools::LocalSpace_PCC_2D::LocalSpace_Data &local_space_data);
} // namespace LocalSpace_PCC_2D
} // namespace PDETools
} // namespace Polydim

#endif
