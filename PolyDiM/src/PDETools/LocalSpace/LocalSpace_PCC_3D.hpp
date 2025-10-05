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

#ifndef __LocalSpace_PCC_3D_H
#define __LocalSpace_PCC_3D_H

#include "DOFsManager.hpp"
#include "FEM_PCC_3D_LocalSpace.hpp"
#include "I_VEM_PCC_3D_ReferenceElement.hpp"
#include "MeshMatricesDAO.hpp"
#include "MeshUtilities.hpp"
#include "QuadratureData.hpp"
#include "VEM_PCC_3D_Creator.hpp"
#include "VEM_PCC_3D_LocalSpace_Data.hpp"
#include "VEM_PCC_PerformanceAnalysis.hpp"

namespace Polydim
{
namespace PDETools
{
namespace LocalSpace_PCC_3D
{
enum struct MethodTypes
{
    FEM_PCC = 0,
    VEM_PCC = 1,
    VEM_PCC_Inertia = 2,
    VEM_PCC_Ortho = 3
};

struct ReferenceElement_Data final
{
    Polydim::PDETools::LocalSpace_PCC_3D::MethodTypes Method_Type;
    unsigned int Order;

    std::unique_ptr<Polydim::VEM::PCC::I_VEM_PCC_2D_ReferenceElement> VEM_ReferenceElement_2D;
    Polydim::VEM::PCC::VEM_PCC_2D_ReferenceElement_Data VEM_ReferenceElement_Data_2D;
    std::unique_ptr<Polydim::VEM::PCC::I_VEM_PCC_3D_ReferenceElement> VEM_ReferenceElement_3D;
    Polydim::VEM::PCC::VEM_PCC_3D_ReferenceElement_Data VEM_ReferenceElement_Data_3D;
    VEM::PCC::VEM_PCC_3D_LocalSpace_Types VEM_Type;
    std::unique_ptr<VEM::PCC::I_VEM_PCC_3D_LocalSpace> VEM_LocalSpace;

    std::unique_ptr<Polydim::FEM::PCC::FEM_PCC_3D_ReferenceElement> FEM_ReferenceElement_3D;
    Polydim::FEM::PCC::FEM_PCC_3D_ReferenceElement_Data FEM_ReferenceElement_Data_3D;
    std::unique_ptr<FEM::PCC::FEM_PCC_3D_LocalSpace> FEM_LocalSpace;
};

struct LocalSpace_Data final
{
    struct VEM_Geometry final
    {
        std::vector<Polydim::VEM::PCC::VEM_PCC_2D_Polygon_Geometry> PolygonalFaces;
        Polydim::VEM::PCC::VEM_PCC_3D_Polyhedron_Geometry Polyhedron;
    };

    Polydim::PDETools::LocalSpace_PCC_3D::LocalSpace_Data::VEM_Geometry VEM_Geometry;
    Polydim::VEM::PCC::VEM_PCC_3D_LocalSpace_Data VEM_LocalSpace_Data;

    Polydim::FEM::PCC::FEM_PCC_3D_Polyhedron_Geometry FEM_Geometry;
    Polydim::FEM::PCC::FEM_PCC_3D_LocalSpace_Data FEM_LocalSpace_Data;
};

struct Performance_Data final
{
    struct Cell3D_Performance final
    {
        unsigned int NumBoundaryQuadraturePoints = 0;
        unsigned int NumInternalQuadraturePoints = 0;
        Polydim::VEM::PCC::VEM_PCC_PerformanceAnalysis_Data Analysis;
    };

    Polydim::PDETools::LocalSpace_PCC_3D::Performance_Data::Cell3D_Performance VEM_Performance_Data;
};

Polydim::PDETools::LocalSpace_PCC_3D::ReferenceElement_Data CreateReferenceElement(const Polydim::PDETools::LocalSpace_PCC_3D::MethodTypes &method_type,
                                                                                   const unsigned int method_order);

PDETools::DOFs::DOFsManager::MeshDOFsInfo SetMeshDOFsInfo(
    const Polydim::PDETools::LocalSpace_PCC_3D::ReferenceElement_Data &reference_element_data,
    const Gedim::MeshMatricesDAO &mesh,
    const std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> &boundary_info);

Polydim::PDETools::LocalSpace_PCC_3D::LocalSpace_Data CreateLocalSpace(
    const double &geometric_tolerance_1D,
    const double &geometric_tolerance_2D,
    const double &geometric_tolerance_3D,
    const Gedim::MeshUtilities::MeshGeometricData3D &mesh_geometric_data,
    const unsigned int cell3D_index,
    const Polydim::PDETools::LocalSpace_PCC_3D::ReferenceElement_Data &reference_element_data);

Eigen::MatrixXd BasisFunctionsValues(const Polydim::PDETools::LocalSpace_PCC_3D::ReferenceElement_Data &reference_element_data,
                                     const Polydim::PDETools::LocalSpace_PCC_3D::LocalSpace_Data &local_space_data,
                                     const Polydim::VEM::PCC::ProjectionTypes &projectionType = Polydim::VEM::PCC::ProjectionTypes::Pi0km1);

Eigen::MatrixXd BasisFunctionsValuesOnFace(const unsigned int &face_local_index,
                                           const Polydim::PDETools::LocalSpace_PCC_3D::ReferenceElement_Data &reference_element_data,
                                           const Polydim::PDETools::LocalSpace_PCC_3D::LocalSpace_Data &local_space_data,
                                           const Eigen::MatrixXd &quadrature_points);

std::vector<Eigen::MatrixXd> BasisFunctionsDerivativeValues(
    const Polydim::PDETools::LocalSpace_PCC_3D::ReferenceElement_Data &reference_element_data,
    const Polydim::PDETools::LocalSpace_PCC_3D::LocalSpace_Data &local_space_data,
    const Polydim::VEM::PCC::ProjectionTypes &projectionType = Polydim::VEM::PCC::ProjectionTypes::Pi0km1Der);

Eigen::MatrixXd StabilizationMatrix(const Polydim::PDETools::LocalSpace_PCC_3D::ReferenceElement_Data &reference_element_data,
                                    const Polydim::PDETools::LocalSpace_PCC_3D::LocalSpace_Data &local_space_data,
                                    const Polydim::VEM::PCC::ProjectionTypes &projectionType = Polydim::VEM::PCC::ProjectionTypes::PiNabla);

Eigen::MatrixXd EdgeDofsCoordinates(const Polydim::PDETools::LocalSpace_PCC_3D::ReferenceElement_Data &reference_element_data,
                                    const Polydim::PDETools::LocalSpace_PCC_3D::LocalSpace_Data &local_space_data,
                                    const unsigned int edge_local_index);

Gedim::Quadrature::QuadratureData FaceDofsCoordinates(const Polydim::PDETools::LocalSpace_PCC_3D::ReferenceElement_Data &reference_element_data,
                                                      const Polydim::PDETools::LocalSpace_PCC_3D::LocalSpace_Data &local_space_data,
                                                      const unsigned int face_local_index,
                                                      unsigned int &quadrature_offset);

Eigen::VectorXd FaceDofs(const Polydim::PDETools::LocalSpace_PCC_3D::ReferenceElement_Data &reference_element_data,
                         const Polydim::PDETools::LocalSpace_PCC_3D::LocalSpace_Data &local_space_data,
                         const unsigned int face_local_index,
                         const Eigen::VectorXd &strong_values,
                         const Gedim::Quadrature::QuadratureData &quadrature_offset);

Gedim::Quadrature::QuadratureData FaceQuadrature(const Polydim::PDETools::LocalSpace_PCC_3D::ReferenceElement_Data &reference_element_data,
                                                 const Polydim::PDETools::LocalSpace_PCC_3D::LocalSpace_Data &local_space_data,
                                                 const unsigned int face_local_index,
                                                 unsigned int &quadrature_offset);

Gedim::Quadrature::QuadratureData InternalQuadrature(const Polydim::PDETools::LocalSpace_PCC_3D::ReferenceElement_Data &reference_element_data,
                                                     const Polydim::PDETools::LocalSpace_PCC_3D::LocalSpace_Data &local_space_data);

unsigned int Size(const Polydim::PDETools::LocalSpace_PCC_3D::ReferenceElement_Data &reference_element_data,
                  const Polydim::PDETools::LocalSpace_PCC_3D::LocalSpace_Data &local_space_data);

Polydim::PDETools::LocalSpace_PCC_3D::Performance_Data ComputePerformance(
    const Polydim::PDETools::LocalSpace_PCC_3D::ReferenceElement_Data &reference_element_data,
    const Polydim::PDETools::LocalSpace_PCC_3D::LocalSpace_Data &local_space_data);
} // namespace LocalSpace_PCC_3D
} // namespace PDETools
} // namespace Polydim

#endif
