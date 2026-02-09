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

#ifndef __TEST_assembler_PCC_2D_H
#define __TEST_assembler_PCC_2D_H

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "MeshMatricesDAO_mesh_connectivity_data.hpp"
#include "PDE_Mesh_Utilities.hpp"
#include "LocalSpace_PCC_2D.hpp"
#include "MeshUtilities.hpp"
#include "assembler_PCC_2D_functions.hpp"

namespace Polydim
{
  namespace UnitTesting
  {

    TEST(TEST_assembler_PCC_2D, TEST_assembler_PCC_2D_forcing_term)
    {
      Gedim::GeometryUtilitiesConfig geometry_utilities_config;
      geometry_utilities_config.Tolerance1D = 1.0e-8;
      geometry_utilities_config.Tolerance2D = 1.0e-12;
      Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

      Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D domain;
      domain.area = 1.0;
      domain.vertices = Eigen::MatrixXd::Zero(3, 4);
      domain.vertices.row(0) << 0.0, 1.0, 1.0, 0.0;
      domain.vertices.row(1) << 0.0, 0.0, 1.0, 1.0;
      domain.shape_type = Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D::Domain_Shape_Types::Parallelogram;


      Gedim::MeshUtilities mesh_utilities;

      Gedim::MeshMatrices meshData;
      Gedim::MeshMatricesDAO mesh(meshData);

      Polydim::PDETools::Mesh::PDE_Mesh_Utilities::create_mesh_2D(geometry_utilities,
                                                                  mesh_utilities,
                                                                  Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::Triangular,
                                                                  domain,
                                                                  0.1,
                                                                  mesh);

      const unsigned int method_order = 1;
      const auto reference_element_data =
          Polydim::PDETools::LocalSpace_PCC_2D::CreateReferenceElement(Polydim::PDETools::LocalSpace_PCC_2D::MethodTypes::FEM_PCC,
                                                                       method_order);

      const auto mesh_geometric_data = Polydim::PDETools::Mesh::PDE_Mesh_Utilities::compute_mesh_2D_geometry_data(geometry_utilities,
                                                                                                                mesh_utilities,
                                                                                                                mesh,
                                                                                                                Polydim::PDETools::LocalSpace_PCC_2D::MeshGeometricDataConfigiguration(reference_element_data));


      Polydim::PDETools::Mesh::MeshMatricesDAO_mesh_connectivity_data mesh_connectivity_data(mesh);
      Polydim::PDETools::DOFs::DOFsManager dofManager;

      std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> boundary_info =
      {{0, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
       {1, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
       {2, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
       {3, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 1}},
       {4, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
       {5, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
       {6, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 1}},
       {7, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 1}},
       {8, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}}};

      const auto mesh_dofs_info = Polydim::PDETools::LocalSpace_PCC_2D::SetMeshDOFsInfo(reference_element_data,
                                                                                      mesh,
                                                                                      boundary_info);
      const auto dofs_data = dofManager.CreateDOFs_2D(mesh_dofs_info,
                                                      mesh_connectivity_data);

      auto source_term = [&method_order](const Eigen::MatrixXd &points)
      {
        Eigen::VectorXd source_term = Eigen::VectorXd::Constant(points.cols(),
                                                                2.0 * method_order * (method_order - 1));
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + 0.5;

        const int max_order = method_order - 2;
        for (int i = 0; i < max_order; ++i)
            source_term.array() *= polynomial;

        return -source_term;
      };

      const auto forcing_term = PDETools::Assembler_Utilities::assembler_forcing_term(geometry_utilities,
                                                                                      mesh,
                                                                                      mesh_geometric_data,
                                                                                      mesh_dofs_info,
                                                                                      dofs_data,
                                                                                      reference_element_data,
                                                                                      source_term);

    }

  } // namespace UnitTesting
} // namespace Polydim

#endif
