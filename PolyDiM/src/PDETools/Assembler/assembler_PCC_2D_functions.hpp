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

#ifndef __PDETOOLS_ASSEMBLER_assembler_PCC_2D_functions_HPP
#define __PDETOOLS_ASSEMBLER_assembler_PCC_2D_functions_HPP

#include "Assembler_Utilities.hpp"
#include "Eigen_Array.hpp"
#include "EllipticEquation.hpp"
#include "MeshUtilities.hpp"
#include "DOFsManager.hpp"
#include "LocalSpace_PCC_2D.hpp"

namespace Polydim
{
namespace PDETools
{
namespace Assembler_Utilities
{
  Eigen::VectorXd source_term_evaluation(const Eigen::MatrixXd &points,
                                          const std::function<double (const double&, const double&, const double&, const Eigen::VectorXd&)> source_term_function)
  {
      Eigen::VectorXd source_term(points.cols());

      for (int i = 0; i < points.cols(); ++i)
      {
          source_term[i] = source_term_function(points(0, i),
                                                points(1, i),
                                                points(2, i),
                                                source_term);
      }

      return source_term;
  };

  // ***************************************************************************
  Eigen::VectorXd assembler_source_term(
      const Gedim::GeometryUtilities& geometry_utilities,
      const Gedim::MeshMatricesDAO &mesh,
      const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
      const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
      const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
      const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
      const std::function<double (const double&, const double&, const double&, const Eigen::VectorXd&)> source_term_function)
  {
      Gedim::Eigen_Array<> forcing_term;

      forcing_term.SetSize(dofs_data.NumberDOFs);

      Polydim::PDETools::Equations::EllipticEquation equation;

      for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); ++c)
      {
          const auto local_space_data = Polydim::PDETools::LocalSpace_PCC_2D::CreateLocalSpace(geometry_utilities.Tolerance1D(),
                                                                                               geometry_utilities.Tolerance2D(),
                                                                                               mesh_geometric_data,
                                                                                               c,
                                                                                               reference_element_data);

          const auto basis_functions_values =
              Polydim::PDETools::LocalSpace_PCC_2D::BasisFunctionsValues(reference_element_data, local_space_data);

          const auto basis_functions_derivative_values =
              Polydim::PDETools::LocalSpace_PCC_2D::BasisFunctionsDerivativeValues(reference_element_data, local_space_data);

          const auto cell2D_internal_quadrature =
              Polydim::PDETools::LocalSpace_PCC_2D::InternalQuadrature(reference_element_data, local_space_data);

          const Eigen::VectorXd source_term_values = source_term_evaluation(cell2D_internal_quadrature.Points,
                                                                            source_term_function);

          Eigen::VectorXd local_rhs =
              equation.ComputeCellForcingTerm(source_term_values, basis_functions_values, cell2D_internal_quadrature.Weights);

          const auto &global_dofs = dofs_data.CellsGlobalDOFs[2].at(c);

          assert(Polydim::PDETools::LocalSpace_PCC_2D::Size(reference_element_data, local_space_data) == global_dofs.size());

          Polydim::PDETools::Assembler_Utilities::local_matrix_to_global_matrix_dofs_data local_matrix_to_global_matrix_dofs_data =
              {{std::cref(dofs_data)}, {0}, {0}, {0}};

          Polydim::PDETools::Assembler_Utilities::assemble_local_matrix_to_global_matrix<2>(c,
                                                                                            local_matrix_to_global_matrix_dofs_data,
                                                                                            local_rhs,
                                                                                            forcing_term);
      }

      return static_cast<Eigen::VectorXd&>(forcing_term);
  }

// ***************************************************************************
} // namespace Assembler_Utilities
} // namespace PDETools
} // namespace Polydim

#endif
