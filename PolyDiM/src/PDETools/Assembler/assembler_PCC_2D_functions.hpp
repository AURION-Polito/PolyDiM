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
#include "Eigen_SparseArray.hpp"
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
      // ***************************************************************************
      struct Sparse_Matrix_Triplet final
      {
          unsigned int i;
          unsigned int j;
          double value;
      };
      // ***************************************************************************
      struct Sparse_Matrix_Data final
      {
          std::array<unsigned int, 2> size;
          std::vector<unsigned int> rows;
          std::vector<unsigned int> cols;
          std::vector<double> values;
      };
      // ***************************************************************************
      struct Variational_Operator final
      {
          Sparse_Matrix_Data A;
          Sparse_Matrix_Data A_Strong;
      };
      // ***************************************************************************
      std::list<Eigen::Triplet<double>> to_triplets(const Eigen::SparseMatrix<double> & M)
      {
          std::list<Eigen::Triplet<double>> v;

          for(int i = 0; i < M.outerSize(); i++)
          {
              for(typename Eigen::SparseMatrix<double>::InnerIterator it(M,i); it; ++it)
                  v.push_back(Eigen::Triplet<double>(it.row(),it.col(),it.value()));
          }

          return v;
      }
      // ***************************************************************************
      Gedim::Eigen_Array<> to_Eigen_Array(const Eigen::VectorXd& v)
      {
        return Gedim::Eigen_Array<>(v);
      }
      // ***************************************************************************
      Gedim::Eigen_SparseArray<> to_Eigen_SparseArray(const Sparse_Matrix_Data& A)
      {
        Gedim::Eigen_SparseArray<> eigen_A;
        eigen_A.SetSize(A.size.at(0),
                        A.size.at(1),
                        Gedim::ISparseArray::SparseArrayTypes::None);

        eigen_A.Triplets(A.rows,
                         A.cols,
                         A.values);
        return eigen_A;
      }
      // ***************************************************************************
      Sparse_Matrix_Data to_Sparse_Matrix_Data(const Gedim::Eigen_SparseArray<>& A)
      {
        Sparse_Matrix_Data result;

          const auto& eigen_A = static_cast<const Eigen::SparseMatrix<double>&>(A);
          const auto A_triplets = to_triplets(eigen_A);
          const auto num_triplets = A_triplets.size();

          result.size = { static_cast<unsigned int>(eigen_A.rows()),
                            static_cast<unsigned int>(eigen_A.cols()) };
          result.rows.resize(num_triplets);
          result.cols.resize(num_triplets);
          result.values.resize(num_triplets);

          unsigned int t = 0;
          for (const auto& triplet : A_triplets)
          {
            result.rows[t] = triplet.row();
            result.cols[t] = triplet.col();
            result.values[t] = triplet.value();
            t++;
          }

          return result;
      }
      // ***************************************************************************
      Variational_Operator convert_operator(const Gedim::Eigen_SparseArray<>& elliptic_matrix,
                                         const Gedim::Eigen_SparseArray<>& elliptic_strong_matrix)
      {
        Variational_Operator result;

        result.A = to_Sparse_Matrix_Data(elliptic_matrix);
        result.A_Strong = to_Sparse_Matrix_Data(elliptic_strong_matrix);

        return result;
      }
      // ***************************************************************************
      Eigen::VectorXd function_evaluation(const Eigen::MatrixXd &points,
                                          const std::function<double (const double&, const double&, const double&, const Eigen::VectorXd&)> f)
      {
        Eigen::VectorXd function_values(points.cols());

        for (int i = 0; i < points.cols(); ++i)
        {
          function_values[i] = f(points(0, i),
                             points(1, i),
                             points(2, i),
                             function_values);
        }

        return function_values;
      };
      // ***************************************************************************
      Eigen::VectorXd function_evaluation(const unsigned int marker,
                                          const Eigen::MatrixXd &points,
                                          const std::function<double (const unsigned int, const double&, const double&, const double&)> f)
      {
        Eigen::VectorXd function_values(points.cols());

        for (int i = 0; i < points.cols(); ++i)
        {
          function_values[i] = f(marker,
                                 points(0, i),
                                 points(1, i),
                                 points(2, i));
        }

        return function_values;
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

        Polydim::PDETools::Assembler_Utilities::local_matrix_to_global_matrix_dofs_data local_matrix_to_global_matrix_dofs_data =
        {{std::cref(dofs_data)}, {0}, {0}, {0}};

        for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); ++c)
        {
          if (!mesh.Cell2DIsActive(c))
            continue;

          const auto local_space_data = Polydim::PDETools::LocalSpace_PCC_2D::CreateLocalSpace(geometry_utilities.Tolerance1D(),
                                                                                               geometry_utilities.Tolerance2D(),
                                                                                               mesh_geometric_data,
                                                                                               c,
                                                                                               reference_element_data);

          const auto basis_functions_values =
              Polydim::PDETools::LocalSpace_PCC_2D::BasisFunctionsValues(reference_element_data, local_space_data);

          const auto cell2D_internal_quadrature =
              Polydim::PDETools::LocalSpace_PCC_2D::InternalQuadrature(reference_element_data, local_space_data);

          const Eigen::VectorXd source_term_values = function_evaluation(cell2D_internal_quadrature.Points,
                                                                            source_term_function);

          const Eigen::VectorXd local_rhs =
              equation.ComputeCellForcingTerm(source_term_values, basis_functions_values, cell2D_internal_quadrature.Weights);

          assert(Polydim::PDETools::LocalSpace_PCC_2D::Size(reference_element_data, local_space_data) == dofs_data.CellsGlobalDOFs[2].at(c).size());

          Polydim::PDETools::Assembler_Utilities::assemble_local_matrix_to_global_matrix<2>(c,
                                                                                            local_matrix_to_global_matrix_dofs_data,
                                                                                            local_rhs,
                                                                                            forcing_term);
        }

        return static_cast<Eigen::VectorXd&>(forcing_term);
      }
      // ***************************************************************************
      Variational_Operator assembler_elliptic_operator(
          const Gedim::GeometryUtilities& geometry_utilities,
          const Gedim::MeshMatricesDAO &mesh,
          const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
          const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
          const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
          const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
          const std::function<double (const double&, const double&, const double&, const Eigen::VectorXd&)> diffusion_term_function)
      {
        Gedim::Eigen_SparseArray<> elliptic_matrix;
        Gedim::Eigen_SparseArray<> elliptic_strong_matrix;

        elliptic_matrix.SetSize(dofs_data.NumberDOFs,
                                dofs_data.NumberDOFs,
                                Gedim::ISparseArray::SparseArrayTypes::None);
        elliptic_strong_matrix.SetSize(dofs_data.NumberDOFs,
                                       dofs_data.NumberStrongs);

        Polydim::PDETools::Equations::EllipticEquation equation;

        Polydim::PDETools::Assembler_Utilities::local_matrix_to_global_matrix_dofs_data local_matrix_to_global_matrix_dofs_data =
        {{std::cref(dofs_data)}, {0}, {0}, {0}};

        for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); ++c)
        {
          if (!mesh.Cell2DIsActive(c))
            continue;

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

          const Eigen::VectorXd diffusion_term_values = function_evaluation(cell2D_internal_quadrature.Points,
                                                                            diffusion_term_function);

          const Eigen::MatrixXd local_A = equation.ComputeCellDiffusionMatrix(diffusion_term_values,
                                                                              basis_functions_derivative_values,
                                                                              cell2D_internal_quadrature.Weights);

          const double k_max = diffusion_term_values.cwiseAbs().maxCoeff();
          const Eigen::MatrixXd local_A_stab =
              k_max *
              Polydim::PDETools::LocalSpace_PCC_2D::StabilizationMatrix(reference_element_data,
                                                                        local_space_data);


          assert(Polydim::PDETools::LocalSpace_PCC_2D::Size(reference_element_data, local_space_data) == dofs_data.CellsGlobalDOFs[2].at(c).size());


          Polydim::PDETools::Assembler_Utilities::assemble_local_matrix_to_global_matrix<2>(c,
                                                                                            local_matrix_to_global_matrix_dofs_data,
                                                                                            local_matrix_to_global_matrix_dofs_data,
                                                                                            local_A + local_A_stab,
                                                                                            elliptic_matrix,
                                                                                            elliptic_strong_matrix);

        }

        elliptic_matrix.Create();
        elliptic_strong_matrix.Create();

        return convert_operator(elliptic_matrix,
                                elliptic_strong_matrix);
      }
      // ***************************************************************************
      Eigen::VectorXd assembler_strong_solution(
          const Gedim::GeometryUtilities& geometry_utilities,
          const Gedim::MeshMatricesDAO &mesh,
          const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
          const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
          const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
          const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
          const std::function<double (const unsigned int, const double&, const double&, const double&)> strong_solution_function)
      {
        Gedim::Eigen_Array<> strong_solution;

        strong_solution.SetSize(dofs_data.NumberStrongs);

        Polydim::PDETools::Equations::EllipticEquation equation;

        Polydim::PDETools::Assembler_Utilities::local_matrix_to_global_matrix_dofs_data local_matrix_to_global_matrix_dofs_data =
        {{std::cref(dofs_data)}, {0}, {0}, {0}};

        for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); ++c)
        {
          if (!mesh.Cell2DIsActive(c))
            continue;

          const auto local_space_data = Polydim::PDETools::LocalSpace_PCC_2D::CreateLocalSpace(geometry_utilities.Tolerance1D(),
                                                                                               geometry_utilities.Tolerance2D(),
                                                                                               mesh_geometric_data,
                                                                                               c,
                                                                                               reference_element_data);

          const auto basis_functions_values =
              Polydim::PDETools::LocalSpace_PCC_2D::BasisFunctionsValues(reference_element_data, local_space_data);

          const auto cell2D_internal_quadrature =
              Polydim::PDETools::LocalSpace_PCC_2D::InternalQuadrature(reference_element_data, local_space_data);

          assert(Polydim::PDETools::LocalSpace_PCC_2D::Size(reference_element_data, local_space_data) == dofs_data.CellsGlobalDOFs[2].at(c).size());

          for (unsigned int v = 0; v < mesh.Cell2DNumberVertices(c); ++v)
          {
              const unsigned int cell0D_index = mesh.Cell2DVertex(c, v);
              const auto &boundary_info = mesh_dofs_info.CellsBoundaryInfo.at(0).at(cell0D_index);

              if (boundary_info.Type != Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::Strong)
                  continue;

              const auto coordinates = mesh.Cell0DCoordinates(cell0D_index);

              const auto strong_boundary_values = function_evaluation(boundary_info.Marker,
                                                                      coordinates,
                                                                      strong_solution_function);

              const auto local_dofs = dofs_data.CellsDOFs.at(0).at(cell0D_index);

              assert(local_dofs.size() == strong_boundary_values.size());

              for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
              {
                  const auto &local_dof_i = local_dofs.at(loc_i);

                  switch (local_dof_i.Type)
                  {
                  case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong: {
                      strong_solution.SetValue(local_dof_i.Global_Index, strong_boundary_values[loc_i]);
                  }
                  break;
                  case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                      continue;
                  default:
                      throw std::runtime_error("Unknown DOF Type");
                  }
              }
          }

          for (unsigned int ed = 0; ed < mesh.Cell2DNumberEdges(c); ++ed)
          {
              const unsigned int cell1D_index = mesh.Cell2DEdge(c, ed);

              const auto &boundary_info = mesh_dofs_info.CellsBoundaryInfo.at(1).at(cell1D_index);
              const auto local_dofs = dofs_data.CellsDOFs.at(1).at(cell1D_index);

              if (boundary_info.Type != Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::Strong ||
                  local_dofs.size() == 0)
                  continue;

              const auto edge_dofs_coordinates =
                  Polydim::PDETools::LocalSpace_PCC_2D::EdgeDofsCoordinates(reference_element_data, local_space_data, ed);

              const auto strong_boundary_values = function_evaluation(boundary_info.Marker,
                                                                      edge_dofs_coordinates,
                                                                      strong_solution_function);

              assert(local_dofs.size() == strong_boundary_values.size());

              for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
              {
                  const auto &local_dof_i = local_dofs.at(loc_i);

                  switch (local_dof_i.Type)
                  {
                  case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong: {
                      strong_solution.SetValue(local_dof_i.Global_Index, strong_boundary_values[loc_i]);
                  }
                  break;
                  case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                      continue;
                  default:
                      throw std::runtime_error("Unknown DOF Type");
                  }
              }
          }
        }

        return static_cast<Eigen::VectorXd&>(strong_solution);
      }
      // ***************************************************************************
    } // namespace Assembler_Utilities
  } // namespace PDETools
} // namespace Polydim

#endif
