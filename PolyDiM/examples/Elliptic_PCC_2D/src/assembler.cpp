#include "assembler.hpp"


#include "Quadrature_Gauss1D.hpp"

#include "VEM_PCC_2D_LocalSpace.hpp"
#include "EllipticEquation.hpp"

using namespace std;
using namespace Eigen;

#define DEBUG_ASSEMBLER 1
#undef DEBUG_ASSEMBLER

namespace Elliptic_PCC_2D
{
  // ***************************************************************************
  Assembler::Elliptic_PCC_2D_Problem_Data Assembler::Assemble(const Gedim::GeometryUtilities& geometryUtilities,
                                                              const Gedim::MeshMatricesDAO& mesh,
                                                              const Gedim::MeshUtilities::MeshGeometricData2D& mesh_geometric_data,
                                                              const Polydim::PDETools::DOFs::DOFsManager<2>::MeshDOFsInfo& mesh_dofs_info,
                                                              const Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData& dofs_data,
                                                              const Polydim::VEM::PCC::VEM_PCC_2D_ReferenceElement_Data& reference_element_data,
                                                              const std::function<Eigen::VectorXd(const Eigen::MatrixXd&)>& diffusion_term,
                                                              const std::function<Eigen::VectorXd(const Eigen::MatrixXd&)>& source_term,
                                                              const std::function<Eigen::VectorXd(const unsigned int,
                                                                                                  const Eigen::MatrixXd&)>& strong_boundary_condition,
                                                              const std::function<Eigen::VectorXd(const unsigned int,
                                                                                                  const Eigen::MatrixXd&)>& weak_boundary_condition) const
  {
    Elliptic_PCC_2D_Problem_Data result;

    result.globalMatrixA.SetSize(dofs_data.NumberDOFs,
                                 dofs_data.NumberDOFs,
                                 Gedim::ISparseArray::SparseArrayTypes::Symmetric);
    result.dirichletMatrixA.SetSize(dofs_data.NumberDOFs,
                                    dofs_data.NumberStrongs);
    result.rightHandSide.SetSize(dofs_data.NumberDOFs);
    result.solution.SetSize(dofs_data.NumberDOFs);
    result.solutionDirichlet.SetSize(dofs_data.NumberStrongs);

    Polydim::PDETools::Equations::EllipticEquation equation;

    for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); ++c)
    {
      const Polydim::VEM::PCC::VEM_PCC_2D_Polygon_Geometry polygon =
      {
        mesh_geometric_data.Cell2DsVertices.at(c),
        mesh_geometric_data.Cell2DsCentroids.at(c),
        mesh_geometric_data.Cell2DsAreas.at(c),
        mesh_geometric_data.Cell2DsDiameters.at(c),
        mesh_geometric_data.Cell2DsTriangulations.at(c),
        mesh_geometric_data.Cell2DsEdgeLengths.at(c),
        mesh_geometric_data.Cell2DsEdgeDirections.at(c),
        mesh_geometric_data.Cell2DsEdgeTangents.at(c),
        mesh_geometric_data.Cell2DsEdgeNormals.at(c)
      };

      Polydim::VEM::PCC::VEM_PCC_2D_LocalSpace vem_local_space;

      const auto local_space = vem_local_space.CreateLocalSpace(reference_element_data,
                                                                polygon);

      const auto basis_functions_values = vem_local_space.ComputeBasisFunctionsValues(local_space,
                                                                                      Polydim::VEM::PCC::ProjectionTypes::Pi0km1);

#ifdef DEBUG_ASSEMBLER
      std::cout.precision(2);
      std::cout<< std::scientific<< "u: "<< basis_functions_values<< std::endl;
#endif


      const auto basis_functions_derivative_values = vem_local_space.ComputeBasisFunctionsDerivativeValues(local_space,
                                                                                                           Polydim::VEM::PCC::ProjectionTypes::Pi0km1Der);
#ifdef DEBUG_ASSEMBLER
      std::cout<< std::scientific<< "du_x: "<< basis_functions_derivative_values[0]<< std::endl;
      std::cout<< std::scientific<< "du_y: "<< basis_functions_derivative_values[1]<< std::endl;
#endif

      const auto diffusion_term_values = diffusion_term(local_space.InternalQuadrature.Points);
      const auto source_term_values = source_term(local_space.InternalQuadrature.Points);

#ifdef DEBUG_ASSEMBLER
      std::cout<< std::scientific<< "k: "<< diffusion_term_values.transpose()<< std::endl;
      std::cout<< std::scientific<< "f: "<< source_term_values.transpose()<< std::endl;
#endif

      const auto local_A = equation.ComputeCellDiffusionMatrix(diffusion_term_values,
                                                               basis_functions_derivative_values,
                                                               local_space.InternalQuadrature.Weights);
#ifdef DEBUG_ASSEMBLER
      std::cout<< std::scientific<< "A: "<< local_A<< std::endl;
#endif

      const auto local_stab_A = diffusion_term_values.cwiseAbs().maxCoeff() *
                                local_space.StabMatrix;

#ifdef DEBUG_ASSEMBLER
      std::cout<< std::scientific<< "S: "<< local_stab_A<< std::endl;
#endif

      const auto local_rhs = equation.ComputeCellForcingTerm(source_term_values,
                                                             basis_functions_values,
                                                             local_space.InternalQuadrature.Weights);

#ifdef DEBUG_ASSEMBLER
      std::cout<< std::scientific<< "rhs: "<< local_rhs.transpose()<< std::endl;
#endif

      const auto& global_dofs = dofs_data.CellsGlobalDOFs[2].at(c);

      assert(local_space.NumBasisFunctions ==  global_dofs.size());

      for (unsigned int loc_i = 0; loc_i < global_dofs.size(); ++loc_i)
      {
        const auto& global_dof_i = global_dofs.at(loc_i);
        const auto& local_dof_i = dofs_data.CellsDOFs.at(global_dof_i.Dimension).at(global_dof_i.CellIndex).at(global_dof_i.DOFIndex);

        switch (local_dof_i.Type)
        {
          case Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData::DOF::Types::Strong:
            continue;
          case Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData::DOF::Types::DOF:
            break;
          default:
            throw std::runtime_error("Unknown DOF Type");
        }

        const unsigned int global_index_i = local_dof_i.Global_Index;

        result.rightHandSide.AddValue(global_index_i,
                                      local_rhs[loc_i]);

#ifdef DEBUG_ASSEMBLER
        std::cout<< std::scientific<< "g_rhs["<< global_index_i<< "] = "<< local_rhs[loc_i]<< std::endl;
#endif

        for (unsigned int loc_j = 0; loc_j < global_dofs.size(); ++loc_j)
        {
          const auto& global_dof_j = global_dofs.at(loc_j);
          const auto& local_dof_j = dofs_data.CellsDOFs.at(global_dof_j.Dimension).at(global_dof_j.CellIndex).at(global_dof_j.DOFIndex);

          const unsigned int global_index_j = local_dof_j.Global_Index;
          const double loc_A_element =  local_A(loc_i,
                                                loc_j) +
                                        local_stab_A(loc_i,
                                                     loc_j);

          switch (local_dof_j.Type)
          {
            case Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData::DOF::Types::Strong:
              result.dirichletMatrixA.Triplet(global_index_i,
                                              global_index_j,
                                              loc_A_element);

#ifdef DEBUG_ASSEMBLER
              std::cout<< std::scientific<< "g_AD("<< global_index_i<< ","<< global_index_j<< ") = "<< loc_A_element<< std::endl;
#endif
              break;
            case Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData::DOF::Types::DOF:
              result.globalMatrixA.Triplet(global_index_i,
                                           global_index_j,
                                           loc_A_element);

#ifdef DEBUG_ASSEMBLER
              std::cout<< std::scientific<< "g_A("<< global_index_i<< ","<< global_index_j<< ") = "<< loc_A_element<< std::endl;
#endif
              break;
            default:
              throw std::runtime_error("Unknown DOF Type");
          }
        }
      }

      ComputeWeakTerm(c,
                      mesh,
                      polygon,
                      mesh_dofs_info,
                      dofs_data,
                      reference_element_data,
                      weak_boundary_condition,
                      result);
    }

    ComputeStrongTerm(geometryUtilities,
                      mesh,
                      mesh_geometric_data,
                      mesh_dofs_info,
                      dofs_data,
                      reference_element_data,
                      strong_boundary_condition,
                      result);

    result.rightHandSide.Create();
    result.solutionDirichlet.Create();
    result.globalMatrixA.Create();
    result.dirichletMatrixA.Create();

    if (dofs_data.NumberStrongs > 0)
      result.rightHandSide.SubtractionMultiplication(result.dirichletMatrixA,
                                                     result.solutionDirichlet);

    return result;
  }
  // ***************************************************************************
  void Assembler::ComputeStrongTerm(const Gedim::GeometryUtilities& geometryUtilities,
                                    const Gedim::MeshMatricesDAO& mesh,
                                    const Gedim::MeshUtilities::MeshGeometricData2D& mesh_geometric_data,
                                    const Polydim::PDETools::DOFs::DOFsManager<2>::MeshDOFsInfo& mesh_dofs_info,
                                    const Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData& dofs_data,
                                    const Polydim::VEM::PCC::VEM_PCC_2D_ReferenceElement_Data& reference_element_data,
                                    const std::function<Eigen::VectorXd(const unsigned int,
                                                                        const Eigen::MatrixXd&)>& strong_boundary_condition,
                                    Elliptic_PCC_2D_Problem_Data& assembler_data) const
  {
    // Assemble strong boundary condition on Cell0Ds
    for (unsigned int p = 0; p < mesh.Cell0DTotalNumber(); ++p)
    {
      const auto& boundary_info = mesh_dofs_info.CellsBoundaryInfo.at(0).at(p);

      if (boundary_info.Type !=
          Polydim::PDETools::DOFs::DOFsManager<2>::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::Strong)
        continue;

      const auto coordinates = mesh.Cell0DCoordinates(p);

      const auto strong_boundary_values = strong_boundary_condition(boundary_info.Marker,
                                                                    coordinates);

      const auto local_dofs = dofs_data.CellsDOFs.at(0).at(p);

      assert(local_dofs.size() == strong_boundary_values.size());

      for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
      {
        const auto& local_dof_i = local_dofs.at(loc_i);

        switch (local_dof_i.Type)
        {
          case Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData::DOF::Types::Strong:
          {
            assembler_data.solutionDirichlet.SetValue(local_dof_i.Global_Index,
                                                      strong_boundary_values[loc_i]);
          }
            break;
          case Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData::DOF::Types::DOF:
            continue;
          default:
            throw std::runtime_error("Unknown DOF Type");
        }
      }
    }

    // Assemble strong boundary condition on Cell1Ds
    const auto& referenceSegmentInternalPoints = reference_element_data.Quadrature.ReferenceSegmentInternalPoints;
    const unsigned int numReferenceSegmentInternalPoints = referenceSegmentInternalPoints.cols();

    for (unsigned int e = 0; e < mesh.Cell1DTotalNumber(); ++e)
    {
      const auto& boundary_info = mesh_dofs_info.CellsBoundaryInfo.at(1).at(e);

      if (boundary_info.Type !=
          Polydim::PDETools::DOFs::DOFsManager<2>::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::Strong)
        continue;

      const auto cell1D_origin = mesh.Cell1DOriginCoordinates(e);
      const auto cell1D_end = mesh.Cell1DEndCoordinates(e);
      const auto cell1D_tangent = geometryUtilities.SegmentTangent(cell1D_origin,
                                                                   cell1D_end);

      Eigen::MatrixXd coordinates = Eigen::MatrixXd::Zero(3, numReferenceSegmentInternalPoints);
      for (unsigned int r = 0; r < numReferenceSegmentInternalPoints; r++)
        coordinates.col(r)<< cell1D_origin +
                             referenceSegmentInternalPoints(0, r) * cell1D_tangent;

      const auto strong_boundary_values = strong_boundary_condition(boundary_info.Marker,
                                                                    coordinates);

      const auto local_dofs = dofs_data.CellsDOFs.at(1).at(e);

      assert(local_dofs.size() == strong_boundary_values.size());

      for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
      {
        const auto& local_dof_i = local_dofs.at(loc_i);

        switch (local_dof_i.Type)
        {
          case Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData::DOF::Types::Strong:
          {
            assembler_data.solutionDirichlet.SetValue(local_dof_i.Global_Index,
                                                      strong_boundary_values[loc_i]);
          }
            break;
          case Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData::DOF::Types::DOF:
            continue;
          default:
            throw std::runtime_error("Unknown DOF Type");
        }
      }
    }
  }
  // ***************************************************************************
  void Assembler::ComputeWeakTerm(const unsigned int cell2DIndex,
                                  const Gedim::MeshMatricesDAO& mesh,
                                  const Polydim::VEM::PCC::VEM_PCC_2D_Polygon_Geometry& polygon,
                                  const Polydim::PDETools::DOFs::DOFsManager<2>::MeshDOFsInfo& mesh_dofs_info,
                                  const Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData& dofs_data,
                                  const Polydim::VEM::PCC::VEM_PCC_2D_ReferenceElement_Data& reference_element_data,
                                  const std::function<Eigen::VectorXd(const unsigned int,
                                                                      const Eigen::MatrixXd&)>& weak_boundary_condition,
                                  Elliptic_PCC_2D_Problem_Data& assembler_data) const
  {
    const unsigned numVertices = polygon.Vertices.cols();

    for(unsigned int ed = 0; ed < numVertices; ed ++)
    {
      const unsigned int cell1D_index = mesh.Cell2DEdge(cell2DIndex,
                                                        ed);

      const auto& boundary_info = mesh_dofs_info.CellsBoundaryInfo.at(1).at(cell1D_index);

      if (boundary_info.Type !=
          Polydim::PDETools::DOFs::DOFsManager<2>::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::Weak)
        continue;

      Polydim::VEM::PCC::VEM_PCC_2D_LocalSpace vem_local_space;

      // compute vem values
      const auto weakReferenceSegment = Gedim::Quadrature::Quadrature_Gauss1D::FillPointsAndWeights(2 * reference_element_data.Order);

      const VectorXd pointsCurvilinearCoordinates = weakReferenceSegment.Points.row(0).transpose();

      const auto weak_basis_function_values = vem_local_space.ComputeValuesOnEdge(reference_element_data,
                                                                                  pointsCurvilinearCoordinates);

      // map edge internal quadrature points
      const Vector3d& edgeStart = polygon.EdgesDirection[ed] ?
                                    polygon.Vertices.col(ed) :
                                    polygon.Vertices.col((ed + 1) % numVertices);

      const Vector3d& edgeTangent = polygon.EdgesTangent.col(ed);
      const double direction = polygon.EdgesDirection[ed] ? 1.0 : -1.0;

      const unsigned int numEdgeWeakQuadraturePoints = weakReferenceSegment.Points.cols();
      MatrixXd weakQuadraturePoints(3, numEdgeWeakQuadraturePoints);
      for (unsigned int q = 0; q < numEdgeWeakQuadraturePoints; q++)
      {
        weakQuadraturePoints.col(q) = edgeStart + direction *
                                      weakReferenceSegment.Points(0, q) *
                                      edgeTangent;
      }
      const double absMapDeterminant = std::abs(polygon.EdgesLength[ed]);
      const MatrixXd weakQuadratureWeights = weakReferenceSegment.Weights *
                                             absMapDeterminant;

      const VectorXd neumannValues = weak_boundary_condition(boundary_info.Marker,
                                                             weakQuadraturePoints);

      // compute values of Neumann condition
      const VectorXd neumannContributions = weak_basis_function_values.transpose() *
                                            weakQuadratureWeights.asDiagonal() *
                                            neumannValues;

#ifdef DEBUG_ASSEMBLER
      std::cout.precision(2);
      std::cout<< std::scientific<< "u_n: "<< neumannContributions.transpose()<< std::endl;
#endif

      for (unsigned int p = 0; p < 2; ++p)
      {
        const unsigned int cell0D_index = mesh.Cell1DVertex(cell1D_index,
                                                            p);

        const auto local_dofs = dofs_data.CellsDOFs.at(0).at(cell0D_index);

        for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
        {
          const auto& local_dof_i = local_dofs.at(loc_i);

          switch (local_dof_i.Type)
          {
            case Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData::DOF::Types::Strong:
              continue;
            case Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData::DOF::Types::DOF:
            {
              assembler_data.rightHandSide.AddValue(local_dof_i.Global_Index,
                                                    neumannContributions[p]);

#ifdef DEBUG_ASSEMBLER
              std::cout<< std::scientific<< "rhs["<< local_dof_i.Global_Index<< "]: "<< neumannContributions(p)<< " -> "<< assembler_data.rightHandSide.GetValue(local_dof_i.Global_Index)<< std::endl;
#endif
            }
              break;
            default:
              throw std::runtime_error("Unknown DOF Type");
          }
        }
      }

      const auto local_dofs = dofs_data.CellsDOFs.at(1).at(cell1D_index);

      for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
      {
        const auto& local_dof_i = local_dofs.at(loc_i);

        const unsigned int localIndex = polygon.EdgesDirection[ed] ? loc_i :
                                                                     local_dofs.size() - 1 - loc_i;


        switch (local_dof_i.Type)
        {
          case Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData::DOF::Types::Strong:
            continue;
          case Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData::DOF::Types::DOF:
          {
            assembler_data.rightHandSide.AddValue(local_dof_i.Global_Index,
                                                  neumannContributions[localIndex + 2]);

#ifdef DEBUG_ASSEMBLER
            std::cout<< std::scientific<< "rhs_ed["<< local_dof_i.Global_Index<< "]: "<< neumannContributions(localIndex + 2)<< " -> "<< assembler_data.rightHandSide.GetValue(local_dof_i.Global_Index)<< std::endl;
#endif
          }
            break;
          default:
            throw std::runtime_error("Unknown DOF Type");
        }
      }
    }
  }
  // ***************************************************************************
  Assembler::VEM_Performance_Result Assembler::ComputeVemPerformance(const Gedim::GeometryUtilities& geometryUtilities,
                                                                     const Gedim::MeshMatricesDAO& mesh,
                                                                     const Gedim::MeshUtilities::MeshGeometricData2D& mesh_geometric_data,
                                                                     const Polydim::VEM::PCC::VEM_PCC_2D_ReferenceElement_Data& reference_element_data) const
  {
    Assembler::VEM_Performance_Result result;
    result.Cell2DsPerformance.resize(mesh.Cell2DTotalNumber());

    // Assemble equation elements
    for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); c++)
    {
      const Polydim::VEM::PCC::VEM_PCC_2D_Polygon_Geometry polygon =
      {
        mesh_geometric_data.Cell2DsVertices.at(c),
        mesh_geometric_data.Cell2DsCentroids.at(c),
        mesh_geometric_data.Cell2DsAreas.at(c),
        mesh_geometric_data.Cell2DsDiameters.at(c),
        mesh_geometric_data.Cell2DsTriangulations.at(c),
        mesh_geometric_data.Cell2DsEdgeLengths.at(c),
        mesh_geometric_data.Cell2DsEdgeDirections.at(c),
        mesh_geometric_data.Cell2DsEdgeTangents.at(c),
        mesh_geometric_data.Cell2DsEdgeNormals.at(c)
      };

      Polydim::VEM::PCC::VEM_PCC_2D_LocalSpace vem_local_space;

      const auto local_space = vem_local_space.CreateLocalSpace(reference_element_data,
                                                                polygon);

      Polydim::VEM::PCC::VEM_PCC_PerformanceAnalysis performanceAnalysis;

      result.Cell2DsPerformance[c].Analysis = performanceAnalysis.Compute(polygon.Measure,
                                                                          polygon.Diameter,
                                                                          Polydim::VEM::Monomials::VEM_Monomials_2D(),
                                                                          reference_element_data.Monomials,
                                                                          vem_local_space,
                                                                          local_space);

      result.Cell2DsPerformance[c].NumInternalQuadraturePoints = local_space.InternalQuadrature.Weights.size();
      result.Cell2DsPerformance[c].NumBoundaryQuadraturePoints = local_space.BoundaryQuadrature.Quadrature.Weights.size();
    }

    return result;
  }
  // ***************************************************************************
  Assembler::PostProcess_Data Assembler::PostProcessSolution(const Gedim::GeometryUtilities& geometryUtilities,
                                                             const Gedim::MeshMatricesDAO& mesh,
                                                             const Gedim::MeshUtilities::MeshGeometricData2D& mesh_geometric_data,
                                                             const Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData& dofs_data,
                                                             const Polydim::VEM::PCC::VEM_PCC_2D_ReferenceElement_Data& reference_element_data,
                                                             const Elliptic_PCC_2D_Problem_Data& assembler_data,
                                                             const std::function<Eigen::VectorXd(const Eigen::MatrixXd&)>& exact_solution,
                                                             const std::function<std::array<Eigen::VectorXd, 3>(const Eigen::MatrixXd&)>& exact_derivative_solution) const
  {
    PostProcess_Data result;

    result.residual_norm = 0.0;
    if (dofs_data.NumberDOFs > 0)
    {
      Gedim::Eigen_Array<> residual;
      residual.SetSize(dofs_data.NumberDOFs);
      residual.SumMultiplication(assembler_data.globalMatrixA,
                                 assembler_data.solution);
      residual -= assembler_data.rightHandSide;

      result.residual_norm = residual.Norm();
    }

    result.cell0Ds_numeric.setZero(mesh.Cell0DTotalNumber());
    result.cell0Ds_exact.setZero(mesh.Cell0DTotalNumber());

    for (unsigned int p = 0; p < mesh.Cell0DTotalNumber(); p++)
    {
      result.cell0Ds_exact[p] = exact_solution(mesh.Cell0DCoordinates(p))[0];

      const auto local_dofs = dofs_data.CellsDOFs.at(0).at(p);

      for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
      {
        const auto& local_dof_i = local_dofs.at(loc_i);

        switch (local_dof_i.Type)
        {
          case Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData::DOF::Types::Strong:
            result.cell0Ds_numeric[p] = assembler_data.solutionDirichlet.GetValue(local_dof_i.Global_Index);
            break;
          case Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData::DOF::Types::DOF:
            result.cell0Ds_numeric[p] = assembler_data.solution.GetValue(local_dof_i.Global_Index);
            break;
          default:
            throw std::runtime_error("Unknown DOF Type");
        }
      }
    }

    result.cell2Ds_error_L2.setZero(mesh.Cell2DTotalNumber());
    result.cell2Ds_norm_L2.setZero(mesh.Cell2DTotalNumber());
    result.cell2Ds_error_H1.setZero(mesh.Cell2DTotalNumber());
    result.cell2Ds_norm_H1.setZero(mesh.Cell2DTotalNumber());
    result.error_L2 = 0.0;
    result.norm_L2 = 0.0;
    result.error_H1 = 0.0;
    result.norm_H1 = 0.0;
    result.mesh_size = 0.0;

    for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); c++)
    {
      const Polydim::VEM::PCC::VEM_PCC_2D_Polygon_Geometry polygon =
      {
        mesh_geometric_data.Cell2DsVertices.at(c),
        mesh_geometric_data.Cell2DsCentroids.at(c),
        mesh_geometric_data.Cell2DsAreas.at(c),
        mesh_geometric_data.Cell2DsDiameters.at(c),
        mesh_geometric_data.Cell2DsTriangulations.at(c),
        mesh_geometric_data.Cell2DsEdgeLengths.at(c),
        mesh_geometric_data.Cell2DsEdgeDirections.at(c),
        mesh_geometric_data.Cell2DsEdgeTangents.at(c),
        mesh_geometric_data.Cell2DsEdgeNormals.at(c)
      };

      Polydim::VEM::PCC::VEM_PCC_2D_LocalSpace vem_local_space;

      const auto local_space = vem_local_space.CreateLocalSpace(reference_element_data,
                                                                polygon);

      const auto basis_functions_values = vem_local_space.ComputeBasisFunctionsValues(local_space,
                                                                                      Polydim::VEM::PCC::ProjectionTypes::Pi0k);

      const auto basis_functions_derivative_values = vem_local_space.ComputeBasisFunctionsDerivativeValues(local_space,
                                                                                                           Polydim::VEM::PCC::ProjectionTypes::Pi0km1Der);


      const auto exact_solution_values = exact_solution(local_space.InternalQuadrature.Points);
      const auto exact_derivative_solution_values = exact_derivative_solution(local_space.InternalQuadrature.Points);

      const auto& global_dofs = dofs_data.CellsGlobalDOFs[2].at(c);
      Eigen::VectorXd dofs_values = Eigen::VectorXd::Zero(global_dofs.size());

      for (unsigned int loc_i = 0; loc_i < global_dofs.size(); ++loc_i)
      {
        const auto& global_dof_i = global_dofs.at(loc_i);
        const auto& local_dof_i = dofs_data.CellsDOFs.at(global_dof_i.Dimension).at(global_dof_i.CellIndex).at(global_dof_i.DOFIndex);

        switch (local_dof_i.Type)
        {
          case Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData::DOF::Types::Strong:
            dofs_values[loc_i] = assembler_data.solutionDirichlet.GetValue(local_dof_i.Global_Index);
            break;
          case Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData::DOF::Types::DOF:
            dofs_values[loc_i] = assembler_data.solution.GetValue(local_dof_i.Global_Index);
            break;
          default:
            throw std::runtime_error("Unknown DOF Type");
        }
      }

      const Eigen::VectorXd local_error_L2 = (basis_functions_values * dofs_values -
                                              exact_solution_values).array().square();
      const Eigen::VectorXd local_norm_L2 = (basis_functions_values * dofs_values).array().square();

      result.cell2Ds_error_L2[c] = local_space.InternalQuadrature.Weights.transpose() *
                                   local_error_L2;
      result.cell2Ds_norm_L2[c] = local_space.InternalQuadrature.Weights.transpose() *
                                  local_norm_L2;


      const Eigen::VectorXd local_error_H1 =
          (basis_functions_derivative_values[0] * dofs_values -
          exact_derivative_solution_values[0]).array().square() +
          (basis_functions_derivative_values[1] * dofs_values -
          exact_derivative_solution_values[1]).array().square();
      const Eigen::VectorXd local_norm_H1 =
          (basis_functions_derivative_values[0] * dofs_values).array().square() +
          (basis_functions_derivative_values[1] * dofs_values).array().square();

      result.cell2Ds_error_H1[c] = local_space.InternalQuadrature.Weights.transpose() *
                                   local_error_H1;
      result.cell2Ds_norm_H1[c] = local_space.InternalQuadrature.Weights.transpose() *
                                  local_norm_H1;

      if (mesh_geometric_data.Cell2DsDiameters.at(c) > result.mesh_size)
        result.mesh_size = mesh_geometric_data.Cell2DsDiameters.at(c);
    }

    result.error_L2 = std::sqrt(result.cell2Ds_error_L2.sum());
    result.norm_L2 = std::sqrt(result.cell2Ds_norm_L2.sum());
    result.error_H1 = std::sqrt(result.cell2Ds_error_H1.sum());
    result.norm_H1 = std::sqrt(result.cell2Ds_norm_H1.sum());

    return result;
  }
  // ***************************************************************************
}
