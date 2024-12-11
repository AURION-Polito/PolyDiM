#include "assembler.hpp"

#include "ranges"

#include "VEM_MCC_2D_VelocityLocalSpace.hpp"
#include "EllipticEquation.hpp"
#include "VEM_MCC_VelocityLocalSpace_Data.hpp"

using namespace std;
using namespace Eigen;


namespace Elliptic_MCC_2D
{
// ***************************************************************************
Assembler::Elliptic_MCC_2D_Problem_Data Assembler::Assemble(const Gedim::GeometryUtilities& geometryUtilities,
                                                            const Gedim::MeshMatricesDAO& mesh,
                                                            const Gedim::MeshUtilities::MeshGeometricData2D& mesh_geometric_data,
                                                            const std::vector<Polydim::PDETools::DOFs::DOFsManager<2>::MeshDOFsInfo>& mesh_dofs_info,
                                                            const std::vector<Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData>& dofs_data,
                                                            const Polydim::VEM::MCC::VEM_MCC_2D_Velocity_ReferenceElement_Data& velocity_reference_element_data,
                                                            const Polydim::VEM::MCC::VEM_MCC_2D_Pressure_ReferenceElement_Data& pressure_reference_element_data,
                                                            const std::function<std::array<Eigen::VectorXd, 3>(const Eigen::MatrixXd&)>& advection_term,
                                                            const std::function<Eigen::VectorXd(const Eigen::MatrixXd&)>& reaction_term,
                                                            const std::function<std::array<Eigen::VectorXd, 9>(const Eigen::MatrixXd&)>& diffusion_term,
                                                            const std::function<Eigen::VectorXd(const Eigen::MatrixXd&)>& source_term,
                                                            const std::function<Eigen::VectorXd(const unsigned int,
                                                                                                const Eigen::MatrixXd&)>& weak_boundary_condition) const
{
    Elliptic_MCC_2D_Problem_Data result;
    const unsigned int numDOFHandler = mesh_dofs_info.size();
    unsigned int numberDOFs = 0;
    unsigned int numberStrongs = 0;
    std::vector<unsigned int> offsetDOFs = {0};
    std::vector<unsigned int> offsetStrongs = {0};
    for(unsigned int i = 0; i < numDOFHandler; i++)
    {
        numberDOFs += dofs_data[i].NumberDOFs;
        numberStrongs += dofs_data[i].NumberStrongs;
    }

    result.globalMatrixA.SetSize(numberDOFs, numberDOFs,
                                 Gedim::ISparseArray::SparseArrayTypes::Symmetric);
    result.neumannMatrixA.SetSize(numberDOFs,
                                  numberStrongs);
    result.rightHandSide.SetSize(numberDOFs);
    result.solution.SetSize(numberDOFs);
    result.solutionNeumann.SetSize(numberStrongs);

    Polydim::PDETools::Equations::EllipticEquation equation;

    for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); c++)
    {
        const Polydim::VEM::MCC::VEM_MCC_2D_Polygon_Geometry polygon =
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

        Polydim::VEM::MCC::VEM_MCC_2D_VelocityLocalSpace vem_local_space;

        const auto local_space = vem_local_space.CreateLocalSpace(velocity_reference_element_data,
                                                                  polygon);

        const auto velocity_basis_functions_values = vem_local_space.ComputeBasisFunctionsValues(local_space);
        const auto velocity_basis_functions_divergence_values = vem_local_space.ComputeBasisFunctionsDivergenceValues(local_space);
        const auto pressure_basis_functions_values = monomials.Vander(pressure_reference_element_data.Monomials,
                                                                      local_space.InternalQuadrature.Points,
                                                                      mesh_geometric_data.Cell2DsCentroids.at(c),
                                                                      mesh_geometric_data.Cell2DsDiameters.at(c));

        const auto reaction_term_values = reaction_term(local_space.InternalQuadrature.Points);
        const auto advection_term_values = advection_term(local_space.InternalQuadrature.Points);
        const auto diffusion_term_values = diffusion_term(local_space.InternalQuadrature.Points);
        const auto source_term_values = source_term(local_space.InternalQuadrature.Points);

        auto local_A = equation.ComputeCellDiffusionMatrix(diffusion_term_values,
                                                           velocity_basis_functions_values,
                                                           local_space.InternalQuadrature.Weights);

        double kmax = 0.0;
        std::ranges::for_each(diffusion_term_values |
                                  std::ranges::views::transform([](const auto m){ return m.cwiseAbs().maxCoeff(); }),
                              [&kmax](const auto m){ kmax = kmax < m ? m : kmax;});

        local_A += kmax * local_space.StabMatrix;

        const auto local_M = equation.ComputeCellReactionMatrix(reaction_term_values,
                                                                pressure_basis_functions_values,
                                                                local_space.InternalQuadrature.Weights);

        const auto local_T = equation.ComputeCellAdvectionMatrix(advection_term_values,
                                                                 pressure_basis_functions_values,
                                                                 velocity_basis_functions_values,
                                                                 local_space.InternalQuadrature.Weights);

        const Eigen::MatrixXd local_B = pressure_basis_functions_values.transpose()
                                        * local_space.InternalQuadrature.Weights.asDiagonal()
                                        * velocity_basis_functions_divergence_values;

        const auto local_rhs = equation.ComputeCellForcingTerm(source_term_values,
                                                               pressure_basis_functions_values,
                                                               local_space.InternalQuadrature.Weights);

        const std::vector<size_t> offset_global_dofs = {0lu, dofs_data[0].CellsGlobalDOFs[2].at(c).size()};
        const std::vector<std::vector<Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData::GlobalCell_DOF>> global_dofs
            = {dofs_data[0].CellsGlobalDOFs[2].at(c), dofs_data[1].CellsGlobalDOFs[2].at(c)};

        Eigen::MatrixXd elemental_matrix = MatrixXd::Zero(global_dofs[0].size() + global_dofs[1].size(), global_dofs[0].size() + global_dofs[1].size());
        Eigen::VectorXd elemental_rhs = VectorXd::Zero(global_dofs[0].size() + global_dofs[1].size());
        elemental_matrix << local_A, -local_B.transpose()+local_T,
            local_B, local_M;
        elemental_rhs <<  VectorXd::Zero(global_dofs[0].size()), local_rhs;

        assert(local_space.NumBasisFunctions ==  global_dofs[0].size());

        for(unsigned h1 = 0; h1 < numDOFHandler; h1++)
        {
            for (unsigned int loc_i = 0; loc_i < global_dofs[h1].size(); loc_i++)
            {
                const auto global_dof_i = global_dofs[h1].at(loc_i);
                const auto local_dof_i = dofs_data[h1].CellsDOFs.at(global_dof_i.Dimension).at(global_dof_i.CellIndex).at(global_dof_i.DOFIndex);

                switch (local_dof_i.Type)
                {
                case Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData::DOF::Types::Strong:
                    continue;
                case Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData::DOF::Types::DOF:
                    break;
                default:
                    throw std::runtime_error("Unknown DOF Type");
                }

                const unsigned int global_index_i = local_dof_i.Global_Index + offsetDOFs[h1];

                result.rightHandSide.AddValue(global_index_i,
                                              elemental_rhs[loc_i + offset_global_dofs[h1]]);


                for(unsigned int h2 = 0; h2 < numDOFHandler; h2++)
                {
                    for (unsigned int loc_j = 0; loc_j < global_dofs[h2].size(); loc_j++)
                    {
                        const auto& global_dof_j = global_dofs[h2].at(loc_j);
                        const auto& local_dof_j = dofs_data[h2].CellsDOFs.at(global_dof_j.Dimension).at(global_dof_j.CellIndex).at(global_dof_j.DOFIndex);

                        const unsigned int global_index_j = local_dof_j.Global_Index;
                        const double loc_A_element =  elemental_matrix(loc_i + offset_global_dofs[h1], loc_j + offset_global_dofs[h2]);

                        switch (local_dof_j.Type)
                        {
                        case Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData::DOF::Types::Strong:
                            result.neumannMatrixA.Triplet(global_index_i,
                                                          global_index_j + offsetStrongs[h2],
                                                          loc_A_element);
                            break;
                        case Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData::DOF::Types::DOF:
                            result.globalMatrixA.Triplet(global_index_i,
                                                         global_index_j + offsetDOFs[h2],
                                                         loc_A_element);
                            break;
                        default:
                            throw std::runtime_error("Unknown DOF Type");
                        }
                    }
                }
            }
        }
    }

    //    ComputeStrongTerm(geometryUtilities,
    //                      mesh,
    //                      mesh_geometric_data,
    //                      mesh_dofs_info,
    //                      dofs_data,
    //                      reference_element_data,
    //                      strong_boundary_condition,
    //                      result);

    result.rightHandSide.Create();
    result.solutionNeumann.Create();
    result.globalMatrixA.Create();
    result.neumannMatrixA.Create();

    if (numberStrongs > 0)
        result.rightHandSide.SubtractionMultiplication(result.neumannMatrixA,
                                                       result.solutionNeumann);

    return result;
}
// ***************************************************************************
//void Assembler::ComputeStrongTerm(const Gedim::GeometryUtilities& geometryUtilities,
//                                  const Gedim::MeshMatricesDAO& mesh,
//                                  const Gedim::MeshUtilities::MeshGeometricData2D& mesh_geometric_data,
//                                  const Polydim::PDETools::DOFs::DOFsManager<2>::MeshDOFsInfo& mesh_dofs_info,
//                                  const Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData& dofs_data,
//                                  const Polydim::VEM::MCC::VEM_MCC_2D_Velocity_ReferenceElement_Data& reference_element_data,
//                                  const std::function<Eigen::VectorXd(const unsigned int,
//                                                                      const Eigen::MatrixXd&)>& strong_boundary_condition,
//                                  Elliptic_MCC_2D_Problem_Data& assembler_data) const
//{
//    // Assemble strong boundary condition on Cell0Ds
//    for (unsigned int p = 0; p < mesh.Cell0DTotalNumber(); ++p)
//    {
//        const auto& boundary_info = mesh_dofs_info.CellsBoundaryInfo.at(0).at(p);

//        if (boundary_info.Type !=
//            Polydim::PDETools::DOFs::DOFsManager<2>::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::Strong)
//            continue;

//        const auto coordinates = mesh.Cell0DCoordinates(p);

//        const auto strong_boundary_values = strong_boundary_condition(boundary_info.Marker,
//                                                                      coordinates);

//        const auto local_dofs = dofs_data.CellsDOFs.at(0).at(p);

//        assert(local_dofs.size() == strong_boundary_values.size());

//        for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
//        {
//            const auto& local_dof_i = local_dofs.at(loc_i);

//            switch (local_dof_i.Type)
//            {
//            case Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData::DOF::Types::Strong:
//            {
//                assembler_data.solutionNeumann.SetValue(local_dof_i.Global_Index,
//                                                        strong_boundary_values[loc_i]);
//            }
//            break;
//            case Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData::DOF::Types::DOF:
//                continue;
//            default:
//                throw std::runtime_error("Unknown DOF Type");
//            }
//        }
//    }

//    // Assemble strong boundary condition on Cell1Ds
//    const auto& referenceSegmentInternalPoints = reference_element_data.Quadrature.ReferenceSegmentInternalPoints;
//    const unsigned int numReferenceSegmentInternalPoints = referenceSegmentInternalPoints.cols();

//    for (unsigned int e = 0; e < mesh.Cell1DTotalNumber(); ++e)
//    {
//        const auto& boundary_info = mesh_dofs_info.CellsBoundaryInfo.at(1).at(e);

//        if (boundary_info.Type !=
//            Polydim::PDETools::DOFs::DOFsManager<2>::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::Strong)
//            continue;

//        const auto cell1D_origin = mesh.Cell1DOriginCoordinates(e);
//        const auto cell1D_end = mesh.Cell1DEndCoordinates(e);
//        const auto cell1D_tangent = geometryUtilities.SegmentTangent(cell1D_origin,
//                                                                     cell1D_end);

//        Eigen::MatrixXd coordinates = Eigen::MatrixXd::Zero(3, numReferenceSegmentInternalPoints);
//        for (unsigned int r = 0; r < numReferenceSegmentInternalPoints; r++)
//            coordinates.col(r)<< cell1D_origin +
//                                      referenceSegmentInternalPoints(0, r) * cell1D_tangent;

//        const auto strong_boundary_values = strong_boundary_condition(boundary_info.Marker,
//                                                                      coordinates);

//        const auto local_dofs = dofs_data.CellsDOFs.at(1).at(e);

//        assert(local_dofs.size() == strong_boundary_values.size());

//        for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
//        {
//            const auto& local_dof_i = local_dofs.at(loc_i);

//            switch (local_dof_i.Type)
//            {
//            case Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData::DOF::Types::Strong:
//            {
//                assembler_data.solutionNeumann.SetValue(local_dof_i.Global_Index,
//                                                        strong_boundary_values[loc_i]);
//            }
//            break;
//            case Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData::DOF::Types::DOF:
//                continue;
//            default:
//                throw std::runtime_error("Unknown DOF Type");
//            }
//        }
//    }
//}
//// ***************************************************************************
////  void Assembler::ComputeWeakTerm(const Gedim::IMeshDAO& mesh,
////                                  const unsigned int& cell2DIndex,
////                                  const VectorXd& cell2DEdgeLengths,
////                                  const MatrixXd& cell2DEdgeTangents,
////                                  const MatrixXd& cell2DEdgeNormals,
////                                  const MatrixXd& cell2DVertices,
////                                  const vector<bool>& cell2DEdgeDirections,
////                                  const Gedim::IDOFManagement& dofManager,
////                                  const Gedim::VEM_IValues_MCC_2D& vemValues,
////                                  const Gedim::VEM_ValuesData& vemLocalSpace,
////                                  const Gedim::VEM_IQuadrature2D& vemQuadrature,
////                                  const Gedim::IWeakBoundaryCondition& weakBoundaryCondition,
////                                  Gedim::IArray& rightHandSide) const
////  {
////    if (dofManager.NumberLocalWeaks(cell2DIndex) == 0)
////      return;

////    const unsigned numVertices = cell2DVertices.cols();

////    for(unsigned int ed = 0; ed < numVertices; ed ++)
////    {
////      const unsigned int edgeMeshIndex = mesh.Cell2DEdge(cell2DIndex,
////                                                         ed);

////      if (!dofManager.IsCellWeakBoundaryCondition(edgeMeshIndex, 1))
////        continue;

////      const unsigned int edgeMarker = dofManager.CellMarker(edgeMeshIndex, 1);

////      // compute vem values
////      const unsigned int quadratureOrder = max(weakBoundaryCondition.QuadratureOrder(edgeMarker),
////                                               2 * vemLocalSpace.Order);

////      Eigen::MatrixXd weakReferenceSegmentPoints;
////      Eigen::VectorXd weakReferenceSegmentWeights;
////      Gedim::Quadrature_Gauss1D::FillPointsAndWeights(quadratureOrder,
////                                                      weakReferenceSegmentPoints,
////                                                      weakReferenceSegmentWeights);

////      const VectorXd pointsCurvilinearCoordinates = weakReferenceSegmentPoints.row(0).transpose();

////      VectorXd edgeInternalPoints; // edge DOF: Gauss Lobatto quadrature points
////      const MatrixXd vemRefSegmentInternalPoints = vemQuadrature.ReferenceSegmentInternalPoints();
////      if (vemRefSegmentInternalPoints.rows() > 0)
////        edgeInternalPoints = vemRefSegmentInternalPoints.row(0).transpose();

////      const MatrixXd weakBasisFunctionsValues = vemValues.ComputeValuesOnEdge(edgeInternalPoints,
////                                                                              pointsCurvilinearCoordinates);

////      // map edge internal quadrature points
////      const Vector3d& edgeStart = cell2DEdgeDirections[ed] ? cell2DVertices.col(ed) :
////                                                             cell2DVertices.col((ed + 1) % numVertices);
////      const Vector3d& edgeTangent = cell2DEdgeTangents.col(ed);
////      const double direction = cell2DEdgeDirections[ed] ? 1.0 : -1.0;

////      const unsigned int numEdgeWeakQuadraturePoints = weakReferenceSegmentPoints.cols();
////      MatrixXd weakQuadraturePoints(3, numEdgeWeakQuadraturePoints);
////      for (unsigned int q = 0; q < numEdgeWeakQuadraturePoints; q++)
////      {
////        weakQuadraturePoints.col(q) = edgeStart + direction *
////                                      weakReferenceSegmentPoints(0, q) *
////                                      edgeTangent;
////      }
////      const double absMapDeterminant = std::abs(cell2DEdgeLengths[ed]);
////      const MatrixXd weakQuadratureWeights = weakReferenceSegmentWeights *
////                                             absMapDeterminant;

////      const VectorXd neumannValues = weakBoundaryCondition.Evaluate(edgeMarker,
////                                                                    weakQuadraturePoints);

////      // compute values of Neumann condition
////      const VectorXd neumannContributions = weakBasisFunctionsValues.transpose() *
////                                            weakQuadratureWeights.asDiagonal() *
////                                            neumannValues;

////      // add contributions relative to edge extrema.
////      for(unsigned int p = 0; p < 2; ++p)
////      {
////        const unsigned int vertexGlobalIndex = mesh.Cell1DVertex(edgeMeshIndex,
////                                                                 p);

////        //        if (dofManager.CellMarker(vertexGlobalIndex, 0) != edgeMarker)
////        //          continue;

////        const unsigned int numCell0DLocals = dofManager.NumberLocals(vertexGlobalIndex,
////                                                                     0);
////        for (unsigned int l = 0; l < numCell0DLocals; l++)
////        {
////          //                    if (!dofManager.IsWeakBoundaryCondition(vertexGlobalIndex, l, 0))
////          //                      continue;
////          if (dofManager.IsStrongBoundaryCondition(vertexGlobalIndex, l, 0))
////            continue;

////          const int globalNeumann_i = dofManager.GlobalIndex(vertexGlobalIndex,
////                                                             l,
////                                                             0);
////          rightHandSide.AddValue(globalNeumann_i, neumannContributions(p));
////        }
////      }

////      const unsigned int numCell1DLocals = dofManager.NumberLocals(edgeMeshIndex,
////                                                                   1);

////      // add contributions relative to edge internal dofs (here we assume we are dealing with
////      // conforming VEM)
////      for (unsigned int l = 0; l < numCell1DLocals; l++)
////      {
////        const unsigned int localIndex = cell2DEdgeDirections[ed] ? l :
////                                                                   numCell1DLocals - 1 - l;

////        const int globalNeumann_i = dofManager.GlobalIndex(edgeMeshIndex, localIndex, 1);

////        rightHandSide.AddValue(globalNeumann_i, neumannContributions(localIndex + 2));
////      }
////    }
////  }
// ***************************************************************************
Assembler::VEM_Performance_Result Assembler::ComputeVemPerformance(const Gedim::GeometryUtilities& geometryUtilities,
                                                                   const Gedim::MeshMatricesDAO& mesh,
                                                                   const Gedim::MeshUtilities::MeshGeometricData2D& mesh_geometric_data,
                                                                   const Polydim::VEM::MCC::VEM_MCC_2D_Velocity_ReferenceElement_Data& reference_element_data) const
{
    Assembler::VEM_Performance_Result result;
    result.Cell2DsPerformance.resize(mesh.Cell2DTotalNumber());

    // Assemble equation elements
    for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); c++)
    {
        const Polydim::VEM::MCC::VEM_MCC_2D_Polygon_Geometry polygon =
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

        Polydim::VEM::MCC::VEM_MCC_2D_VelocityLocalSpace vem_local_space;

        const auto local_space = vem_local_space.CreateLocalSpace(reference_element_data,
                                                                  polygon);

        Polydim::VEM::MCC::VEM_MCC_PerformanceAnalysis performanceAnalysis;

        result.Cell2DsPerformance[c].Analysis = performanceAnalysis.Compute(polygon.Measure,
                                                                            polygon.Diameter,
                                                                            Polydim::VEM::Monomials::VEM_Monomials_2D(),
                                                                            reference_element_data.MonomialsKp1,
                                                                            vem_local_space,
                                                                            local_space);

        result.Cell2DsPerformance[c].NumInternalQuadraturePoints = local_space.InternalQuadrature.Weights.size();
        result.Cell2DsPerformance[c].NumBoundaryQuadraturePoints = local_space.BoundaryQuadrature.Quadrature.Weights.size();
    }

    return result;
}
// ***************************************************************************
//Assembler::PostProcess_Data Assembler::PostProcessSolution(const Gedim::GeometryUtilities& geometryUtilities,
//                                                           const Gedim::MeshMatricesDAO& mesh,
//                                                           const Gedim::MeshUtilities::MeshGeometricData2D& mesh_geometric_data,
//                                                           const Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData& dofs_data,
//                                                           const Polydim::VEM::MCC::VEM_MCC_2D_ReferenceElement_Data& reference_element_data,
//                                                           const Elliptic_MCC_2D_Problem_Data& assembler_data,
//                                                           const std::function<Eigen::VectorXd(const Eigen::MatrixXd&)>& exact_solution,
//                                                           const std::function<std::array<Eigen::VectorXd, 3>(const Eigen::MatrixXd&)>& exact_derivative_solution) const
//{
//    PostProcess_Data result;

//    result.residual_norm = 0.0;
//    if (dofs_data.NumberDOFs > 0)
//    {
//        Gedim::Eigen_Array<> residual;
//        residual.SetSize(dofs_data.NumberDOFs);
//        residual.SumMultiplication(assembler_data.globalMatrixA,
//                                   assembler_data.solution);
//        residual -= assembler_data.rightHandSide;

//        result.residual_norm = residual.Norm();
//    }

//    result.cell0Ds_numeric.setZero(mesh.Cell0DTotalNumber());
//    result.cell0Ds_exact.setZero(mesh.Cell0DTotalNumber());

//    for (unsigned int p = 0; p < mesh.Cell0DTotalNumber(); p++)
//    {
//        result.cell0Ds_exact[p] = exact_solution(mesh.Cell0DCoordinates(p))[0];

//        const auto local_dofs = dofs_data.CellsDOFs.at(0).at(p);

//        for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
//        {
//            const auto& local_dof_i = local_dofs.at(loc_i);

//            switch (local_dof_i.Type)
//            {
//            case Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData::DOF::Types::Strong:
//                result.cell0Ds_numeric[p] = assembler_data.solutionNeumann.GetValue(local_dof_i.Global_Index);
//                break;
//            case Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData::DOF::Types::DOF:
//                result.cell0Ds_numeric[p] = assembler_data.solution.GetValue(local_dof_i.Global_Index);
//                break;
//            default:
//                throw std::runtime_error("Unknown DOF Type");
//            }
//        }
//    }

//    result.cell2Ds_error_L2.setZero(mesh.Cell2DTotalNumber());
//    result.cell2Ds_norm_L2.setZero(mesh.Cell2DTotalNumber());
//    result.cell2Ds_error_H1.setZero(mesh.Cell2DTotalNumber());
//    result.cell2Ds_norm_H1.setZero(mesh.Cell2DTotalNumber());
//    result.error_L2 = 0.0;
//    result.norm_L2 = 0.0;
//    result.error_H1 = 0.0;
//    result.norm_H1 = 0.0;
//    result.mesh_size = 0.0;

//    for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); c++)
//    {
//        const Polydim::VEM::MCC::VEM_MCC_2D_Polygon_Geometry polygon =
//            {
//                mesh_geometric_data.Cell2DsVertices.at(c),
//                mesh_geometric_data.Cell2DsCentroids.at(c),
//                mesh_geometric_data.Cell2DsAreas.at(c),
//                mesh_geometric_data.Cell2DsDiameters.at(c),
//                mesh_geometric_data.Cell2DsTriangulations.at(c),
//                mesh_geometric_data.Cell2DsEdgeLengths.at(c),
//                mesh_geometric_data.Cell2DsEdgeDirections.at(c),
//                mesh_geometric_data.Cell2DsEdgeTangents.at(c),
//                mesh_geometric_data.Cell2DsEdgeNormals.at(c)
//            };

//        Polydim::VEM::MCC::VEM_MCC_2D_LocalSpace vem_local_space;

//        const auto local_space = vem_local_space.CreateLocalSpace(reference_element_data,
//                                                                  polygon);

//        const auto basis_functions_values = vem_local_space.ComputeBasisFunctionsValues(local_space,
//                                                                                        Polydim::VEM::MCC::ProjectionTypes::Pi0k);

//        const auto basis_functions_derivative_values = vem_local_space.ComputeBasisFunctionsDerivativeValues(local_space,
//                                                                                                             Polydim::VEM::MCC::ProjectionTypes::Pi0km1Der);


//        const auto exact_solution_values = exact_solution(local_space.InternalQuadrature.Points);
//        const auto exact_derivative_solution_values = exact_derivative_solution(local_space.InternalQuadrature.Points);

//        const auto& global_dofs = dofs_data.CellsGlobalDOFs[2].at(c);
//        Eigen::VectorXd dofs_values = Eigen::VectorXd::Zero(global_dofs.size());

//        for (unsigned int loc_i = 0; loc_i < global_dofs.size(); ++loc_i)
//        {
//            const auto& global_dof_i = global_dofs.at(loc_i);
//            const auto& local_dof_i = dofs_data.CellsDOFs.at(global_dof_i.Dimension).at(global_dof_i.CellIndex).at(global_dof_i.DOFIndex);

//            switch (local_dof_i.Type)
//            {
//            case Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData::DOF::Types::Strong:
//                dofs_values[loc_i] = assembler_data.solutionNeumann.GetValue(local_dof_i.Global_Index);
//                break;
//            case Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData::DOF::Types::DOF:
//                dofs_values[loc_i] = assembler_data.solution.GetValue(local_dof_i.Global_Index);
//                break;
//            default:
//                throw std::runtime_error("Unknown DOF Type");
//            }
//        }

//        const Eigen::VectorXd local_error_L2 = (basis_functions_values * dofs_values -
//                                                exact_solution_values).array().square();
//        const Eigen::VectorXd local_norm_L2 = (basis_functions_values * dofs_values).array().square();

//        result.cell2Ds_error_L2[c] = local_space.InternalQuadrature.Weights.transpose() *
//                                     local_error_L2;
//        result.cell2Ds_norm_L2[c] = local_space.InternalQuadrature.Weights.transpose() *
//                                    local_norm_L2;


//        const Eigen::VectorXd local_error_H1 =
//            (basis_functions_derivative_values[0] * dofs_values -
//             exact_derivative_solution_values[0]).array().square() +
//            (basis_functions_derivative_values[1] * dofs_values -
//             exact_derivative_solution_values[1]).array().square();
//        const Eigen::VectorXd local_norm_H1 =
//            (basis_functions_derivative_values[0] * dofs_values).array().square() +
//            (basis_functions_derivative_values[1] * dofs_values).array().square();

//        result.cell2Ds_error_H1[c] = local_space.InternalQuadrature.Weights.transpose() *
//                                     local_error_H1;
//        result.cell2Ds_norm_H1[c] = local_space.InternalQuadrature.Weights.transpose() *
//                                    local_norm_H1;

//        if (mesh_geometric_data.Cell2DsDiameters.at(c) > result.mesh_size)
//            result.mesh_size = mesh_geometric_data.Cell2DsDiameters.at(c);
//    }

//    result.error_L2 = std::sqrt(result.cell2Ds_error_L2.sum());
//    result.norm_L2 = std::sqrt(result.cell2Ds_norm_L2.sum());
//    result.error_H1 = std::sqrt(result.cell2Ds_error_H1.sum());
//    result.norm_H1 = std::sqrt(result.cell2Ds_norm_H1.sum());

//    return result;
//}
// ***************************************************************************
}
