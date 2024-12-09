#include "assembler.hpp"


#include "Quadrature_Gauss1D.hpp"

#include "VEM_PCC_2D_LocalSpace.hpp"
#include "EllipticEquation.hpp"

using namespace std;
using namespace Eigen;

namespace Elliptic_PCC_2D
{
  // ***************************************************************************
  Assembler::Elliptic_PCC_2D_Problem_Data Assembler::Assemble(const Gedim::GeometryUtilities& geometryUtilities,
                                                              const Gedim::MeshMatricesDAO& mesh,
                                                              const Gedim::MeshUtilities::MeshGeometricData2D& mesh_geometric_data,
                                                              const Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData& dofs_data,
                                                              const Polydim::VEM::PCC::VEM_PCC_2D_ReferenceElement_Data& reference_element_data,
                                                              const std::function<Eigen::VectorXd(const Eigen::MatrixXd&)>& diffusionTerm,
                                                              const std::function<Eigen::VectorXd(const Eigen::MatrixXd&)>& sourceTerm) const
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
                                                                                      Polydim::VEM::PCC::ProjectionTypes::Pi0k);


      const auto basis_functions_derivative_values = vem_local_space.ComputeBasisFunctionsDerivativeValues(local_space,
                                                                                                           Polydim::VEM::PCC::ProjectionTypes::Pi0km1Der);


      const auto diffusion_term_values = diffusionTerm(local_space.InternalQuadrature.Points);
      const auto source_term_values = sourceTerm(local_space.InternalQuadrature.Points);


      const auto local_A = equation.ComputeCellDiffusionMatrix(diffusion_term_values,
                                                               basis_functions_derivative_values,
                                                               local_space.InternalQuadrature.Weights);

      const auto local_stab_A = diffusion_term_values.cwiseAbs().maxCoeff() *
                                local_space.StabMatrix;

      const auto local_rhs = equation.ComputeCellForcingTerm(source_term_values,
                                                             basis_functions_values,
                                                             local_space.InternalQuadrature.Weights);

      const auto& global_dofs = dofs_data.CellsGlobalDOFs[2].at(c);

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

        for (unsigned int loc_j = 0; loc_j < global_dofs.size(); ++loc_j)
        {
          const auto& global_dof_j = global_dofs.at(loc_j);
          const auto& local_dof_j = dofs_data.CellsDOFs.at(global_dof_j.Dimension).at(global_dof_j.CellIndex).at(global_dof_j.DOFIndex);

          const unsigned int global_index_j = local_dof_j.Global_Index;
          const double loc_A_element =  local_A(loc_i,
                                                loc_j) +
                                        local_stab_A(loc_i,
                                                     loc_j);

          switch (local_dof_i.Type)
          {
            case Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData::DOF::Types::Strong:
              result.dirichletMatrixA.Triplet(global_index_i,
                                              global_index_j,
                                              loc_A_element);
              break;
            case Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData::DOF::Types::DOF:
              result.globalMatrixA.Triplet(global_index_i,
                                           global_index_j,
                                           loc_A_element);
              break;
            default:
              throw std::runtime_error("Unknown DOF Type");
          }
        }

      }
    }

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
  //  void Assembler::ComputeStrongTerm(const Gedim::GeometryUtilities geometryUtilities,
  //                                    const Gedim::IMeshDAO& mesh,
  //                                    const Gedim::IDOFManagement& dofManager,
  //                                    const Gedim::IStrongBoundaryCondition& strongBoundaryCondition,
  //                                    Gedim::IArray& solutionDirichlet) const
  //  {
  //    // Assemble strong boundary condition on Cell0Ds
  //    for (unsigned int p = 0; p < mesh.Cell0DTotalNumber(); p++)
  //    {
  //      if (!mesh.Cell0DIsActive(p))
  //        continue;

  //      if (!dofManager.IsCellStrongBoundaryCondition(p, 0))
  //        continue;

  //      const Vector3d coordinates = mesh.Cell0DCoordinates(p);

  //      const unsigned int numCell0DLocals = dofManager.NumberLocals(p, 0);

  //      for (unsigned int l = 0; l < numCell0DLocals; l++)
  //      {
  //        if (!dofManager.IsStrongBoundaryCondition(p, l, 0))
  //          continue;

  //        const unsigned int dofMarker = dofManager.Marker(p, l, 0);
  //        const int globalDirichlet_i = dofManager.PartialGlobalIndex(p, l, 0);

  //        VectorXd dirichletTermValue = strongBoundaryCondition.Evaluate(dofMarker,
  //                                                                       coordinates);

  //        solutionDirichlet.SetValue(globalDirichlet_i,
  //                                   dirichletTermValue(0));
  //      }
  //    }

  //    // Assemble strong boundary condition on Cell1Ds
  //    Eigen::MatrixXd referenceSegmentInternalPoints = vemQuadrature.ReferenceSegmentInternalPoints();
  //    const unsigned int numReferenceSegmentInternalPoints = referenceSegmentInternalPoints.cols();
  //    for (unsigned int e = 0; e < mesh.Cell1DTotalNumber(); e++)
  //    {
  //      if (!mesh.Cell1DIsActive(e))
  //        continue;

  //      if (!dofManager.IsCellStrongBoundaryCondition(e, 1))
  //        continue;

  //      const unsigned int edgeMarker = dofManager.CellMarker(e, 1);

  //      const Eigen::Vector3d cell1DOrigin = mesh.Cell1DOriginCoordinates(e);
  //      const Eigen::Vector3d cell1DEnd = mesh.Cell1DEndCoordinates(e);
  //      const Eigen::Vector3d cell1DTangent = geometryUtilities.SegmentTangent(cell1DOrigin,
  //                                                                             cell1DEnd);
  //      MatrixXd coordinates;
  //      coordinates.setZero(3, numReferenceSegmentInternalPoints);
  //      for (unsigned int r = 0; r < numReferenceSegmentInternalPoints; r++)
  //        coordinates.col(r)<< cell1DOrigin + referenceSegmentInternalPoints(0, r) * cell1DTangent;

  //      const unsigned int numCell1DLocals = dofManager.NumberLocals(e, 1);

  //      for (unsigned int l = 0; l < numCell1DLocals; l++)
  //      {
  //        if (!dofManager.IsStrongBoundaryCondition(e, l, 1))
  //          continue;

  //        const int globalDirichlet_i = dofManager.PartialGlobalIndex(e, l, 1);

  //        VectorXd dirichletTermValue = strongBoundaryCondition.Evaluate(edgeMarker,
  //                                                                       coordinates.col(l));

  //        solutionDirichlet.SetValue(globalDirichlet_i,
  //                                   dirichletTermValue(0));
  //      }
  //    }
  //  }
  // ***************************************************************************
  //  void Assembler::ComputeWeakTerm(const Gedim::IMeshDAO& mesh,
  //                                  const unsigned int& cell2DIndex,
  //                                  const VectorXd& cell2DEdgeLengths,
  //                                  const MatrixXd& cell2DEdgeTangents,
  //                                  const MatrixXd& cell2DEdgeNormals,
  //                                  const MatrixXd& cell2DVertices,
  //                                  const vector<bool>& cell2DEdgeDirections,
  //                                  const Gedim::IDOFManagement& dofManager,
  //                                  const Gedim::VEM_IValues_PCC_2D& vemValues,
  //                                  const Gedim::VEM_ValuesData& vemLocalSpace,
  //                                  const Gedim::VEM_IQuadrature2D& vemQuadrature,
  //                                  const Gedim::IWeakBoundaryCondition& weakBoundaryCondition,
  //                                  Gedim::IArray& rightHandSide) const
  //  {
  //    if (dofManager.NumberLocalWeaks(cell2DIndex) == 0)
  //      return;

  //    const unsigned numVertices = cell2DVertices.cols();

  //    for(unsigned int ed = 0; ed < numVertices; ed ++)
  //    {
  //      const unsigned int edgeMeshIndex = mesh.Cell2DEdge(cell2DIndex,
  //                                                         ed);

  //      if (!dofManager.IsCellWeakBoundaryCondition(edgeMeshIndex, 1))
  //        continue;

  //      const unsigned int edgeMarker = dofManager.CellMarker(edgeMeshIndex, 1);

  //      // compute vem values
  //      const unsigned int quadratureOrder = max(weakBoundaryCondition.QuadratureOrder(edgeMarker),
  //                                               2 * vemLocalSpace.Order);

  //      Eigen::MatrixXd weakReferenceSegmentPoints;
  //      Eigen::VectorXd weakReferenceSegmentWeights;
  //      Gedim::Quadrature_Gauss1D::FillPointsAndWeights(quadratureOrder,
  //                                                      weakReferenceSegmentPoints,
  //                                                      weakReferenceSegmentWeights);

  //      const VectorXd pointsCurvilinearCoordinates = weakReferenceSegmentPoints.row(0).transpose();

  //      VectorXd edgeInternalPoints; // edge DOF: Gauss Lobatto quadrature points
  //      const MatrixXd vemRefSegmentInternalPoints = vemQuadrature.ReferenceSegmentInternalPoints();
  //      if (vemRefSegmentInternalPoints.rows() > 0)
  //        edgeInternalPoints = vemRefSegmentInternalPoints.row(0).transpose();

  //      const MatrixXd weakBasisFunctionsValues = vemValues.ComputeValuesOnEdge(edgeInternalPoints,
  //                                                                              pointsCurvilinearCoordinates);

  //      // map edge internal quadrature points
  //      const Vector3d& edgeStart = cell2DEdgeDirections[ed] ? cell2DVertices.col(ed) :
  //                                                             cell2DVertices.col((ed + 1) % numVertices);
  //      const Vector3d& edgeTangent = cell2DEdgeTangents.col(ed);
  //      const double direction = cell2DEdgeDirections[ed] ? 1.0 : -1.0;

  //      const unsigned int numEdgeWeakQuadraturePoints = weakReferenceSegmentPoints.cols();
  //      MatrixXd weakQuadraturePoints(3, numEdgeWeakQuadraturePoints);
  //      for (unsigned int q = 0; q < numEdgeWeakQuadraturePoints; q++)
  //      {
  //        weakQuadraturePoints.col(q) = edgeStart + direction *
  //                                      weakReferenceSegmentPoints(0, q) *
  //                                      edgeTangent;
  //      }
  //      const double absMapDeterminant = std::abs(cell2DEdgeLengths[ed]);
  //      const MatrixXd weakQuadratureWeights = weakReferenceSegmentWeights *
  //                                             absMapDeterminant;

  //      const VectorXd neumannValues = weakBoundaryCondition.Evaluate(edgeMarker,
  //                                                                    weakQuadraturePoints);

  //      // compute values of Neumann condition
  //      const VectorXd neumannContributions = weakBasisFunctionsValues.transpose() *
  //                                            weakQuadratureWeights.asDiagonal() *
  //                                            neumannValues;

  //      // add contributions relative to edge extrema.
  //      for(unsigned int p = 0; p < 2; ++p)
  //      {
  //        const unsigned int vertexGlobalIndex = mesh.Cell1DVertex(edgeMeshIndex,
  //                                                                 p);

  //        //        if (dofManager.CellMarker(vertexGlobalIndex, 0) != edgeMarker)
  //        //          continue;

  //        const unsigned int numCell0DLocals = dofManager.NumberLocals(vertexGlobalIndex,
  //                                                                     0);
  //        for (unsigned int l = 0; l < numCell0DLocals; l++)
  //        {
  //          //                    if (!dofManager.IsWeakBoundaryCondition(vertexGlobalIndex, l, 0))
  //          //                      continue;
  //          if (dofManager.IsStrongBoundaryCondition(vertexGlobalIndex, l, 0))
  //            continue;

  //          const int globalNeumann_i = dofManager.GlobalIndex(vertexGlobalIndex,
  //                                                             l,
  //                                                             0);
  //          rightHandSide.AddValue(globalNeumann_i, neumannContributions(p));
  //        }
  //      }

  //      const unsigned int numCell1DLocals = dofManager.NumberLocals(edgeMeshIndex,
  //                                                                   1);

  //      // add contributions relative to edge internal dofs (here we assume we are dealing with
  //      // conforming VEM)
  //      for (unsigned int l = 0; l < numCell1DLocals; l++)
  //      {
  //        const unsigned int localIndex = cell2DEdgeDirections[ed] ? l :
  //                                                                   numCell1DLocals - 1 - l;

  //        const int globalNeumann_i = dofManager.GlobalIndex(edgeMeshIndex, localIndex, 1);

  //        rightHandSide.AddValue(globalNeumann_i, neumannContributions(localIndex + 2));
  //      }
  //    }
  //  }
  // ***************************************************************************
  Assembler::VEM_Performance_Result Assembler::ComputeVemPerformance(const Gedim::GeometryUtilities& geometryUtilities,
                                                                     const Gedim::IMeshDAO& mesh,
                                                                     const vector<Eigen::MatrixXd>& meshCell2DsVertices,
                                                                     const vector<vector<bool>>& meshCell2DsEdgeDirections,
                                                                     const vector<vector<Eigen::Matrix3d> >& meshCell2DsTriangulations,
                                                                     const vector<double>& meshCell2DsAreas,
                                                                     const vector<Eigen::Vector3d>& meshCell2DsCentroids,
                                                                     const vector<double>& meshCell2DsDiameters,
                                                                     const vector<Eigen::VectorXd>& meshCell2DsEdgeLengths,
                                                                     const vector<Eigen::MatrixXd>& meshCell2DsEdgeTangents,
                                                                     const vector<Eigen::MatrixXd>& meshCell2DsEdgeNormals) const
  {
    //    Assembler::VEM_Performance_Result result;
    //    result.Cell2DsPerformance.resize(mesh.Cell2DTotalNumber());

    //    Gedim::VEM_PerformanceAnalysis performanceAnalysis(vemMonomials,
    //                                                       vemValues);

    //    // Assemble equation elements
    //    for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); c++)
    //    {
    //      if (!mesh.Cell2DIsActive(c))
    //        continue;

    //      const MatrixXd& cell2DVertices = meshCell2DsVertices.at(c);

    //      const vector<bool>& cell2DEdgeDirections = meshCell2DsEdgeDirections.at(c);
    //      const vector<Matrix3d>& cell2DTriangulationPoints = meshCell2DsTriangulations.at(c);
    //      const double& cell2DArea = meshCell2DsAreas.at(c);
    //      const Vector3d& cell2DCentroid = meshCell2DsCentroids.at(c);
    //      const double& cell2DDiameter = meshCell2DsDiameters.at(c);
    //      const VectorXd& cell2DEdgeLengths = meshCell2DsEdgeLengths.at(c);
    //      const MatrixXd& cell2DEdgeTangents = meshCell2DsEdgeTangents.at(c);
    //      const MatrixXd& cell2DEdgeNormals = meshCell2DsEdgeNormals.at(c);

    //      // Compute Internal Quadrature points
    //      MatrixXd internalQuadraturePoints2D;
    //      VectorXd internalQuadratureWeights2D;
    //      vemQuadrature.PolygonInternalQuadrature(cell2DVertices,
    //                                              cell2DTriangulationPoints,
    //                                              internalQuadraturePoints2D,
    //                                              internalQuadratureWeights2D);

    //      // Compute Boundary Quadrature points
    //      Eigen::MatrixXd boundaryQuadraturePoints1D;
    //      Eigen::VectorXd boundaryQuadratureWeights1D;
    //      std::vector<Eigen::VectorXd> boundaryQuadratureWeightsTimesNormal1D;
    //      vemQuadrature.PolygonEdgesQuadrature(cell2DVertices,
    //                                           cell2DEdgeLengths,
    //                                           cell2DEdgeDirections,
    //                                           cell2DEdgeTangents,
    //                                           cell2DEdgeNormals,
    //                                           boundaryQuadraturePoints1D,
    //                                           boundaryQuadratureWeights1D,
    //                                           boundaryQuadratureWeightsTimesNormal1D);

    //      Gedim::VEM_ValuesData localSpace = vemValues.CreateLocalSpace(cell2DVertices,
    //                                                                    cell2DCentroid,
    //                                                                    cell2DArea,
    //                                                                    cell2DDiameter,
    //                                                                    internalQuadraturePoints2D,
    //                                                                    internalQuadratureWeights2D,
    //                                                                    boundaryQuadraturePoints1D,
    //                                                                    boundaryQuadratureWeights1D,
    //                                                                    boundaryQuadratureWeightsTimesNormal1D,
    //                                                                    Gedim::VEM_ValuesData::ProjectionTypes::Pi0km1);

    //      result.Cell2DsPerformance[c].NumInternalQuadraturePoints = internalQuadratureWeights2D.size();
    //      result.Cell2DsPerformance[c].NumBoundaryQuadraturePoints = boundaryQuadratureWeights1D.size();
    //      result.Cell2DsPerformance[c].Analysis = performanceAnalysis.Compute(cell2DArea,
    //                                                                          cell2DDiameter,
    //                                                                          localSpace);
    //    }

    //return result;
  }
  // ***************************************************************************
}
