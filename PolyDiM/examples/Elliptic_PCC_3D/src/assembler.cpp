#include "assembler.hpp"

namespace Polydim
{
namespace examples
{
namespace Elliptic_PCC_3D
{
//***************************************************************************
void Assembler::ComputeStrongTerm(const unsigned int& cell3DIndex,
                                  const Gedim::MeshMatricesDAO& mesh,
                                  const std::vector<Polydim::VEM::PCC::VEM_PCC_2D_Polygon_Geometry> &polygonalFaces,
                                  const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo& mesh_dofs_info,
                                  const Polydim::PDETools::DOFs::DOFsManager::DOFsData& dofs_data,
                                  const Polydim::VEM::PCC::VEM_PCC_2D_ReferenceElement_Data& reference_element_data_2D,
                                  const Polydim::VEM::PCC::VEM_PCC_3D_ReferenceElement_Data& reference_element_data_3D,
                                  const Polydim::VEM::PCC::VEM_PCC_3D_LocalSpace_Data& local_space_data,
                                  const test::I_Test& test,
                                  Elliptic_PCC_3D_Problem_Data& assembler_data) const
{
    // Assemble strong boundary condition on Cell0Ds
    for (unsigned int p = 0; p < mesh.Cell3DNumberVertices(cell3DIndex); ++p)
    {
        const unsigned int cell0D_index = mesh.Cell3DVertex(cell3DIndex, p);

        const auto& boundary_info = mesh_dofs_info.CellsBoundaryInfo.at(0).at(cell0D_index);

        if (boundary_info.Type !=
            Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::Strong)
            continue;

        const auto coordinates = mesh.Cell0DCoordinates(cell0D_index);

        const auto strong_boundary_values = test.strong_boundary_condition(boundary_info.Marker,
                                                                           coordinates);

        const auto local_dofs = dofs_data.CellsDOFs.at(0).at(cell0D_index);

        assert(local_dofs.size() == strong_boundary_values.size());

        for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
        {
            const auto& local_dof_i = local_dofs.at(loc_i);

            switch (local_dof_i.Type)
            {
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
            {
                assembler_data.solutionDirichlet.SetValue(local_dof_i.Global_Index,
                                                          strong_boundary_values[loc_i]);
            }
            break;
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                continue;
            default:
                throw std::runtime_error("Unknown DOF Type");
            }
        }
    }

    // Assemble strong boundary condition on Cell1Ds
    const auto& referenceSegmentInternalPoints = reference_element_data_2D.Quadrature.ReferenceSegmentInternalPoints;
    const unsigned int numReferenceSegmentInternalPoints = referenceSegmentInternalPoints.cols();

    for (unsigned int e = 0; e < mesh.Cell3DNumberEdges(cell3DIndex); ++e)
    {

        const unsigned int cell1D_index = mesh.Cell3DEdge(cell3DIndex, e);
        const auto local_dofs = dofs_data.CellsDOFs.at(1).at(cell1D_index);
        const auto& boundary_info = mesh_dofs_info.CellsBoundaryInfo.at(1).at(cell1D_index);

        if (boundary_info.Type !=
            Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::Strong
            || local_dofs.size() == 0)
            continue;

        const auto cell1D_origin = mesh.Cell1DOriginCoordinates(cell1D_index);
        const auto cell1D_end = mesh.Cell1DEndCoordinates(cell1D_index);
        const auto cell1D_tangent = cell1D_end - cell1D_origin;

        Eigen::MatrixXd coordinates = Eigen::MatrixXd::Zero(3, numReferenceSegmentInternalPoints);
        for (unsigned int r = 0; r < numReferenceSegmentInternalPoints; r++)
            coordinates.col(r)<< cell1D_origin +
                                      referenceSegmentInternalPoints(0, r) * cell1D_tangent;

        const auto strong_boundary_values = test.strong_boundary_condition(boundary_info.Marker,
                                                                           coordinates);

        assert(local_dofs.size() == strong_boundary_values.size());

        for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
        {
            const auto& local_dof_i = local_dofs.at(loc_i);

            switch (local_dof_i.Type)
            {
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
            {
                assembler_data.solutionDirichlet.SetValue(local_dof_i.Global_Index,
                                                          strong_boundary_values[loc_i]);
            }
            break;
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                continue;
            default:
                throw std::runtime_error("Unknown DOF Type");
            }
        }
    }

    // Assemble strong boundary condition on Cell2Ds
    unsigned int quadraturePointOffset = 0;
    for (unsigned int f = 0; f < mesh.Cell3DNumberFaces(cell3DIndex); f++)
    {
        const unsigned int cell2D_index = mesh.Cell3DFace(cell3DIndex,
                                                          f);

        const auto local_dofs = dofs_data.CellsDOFs.at(2).at(cell2D_index);
        const auto& boundary_info = mesh_dofs_info.CellsBoundaryInfo.at(2).at(cell2D_index);

        const unsigned int numFaceQuadraturePoints = local_space_data.facesLocalSpace[f].InternalQuadrature.Points.cols();

        if (boundary_info.Type !=
            Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::Strong || local_dofs.size() == 0)
        {
            quadraturePointOffset += numFaceQuadraturePoints;
            continue;
        }

        const Eigen::MatrixXd facesInternalQuadraturePoints3D = local_space_data.BoundaryQuadrature.Quadrature.Points.block(0,
                                                                                                                            quadraturePointOffset,
                                                                                                                            3,
                                                                                                                            numFaceQuadraturePoints);

        const Eigen::VectorXd dirichletValues = test.strong_boundary_condition(boundary_info.Marker,
                                                                               facesInternalQuadraturePoints3D);

        const Eigen::VectorXd strong_boundary_values = local_space_data.FaceScaledMomentsBasis[f].transpose()
                                                       * local_space_data.facesLocalSpace[f].InternalQuadrature.Weights.asDiagonal()
                                                       * dirichletValues;

        assert(local_dofs.size() == strong_boundary_values.size());

        for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
        {
            const auto& local_dof_i = local_dofs.at(loc_i);

            switch (local_dof_i.Type)
            {
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
            {
                assembler_data.solutionDirichlet.SetValue(local_dof_i.Global_Index,
                                                          strong_boundary_values[loc_i]);
            }
            break;
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                continue;
            default:
                throw std::runtime_error("Unknown DOF Type");
            }
        }

        quadraturePointOffset += numFaceQuadraturePoints;
    }
}
// ***************************************************************************
void Assembler::ComputeWeakTerm(const unsigned int cell2DIndex,
                                const Gedim::MeshMatricesDAO& mesh,
                                const Polydim::VEM::PCC::VEM_PCC_3D_Polyhedron_Geometry& polygon,
                                const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo& mesh_dofs_info,
                                const Polydim::PDETools::DOFs::DOFsManager::DOFsData& dofs_data,
                                const Polydim::VEM::PCC::VEM_PCC_3D_ReferenceElement_Data& reference_element_data,
                                const Polydim::VEM::PCC::I_VEM_PCC_3D_LocalSpace& vem_local_space,
                                const Polydim::examples::Elliptic_PCC_3D::test::I_Test& test,
                                Elliptic_PCC_3D_Problem_Data& assembler_data) const
{
    //    const unsigned numVertices = polygon.Vertices.cols();

    //    for(unsigned int ed = 0; ed < numVertices; ed ++)
    //    {
    //        const unsigned int cell1D_index = mesh.Cell2DEdge(cell2DIndex,
    //                                                          ed);

    //        const auto& boundary_info = mesh_dofs_info.CellsBoundaryInfo.at(1).at(cell1D_index);

    //        if (boundary_info.Type !=
    //            Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::Weak)
    //            continue;

    //        // compute vem values
    //        const auto weakReferenceSegment = Gedim::Quadrature::Quadrature_Gauss1D::FillPointsAndWeights(2 * reference_element_data.Order);

    //        const Eigen::VectorXd pointsCurvilinearCoordinates = weakReferenceSegment.Points.row(0);

    //        const auto weak_basis_function_values = vem_local_space.ComputeValuesOnEdge(reference_element_data,
    //                                                                                    pointsCurvilinearCoordinates);

    //        // map edge internal quadrature points
    //        const Eigen::Vector3d& edgeStart = polygon.EdgesDirection[ed] ?
    //                                               polygon.Vertices.col(ed) :
    //                                               polygon.Vertices.col((ed + 1) % numVertices);

    //        const Eigen::Vector3d& edgeTangent = polygon.EdgesTangent.col(ed);
    //        const double direction = polygon.EdgesDirection[ed] ? 1.0 : -1.0;

    //        const unsigned int numEdgeWeakQuadraturePoints = weakReferenceSegment.Points.cols();
    //        Eigen::MatrixXd weakQuadraturePoints(3, numEdgeWeakQuadraturePoints);
    //        for (unsigned int q = 0; q < numEdgeWeakQuadraturePoints; q++)
    //        {
    //            weakQuadraturePoints.col(q) = edgeStart + direction *
    //                                                          weakReferenceSegment.Points(0, q) *
    //                                                          edgeTangent;
    //        }
    //        const double absMapDeterminant = std::abs(polygon.EdgesLength[ed]);
    //        const Eigen::MatrixXd weakQuadratureWeights = weakReferenceSegment.Weights *
    //                                                      absMapDeterminant;

    //        const Eigen::VectorXd neumannValues = test.weak_boundary_condition(boundary_info.Marker,
    //                                                                           weakQuadraturePoints);

    //        // compute values of Neumann condition
    //        const Eigen::VectorXd neumannContributions = weak_basis_function_values.transpose() *
    //                                                     weakQuadratureWeights.asDiagonal() *
    //                                                     neumannValues;

    //        for (unsigned int p = 0; p < 2; ++p)
    //        {
    //            const unsigned int cell0D_index = mesh.Cell1DVertex(cell1D_index,
    //                                                                p);

    //            const auto local_dofs = dofs_data.CellsDOFs.at(0).at(cell0D_index);

    //            for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
    //            {
    //                const auto& local_dof_i = local_dofs.at(loc_i);

    //                switch (local_dof_i.Type)
    //                {
    //                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
    //                    continue;
    //                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
    //                {
    //                    assembler_data.rightHandSide.AddValue(local_dof_i.Global_Index,
    //                                                          neumannContributions[p]);
    //                }
    //                break;
    //                default:
    //                    throw std::runtime_error("Unknown DOF Type");
    //                }
    //            }
    //        }

    //        const auto local_dofs = dofs_data.CellsDOFs.at(1).at(cell1D_index);
    //        for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
    //        {
    //            const auto& local_dof_i = local_dofs.at(loc_i);

    //            const unsigned int localIndex = polygon.EdgesDirection[ed] ? loc_i :
    //                                                local_dofs.size() - 1 - loc_i;


    //            switch (local_dof_i.Type)
    //            {
    //            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
    //                continue;
    //            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
    //            {
    //                assembler_data.rightHandSide.AddValue(local_dof_i.Global_Index,
    //                                                      neumannContributions[localIndex + 2]);
    //            }
    //            break;
    //            default:
    //                throw std::runtime_error("Unknown DOF Type");
    //            }
    //        }
    //    }
}
// ***************************************************************************
typename Assembler::Elliptic_PCC_3D_Problem_Data Assembler::Assemble(const Polydim::examples::Elliptic_PCC_3D::Program_configuration& config,
                                                                     const Gedim::MeshMatricesDAO& mesh,
                                                                     const Gedim::MeshUtilities::MeshGeometricData3D& mesh_geometric_data,
                                                                     const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo& mesh_dofs_info,
                                                                     const Polydim::PDETools::DOFs::DOFsManager::DOFsData& dofs_data,
                                                                     const Polydim::VEM::PCC::VEM_PCC_2D_ReferenceElement_Data& reference_element_data_2D,
                                                                     const Polydim::VEM::PCC::VEM_PCC_3D_ReferenceElement_Data& reference_element_data_3D,
                                                                     const Polydim::examples::Elliptic_PCC_3D::test::I_Test& test) const
{
    Elliptic_PCC_3D_Problem_Data result;

    result.globalMatrixA.SetSize(dofs_data.NumberDOFs,
                                 dofs_data.NumberDOFs,
                                 Gedim::ISparseArray::SparseArrayTypes::Symmetric);
    result.dirichletMatrixA.SetSize(dofs_data.NumberDOFs,
                                    dofs_data.NumberStrongs);
    result.rightHandSide.SetSize(dofs_data.NumberDOFs);
    result.solution.SetSize(dofs_data.NumberDOFs);
    result.solutionDirichlet.SetSize(dofs_data.NumberStrongs);

    Polydim::PDETools::Equations::EllipticEquation equation;

    for (unsigned int c = 0; c < mesh.Cell3DTotalNumber(); ++c)
    {
        const unsigned int numFaces = mesh_geometric_data.Cell3DsFaces.at(c).size();

        std::vector<Polydim::VEM::PCC::VEM_PCC_2D_Polygon_Geometry> polygonalFaces;
        for (unsigned int f = 0; f < numFaces; f++)
        {
            polygonalFaces.push_back(
                {
                    config.GeometricTolerance1D(),
                    config.GeometricTolerance2D(),
                    mesh_geometric_data.Cell3DsFaces2DVertices.at(c)[f],
                    mesh_geometric_data.Cell3DsFaces2DCentroids.at(c)[f],
                    mesh_geometric_data.Cell3DsFacesAreas.at(c)[f],
                    mesh_geometric_data.Cell3DsFacesDiameters.at(c)[f],
                    mesh_geometric_data.Cell3DsFaces2DTriangulations.at(c)[f],
                    mesh_geometric_data.Cell3DsFacesEdgeLengths.at(c)[f],
                    mesh_geometric_data.Cell3DsFacesEdgeDirections.at(c)[f],
                    mesh_geometric_data.Cell3DsFacesEdge2DTangents.at(c)[f],
                    mesh_geometric_data.Cell3DsFacesEdge2DNormals.at(c)[f]
                });
        }

        const Polydim::VEM::PCC::VEM_PCC_3D_Polyhedron_Geometry polyhedron =
            {
                config.GeometricTolerance1D(),
                config.GeometricTolerance2D(),
                config.GeometricTolerance3D(),
                mesh_geometric_data.Cell3DsVertices.at(c),
                mesh_geometric_data.Cell3DsEdges.at(c),
                mesh_geometric_data.Cell3DsFaces.at(c),
                mesh_geometric_data.Cell3DsCentroids.at(c),
                mesh_geometric_data.Cell3DsVolumes.at(c),
                mesh_geometric_data.Cell3DsDiameters.at(c),
                mesh_geometric_data.Cell3DsTetrahedronPoints.at(c),
                mesh_geometric_data.Cell3DsFacesRotationMatrices.at(c),
                mesh_geometric_data.Cell3DsFacesTranslations.at(c),
                mesh_geometric_data.Cell3DsFacesNormals.at(c),
                mesh_geometric_data.Cell3DsFacesNormalDirections.at(c),
                mesh_geometric_data.Cell3DsEdgeDirections.at(c),
                mesh_geometric_data.Cell3DsEdgeTangents.at(c)
            };

        const auto vem_local_space = Polydim::VEM::PCC::create_VEM_PCC_3D_local_space_3D(config.VemType());
        const auto local_space = vem_local_space->CreateLocalSpace(reference_element_data_2D,
                                                                   reference_element_data_3D,
                                                                   polygonalFaces,
                                                                   polyhedron);

        const auto basis_functions_values = vem_local_space->ComputeBasisFunctionsValues(local_space,
                                                                                         Polydim::VEM::PCC::ProjectionTypes::Pi0km1);

        const auto basis_functions_derivative_values = vem_local_space->ComputeBasisFunctionsDerivativeValues(local_space,
                                                                                                              Polydim::VEM::PCC::ProjectionTypes::Pi0km1Der);

        const auto diffusion_term_values = test.diffusion_term(local_space.InternalQuadrature.Points);
        const auto source_term_values = test.source_term(local_space.InternalQuadrature.Points);

        const auto local_A = equation.ComputeCellDiffusionMatrix(diffusion_term_values,
                                                                 basis_functions_derivative_values,
                                                                 local_space.InternalQuadrature.Weights);

        const auto local_stab_A = diffusion_term_values.cwiseAbs().maxCoeff() *
                                  local_space.StabMatrix;

        const auto local_rhs = equation.ComputeCellForcingTerm(source_term_values,
                                                               basis_functions_values,
                                                               local_space.InternalQuadrature.Weights);

        const auto& global_dofs = dofs_data.CellsGlobalDOFs[3].at(c);

        assert(local_space.NumBasisFunctions ==  global_dofs.size());


        Polydim::PDETools::Assembler_Utilities::local_matrix_to_global_matrix_dofs_data
            local_matrix_to_global_matrix_dofs_data =
            {
                { std::cref(dofs_data) },
                { 0 },
                { 0 },
                { 0 }
            };

        Polydim::PDETools::Assembler_Utilities::assemble_local_matrix_to_global_matrix<3>(c,
                                                                                          local_matrix_to_global_matrix_dofs_data,
                                                                                          local_matrix_to_global_matrix_dofs_data,
                                                                                          local_A + local_stab_A,
                                                                                          local_rhs,
                                                                                          result.globalMatrixA,
                                                                                          result.dirichletMatrixA,
                                                                                          result.rightHandSide);

        //        ComputeWeakTerm(c,
        //                        mesh,
        //                        polygon,
        //                        mesh_dofs_info,
        //                        dofs_data,
        //                        reference_element_data,
        //                        *vem_local_space,
        //                        test,
        //                        result);

        ComputeStrongTerm(c,
                          mesh,
                          polygonalFaces,
                          mesh_dofs_info,
                          dofs_data,
                          reference_element_data_2D,
                          reference_element_data_3D,
                          local_space,
                          test,
                          result);
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
Assembler::VEM_Performance_Result Assembler::ComputeVemPerformance(const Polydim::examples::Elliptic_PCC_3D::Program_configuration& config,
                                                                   const Gedim::MeshMatricesDAO& mesh,
                                                                   const Gedim::MeshUtilities::MeshGeometricData3D& mesh_geometric_data,
                                                                   const Polydim::VEM::PCC::VEM_PCC_2D_ReferenceElement_Data& reference_element_data_2D,
                                                                   const Polydim::VEM::PCC::VEM_PCC_3D_ReferenceElement_Data& reference_element_data_3D) const
{

    Assembler::VEM_Performance_Result result;
    result.Cell3DsPerformance.resize(mesh.Cell3DTotalNumber());

    // Assemble equation elements
    for (unsigned int c = 0; c < mesh.Cell3DTotalNumber(); c++)
    {
        const unsigned int numFaces = mesh_geometric_data.Cell3DsFaces.at(c).size();

        std::vector<Polydim::VEM::PCC::VEM_PCC_2D_Polygon_Geometry> polygonalFaces;
        for (unsigned int f = 0; f < numFaces; f++)
        {
            polygonalFaces.push_back(
                {
                    config.GeometricTolerance1D(),
                    config.GeometricTolerance2D(),
                    mesh_geometric_data.Cell3DsFaces2DVertices.at(c)[f],
                    mesh_geometric_data.Cell3DsFaces2DCentroids.at(c)[f],
                    mesh_geometric_data.Cell3DsFacesAreas.at(c)[f],
                    mesh_geometric_data.Cell3DsFacesDiameters.at(c)[f],
                    mesh_geometric_data.Cell3DsFaces2DTriangulations.at(c)[f],
                    mesh_geometric_data.Cell3DsFacesEdgeLengths.at(c)[f],
                    mesh_geometric_data.Cell3DsFacesEdgeDirections.at(c)[f],
                    mesh_geometric_data.Cell3DsFacesEdge2DTangents.at(c)[f],
                    mesh_geometric_data.Cell3DsFacesEdge2DNormals.at(c)[f]
                });
        }

        const Polydim::VEM::PCC::VEM_PCC_3D_Polyhedron_Geometry polyhedron =
            {
                config.GeometricTolerance1D(),
                config.GeometricTolerance2D(),
                config.GeometricTolerance3D(),
                mesh_geometric_data.Cell3DsVertices.at(c),
                mesh_geometric_data.Cell3DsEdges.at(c),
                mesh_geometric_data.Cell3DsFaces.at(c),
                mesh_geometric_data.Cell3DsCentroids.at(c),
                mesh_geometric_data.Cell3DsVolumes.at(c),
                mesh_geometric_data.Cell3DsDiameters.at(c),
                mesh_geometric_data.Cell3DsTetrahedronPoints.at(c),
                mesh_geometric_data.Cell3DsFacesRotationMatrices.at(c),
                mesh_geometric_data.Cell3DsFacesTranslations.at(c),
                mesh_geometric_data.Cell3DsFacesNormals.at(c),
                mesh_geometric_data.Cell3DsFacesNormalDirections.at(c),
                mesh_geometric_data.Cell3DsEdgeDirections.at(c),
                mesh_geometric_data.Cell3DsEdgeTangents.at(c)
            };

        const auto vem_local_space = Polydim::VEM::PCC::create_VEM_PCC_3D_local_space_3D(config.VemType());

        const auto local_space = vem_local_space->CreateLocalSpace(reference_element_data_2D,
                                                                   reference_element_data_3D,
                                                                   polygonalFaces,
                                                                   polyhedron);

        Polydim::VEM::PCC::VEM_PCC_PerformanceAnalysis performanceAnalysis;

        result.Cell3DsPerformance[c].Analysis = performanceAnalysis.Compute(Polydim::VEM::Monomials::VEM_Monomials_3D(),
                                                                            reference_element_data_3D.Monomials,
                                                                            *vem_local_space,
                                                                            local_space);

        result.Cell3DsPerformance[c].NumInternalQuadraturePoints = local_space.InternalQuadrature.Weights.size();
        result.Cell3DsPerformance[c].NumBoundaryQuadraturePoints = local_space.BoundaryQuadrature.Quadrature.Weights.size();
    }

    return result;
}
// ***************************************************************************
Assembler::PostProcess_Data Assembler::PostProcessSolution(const Polydim::examples::Elliptic_PCC_3D::Program_configuration& config,
                                                           const Gedim::MeshMatricesDAO& mesh,
                                                           const Gedim::MeshUtilities::MeshGeometricData3D& mesh_geometric_data,
                                                           const Polydim::PDETools::DOFs::DOFsManager::DOFsData& dofs_data,
                                                           const Polydim::VEM::PCC::VEM_PCC_2D_ReferenceElement_Data& reference_element_data_2D,
                                                           const Polydim::VEM::PCC::VEM_PCC_3D_ReferenceElement_Data& reference_element_data_3D,
                                                           const Elliptic_PCC_3D_Problem_Data& assembler_data,
                                                           const Polydim::examples::Elliptic_PCC_3D::test::I_Test& test) const
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
        result.cell0Ds_exact[p] = test.exact_solution(mesh.Cell0DCoordinates(p))[0];

        const auto local_dofs = dofs_data.CellsDOFs.at(0).at(p);

        for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
        {
            const auto& local_dof_i = local_dofs.at(loc_i);

            switch (local_dof_i.Type)
            {
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
                result.cell0Ds_numeric[p] = assembler_data.solutionDirichlet.GetValue(local_dof_i.Global_Index);
                break;
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                result.cell0Ds_numeric[p] = assembler_data.solution.GetValue(local_dof_i.Global_Index);
                break;
            default:
                throw std::runtime_error("Unknown DOF Type");
            }
        }
    }

    result.cell3Ds_error_L2.setZero(mesh.Cell3DTotalNumber());
    result.cell3Ds_norm_L2.setZero(mesh.Cell3DTotalNumber());
    result.cell3Ds_error_H1.setZero(mesh.Cell3DTotalNumber());
    result.cell3Ds_norm_H1.setZero(mesh.Cell3DTotalNumber());
    result.error_L2 = 0.0;
    result.norm_L2 = 0.0;
    result.error_H1 = 0.0;
    result.norm_H1 = 0.0;
    result.mesh_size = 0.0;

    for (unsigned int c = 0; c < mesh.Cell3DTotalNumber(); c++)
    {
        const unsigned int numFaces = mesh_geometric_data.Cell3DsFaces.at(c).size();

        std::vector<Polydim::VEM::PCC::VEM_PCC_2D_Polygon_Geometry> polygonalFaces;
        for (unsigned int f = 0; f < numFaces; f++)
        {
            polygonalFaces.push_back(
                {
                    config.GeometricTolerance1D(),
                    config.GeometricTolerance2D(),
                    mesh_geometric_data.Cell3DsFaces2DVertices.at(c)[f],
                    mesh_geometric_data.Cell3DsFaces2DCentroids.at(c)[f],
                    mesh_geometric_data.Cell3DsFacesAreas.at(c)[f],
                    mesh_geometric_data.Cell3DsFacesDiameters.at(c)[f],
                    mesh_geometric_data.Cell3DsFaces2DTriangulations.at(c)[f],
                    mesh_geometric_data.Cell3DsFacesEdgeLengths.at(c)[f],
                    mesh_geometric_data.Cell3DsFacesEdgeDirections.at(c)[f],
                    mesh_geometric_data.Cell3DsFacesEdge2DTangents.at(c)[f],
                    mesh_geometric_data.Cell3DsFacesEdge2DNormals.at(c)[f]
                });
        }

        const Polydim::VEM::PCC::VEM_PCC_3D_Polyhedron_Geometry polyhedron =
            {
                config.GeometricTolerance1D(),
                config.GeometricTolerance2D(),
                config.GeometricTolerance3D(),
                mesh_geometric_data.Cell3DsVertices.at(c),
                mesh_geometric_data.Cell3DsEdges.at(c),
                mesh_geometric_data.Cell3DsFaces.at(c),
                mesh_geometric_data.Cell3DsCentroids.at(c),
                mesh_geometric_data.Cell3DsVolumes.at(c),
                mesh_geometric_data.Cell3DsDiameters.at(c),
                mesh_geometric_data.Cell3DsTetrahedronPoints.at(c),
                mesh_geometric_data.Cell3DsFacesRotationMatrices.at(c),
                mesh_geometric_data.Cell3DsFacesTranslations.at(c),
                mesh_geometric_data.Cell3DsFacesNormals.at(c),
                mesh_geometric_data.Cell3DsFacesNormalDirections.at(c),
                mesh_geometric_data.Cell3DsEdgeDirections.at(c),
                mesh_geometric_data.Cell3DsEdgeTangents.at(c)
            };

        const auto vem_local_space = Polydim::VEM::PCC::create_VEM_PCC_3D_local_space_3D(config.VemType());

        const auto local_space = vem_local_space->CreateLocalSpace(reference_element_data_2D,
                                                                   reference_element_data_3D,
                                                                   polygonalFaces,
                                                                   polyhedron);

        const auto basis_functions_values = vem_local_space->ComputeBasisFunctionsValues(local_space,
                                                                                         Polydim::VEM::PCC::ProjectionTypes::Pi0k);

        const auto basis_functions_derivative_values = vem_local_space->ComputeBasisFunctionsDerivativeValues(local_space,
                                                                                                              Polydim::VEM::PCC::ProjectionTypes::Pi0km1Der);


        const auto exact_solution_values = test.exact_solution(local_space.InternalQuadrature.Points);
        const auto exact_derivative_solution_values = test.exact_derivative_solution(local_space.InternalQuadrature.Points);

        const auto& global_dofs = dofs_data.CellsGlobalDOFs[3].at(c);
        Eigen::VectorXd dofs_values = Eigen::VectorXd::Zero(global_dofs.size());

        for (unsigned int loc_i = 0; loc_i < global_dofs.size(); ++loc_i)
        {
            const auto& global_dof_i = global_dofs.at(loc_i);
            const auto& local_dof_i = dofs_data.CellsDOFs.at(global_dof_i.Dimension).at(global_dof_i.CellIndex).at(global_dof_i.DOFIndex);

            switch (local_dof_i.Type)
            {
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
                dofs_values[loc_i] = assembler_data.solutionDirichlet.GetValue(local_dof_i.Global_Index);
                break;
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                dofs_values[loc_i] = assembler_data.solution.GetValue(local_dof_i.Global_Index);
                break;
            default:
                throw std::runtime_error("Unknown DOF Type");
            }
        }

        const Eigen::VectorXd local_error_L2 = (basis_functions_values * dofs_values -
                                                exact_solution_values).array().square();
        const Eigen::VectorXd local_norm_L2 = (basis_functions_values * dofs_values).array().square();

        result.cell3Ds_error_L2[c] = local_space.InternalQuadrature.Weights.transpose() *
                                     local_error_L2;
        result.cell3Ds_norm_L2[c] = local_space.InternalQuadrature.Weights.transpose() *
                                    local_norm_L2;


        const Eigen::VectorXd local_error_H1 =
            (basis_functions_derivative_values[0] * dofs_values -
             exact_derivative_solution_values[0]).array().square() +
            (basis_functions_derivative_values[1] * dofs_values -
             exact_derivative_solution_values[1]).array().square();
        const Eigen::VectorXd local_norm_H1 =
            (basis_functions_derivative_values[0] * dofs_values).array().square() +
            (basis_functions_derivative_values[1] * dofs_values).array().square();

        result.cell3Ds_error_H1[c] = local_space.InternalQuadrature.Weights.transpose() *
                                     local_error_H1;
        result.cell3Ds_norm_H1[c] = local_space.InternalQuadrature.Weights.transpose() *
                                    local_norm_H1;

        if (mesh_geometric_data.Cell3DsDiameters.at(c) > result.mesh_size)
            result.mesh_size = mesh_geometric_data.Cell3DsDiameters.at(c);
    }

    result.error_L2 = std::sqrt(result.cell3Ds_error_L2.sum());
    result.norm_L2 = std::sqrt(result.cell3Ds_norm_L2.sum());
    result.error_H1 = std::sqrt(result.cell3Ds_error_H1.sum());
    result.norm_H1 = std::sqrt(result.cell3Ds_norm_H1.sum());

    return result;
}
}
}
}
