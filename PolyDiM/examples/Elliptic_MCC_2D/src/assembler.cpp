#include "assembler.hpp"

#include "ranges"

#include "VEM_MCC_2D_VelocityLocalSpace.hpp"
#include "EllipticEquation.hpp"

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
                                                                                                const Eigen::MatrixXd&)>& strong_boundary_condition,
                                                            const std::function<Eigen::VectorXd(const unsigned int,
                                                                                                const Eigen::MatrixXd&)>& weak_boundary_condition) const
{

    const unsigned int numDOFHandler = mesh_dofs_info.size();
    unsigned int numberDOFs = 0;
    unsigned int numberStrongs = 0;
    std::vector<unsigned int> offsetDOFs = {0, dofs_data[0].NumberDOFs};
    std::vector<unsigned int> offsetStrongs = {0, dofs_data[0].NumberStrongs};
    for(unsigned int i = 0; i < numDOFHandler; i++)
    {
        numberDOFs += dofs_data[i].NumberDOFs;
        numberStrongs += dofs_data[i].NumberStrongs;
    }

    Elliptic_MCC_2D_Problem_Data result;
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

        ComputeWeakTerm(c,
                        mesh,
                        polygon,
                        mesh_dofs_info[0],
                        dofs_data[0],
                        velocity_reference_element_data,
                        local_space,
                        weak_boundary_condition,
                        result);
    }

    //    ComputeStrongTerm(geometryUtilities,
    //                      mesh,
    //                      mesh_geometric_data,
    //                      mesh_dofs_info[0],
    //                      dofs_data[0],
    //                      velocity_reference_element_data,
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
void Assembler::ComputeStrongTerm(const Gedim::GeometryUtilities& geometryUtilities,
                                  const Gedim::MeshMatricesDAO& mesh,
                                  const Gedim::MeshUtilities::MeshGeometricData2D& mesh_geometric_data,
                                  const Polydim::PDETools::DOFs::DOFsManager<2>::MeshDOFsInfo& mesh_dofs_info,
                                  const Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData& dofs_data,
                                  const Polydim::VEM::MCC::VEM_MCC_2D_Velocity_ReferenceElement_Data& reference_element_data,
                                  const std::function<Eigen::VectorXd(const unsigned int,
                                                                      const Eigen::MatrixXd&)>& strong_boundary_condition,
                                  Elliptic_MCC_2D_Problem_Data& assembler_data) const
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
                assembler_data.solutionNeumann.SetValue(local_dof_i.Global_Index,
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
                assembler_data.solutionNeumann.SetValue(local_dof_i.Global_Index,
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
                                const Polydim::VEM::MCC::VEM_MCC_2D_Polygon_Geometry& polygon,
                                const Polydim::PDETools::DOFs::DOFsManager<2>::MeshDOFsInfo& mesh_dofs_info,
                                const Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData& dofs_data,
                                const Polydim::VEM::MCC::VEM_MCC_2D_Velocity_ReferenceElement_Data& reference_element_data,
                                const Polydim::VEM::MCC::VEM_MCC_VelocityLocalSpace_Data& local_space_data,
                                const std::function<Eigen::VectorXd(const unsigned int,
                                                                    const Eigen::MatrixXd&)>& weak_boundary_condition,
                                Elliptic_MCC_2D_Problem_Data& assembler_data) const
{

    const unsigned int numQuadraturePoints = reference_element_data.Quadrature.ReferenceSegmentInternalPoints.cols();
    for(unsigned int e = 0; e < mesh.Cell2DNumberEdges(cell2DIndex); e ++)
    {
        const unsigned int edgeGlobalId = mesh.Cell2DEdge(cell2DIndex,e);

        const unsigned int cell1D_index = mesh.Cell2DEdge(cell2DIndex,
                                                          e);

        const auto& boundary_info = mesh_dofs_info.CellsBoundaryInfo.at(1).at(cell1D_index);

        if (boundary_info.Type !=
            Polydim::PDETools::DOFs::DOFsManager<2>::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::Weak)
            continue;

        // compute values of Neumann condition
        const double direction = polygon.EdgesDirection[e] ? 1.0 : -1.0;

        const VectorXd dirichletValues = weak_boundary_condition(boundary_info.Marker,
                                                                 local_space_data.BoundaryQuadrature.Quadrature.Points.middleCols(numQuadraturePoints * e,
                                                                                                                                  numQuadraturePoints));

        const VectorXd dirichletContributions = - direction *
                                                dirichletValues.cwiseProduct(local_space_data.BoundaryQuadrature.Quadrature.Weights.segment(numQuadraturePoints * e,
                                                                                                                                            numQuadraturePoints));

        const auto local_dofs = dofs_data.CellsDOFs.at(1).at(cell1D_index);

        for (unsigned int loc_i = 0; loc_i < local_dofs.size(); loc_i++)
        {
            const auto& local_dof_i = local_dofs.at(loc_i);

            switch (local_dof_i.Type)
            {
            case Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData::DOF::Types::Strong:
                continue;
            case Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData::DOF::Types::DOF:
            {
                assembler_data.rightHandSide.AddValue(local_dof_i.Global_Index,
                                                      dirichletContributions(loc_i));
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
Assembler::PostProcess_Data Assembler::PostProcessSolution(const Gedim::GeometryUtilities& geometryUtilities,
                                                           const Gedim::MeshMatricesDAO& mesh,
                                                           const Gedim::MeshUtilities::MeshGeometricData2D& mesh_geometric_data,
                                                           const vector<Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData>& dofs_data,
                                                           const Polydim::VEM::MCC::VEM_MCC_2D_Velocity_ReferenceElement_Data& velocity_reference_element_data,
                                                           const Polydim::VEM::MCC::VEM_MCC_2D_Pressure_ReferenceElement_Data& pressure_reference_element_data,
                                                           const Elliptic_MCC_2D_Problem_Data& assembler_data,
                                                           const std::function<std::array<Eigen::VectorXd, 3>(const Eigen::MatrixXd&)>& exact_velocity,
                                                           const std::function<Eigen::VectorXd(const Eigen::MatrixXd&)>& exact_pressure) const
{
    const unsigned int numDOFHandler = dofs_data.size();
    unsigned int numberDOFs = 0;
    unsigned int numberStrongs = 0;
    std::vector<unsigned int> offsetDOFs = {0, dofs_data[0].NumberDOFs};
    std::vector<unsigned int> offsetStrongs = {0, dofs_data[0].NumberStrongs};
    for(unsigned int i = 0; i < numDOFHandler; i++)
    {
        numberDOFs += dofs_data[i].NumberDOFs;
        numberStrongs += dofs_data[i].NumberStrongs;
    }

    PostProcess_Data result;

    result.residual_norm = 0.0;
    if (numberDOFs > 0)
    {
        Gedim::Eigen_Array<> residual;
        residual.SetSize(numberDOFs);
        residual.SumMultiplication(assembler_data.globalMatrixA,
                                   assembler_data.solution);
        residual -= assembler_data.rightHandSide;

        result.residual_norm = residual.Norm();
    }

    result.cell0Ds_numeric_pressure.resize(mesh.Cell2DTotalNumber());
    result.cell0Ds_exact_pressure.resize(mesh.Cell2DTotalNumber());

    result.cell2Ds_error_L2_pressure.setZero(mesh.Cell2DTotalNumber());
    result.cell2Ds_norm_L2_pressure.setZero(mesh.Cell2DTotalNumber());
    result.cell2Ds_error_L2_velocity.setZero(mesh.Cell2DTotalNumber());
    result.cell2Ds_norm_L2_velocity.setZero(mesh.Cell2DTotalNumber());
    result.error_L2_pressure = 0.0;
    result.norm_L2_pressure = 0.0;
    result.error_L2_velocity = 0.0;
    result.norm_L2_velocity = 0.0;
    result.mesh_size = 0.0;

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
        const auto pressure_basis_functions_values = monomials.Vander(pressure_reference_element_data.Monomials,
                                                                      local_space.InternalQuadrature.Points,
                                                                      mesh_geometric_data.Cell2DsCentroids.at(c),
                                                                      mesh_geometric_data.Cell2DsDiameters.at(c));

        const auto& global_dofs_velocity = dofs_data[0].CellsGlobalDOFs[2].at(c);
        Eigen::VectorXd velocity_dofs_values = Eigen::VectorXd::Zero(global_dofs_velocity.size());
        for (unsigned int loc_i = 0; loc_i < global_dofs_velocity.size(); ++loc_i)
        {
            const auto& global_dof_i = global_dofs_velocity.at(loc_i);
            const auto& local_dof_i = dofs_data[0].CellsDOFs.at(global_dof_i.Dimension).at(global_dof_i.CellIndex).at(global_dof_i.DOFIndex);

            switch (local_dof_i.Type)
            {
            case Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData::DOF::Types::Strong:
                velocity_dofs_values[loc_i] = assembler_data.solutionNeumann.GetValue(local_dof_i.Global_Index);
                break;
            case Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData::DOF::Types::DOF:
                velocity_dofs_values[loc_i] = assembler_data.solution.GetValue(local_dof_i.Global_Index);
                break;
            default:
                throw std::runtime_error("Unknown DOF Type");
            }
        }

        const auto& pressure_global_dofs = dofs_data[1].CellsGlobalDOFs[2].at(c);
        Eigen::VectorXd pressure_dofs_values = Eigen::VectorXd::Zero(pressure_global_dofs.size());
        for (unsigned int loc_i = 0; loc_i < pressure_global_dofs.size(); ++loc_i)
        {
            const auto& global_dof_i = pressure_global_dofs.at(loc_i);
            const auto& local_dof_i = dofs_data[1].CellsDOFs.at(global_dof_i.Dimension).at(global_dof_i.CellIndex).at(global_dof_i.DOFIndex);

            switch (local_dof_i.Type)
            {
            case Polydim::PDETools::DOFs::DOFsManager<2>::DOFsData::DOF::Types::DOF:
                pressure_dofs_values[loc_i] = assembler_data.solution.GetValue(local_dof_i.Global_Index + offsetDOFs[1]);
                break;
            default:
                throw std::runtime_error("Unknown DOF Type");
            }
        }


        const auto exact_pressure_values = exact_pressure(local_space.InternalQuadrature.Points);
        const auto exact_velocity_values = exact_velocity(local_space.InternalQuadrature.Points);

        const Eigen::VectorXd local_error_L2_pressure = (pressure_basis_functions_values * pressure_dofs_values -
                                                         exact_pressure_values).array().square();
        const Eigen::VectorXd local_norm_L2_pressure = (pressure_basis_functions_values * pressure_dofs_values).array().square();

        result.cell2Ds_error_L2_pressure[c] = local_space.InternalQuadrature.Weights.transpose() *
                                              local_error_L2_pressure;
        result.cell2Ds_norm_L2_pressure[c] = local_space.InternalQuadrature.Weights.transpose() *
                                             local_norm_L2_pressure;


        const Eigen::VectorXd local_error_L2_velocity =
            (velocity_basis_functions_values[0] * velocity_dofs_values -
             exact_velocity_values[0]).array().square() +
            (velocity_basis_functions_values[1] * velocity_dofs_values -
             exact_velocity_values[1]).array().square();
        const Eigen::VectorXd local_norm_L2_velocity =
            (velocity_basis_functions_values[0] * velocity_dofs_values).array().square() +
            (velocity_basis_functions_values[1] * velocity_dofs_values).array().square();

        result.cell2Ds_error_L2_velocity[c] = local_space.InternalQuadrature.Weights.transpose() *
                                     local_error_L2_velocity;
        result.cell2Ds_norm_L2_velocity[c] = local_space.InternalQuadrature.Weights.transpose() *
                                    local_norm_L2_velocity;

        if (mesh_geometric_data.Cell2DsDiameters.at(c) > result.mesh_size)
            result.mesh_size = mesh_geometric_data.Cell2DsDiameters.at(c);
    }

    result.error_L2_pressure = std::sqrt(result.cell2Ds_error_L2_pressure.sum());
    result.norm_L2_pressure = std::sqrt(result.cell2Ds_norm_L2_pressure.sum());
    result.error_L2_velocity = std::sqrt(result.cell2Ds_error_L2_velocity.sum());
    result.norm_L2_velocity = std::sqrt(result.cell2Ds_norm_L2_velocity.sum());

    return result;
}
// ***************************************************************************
}
