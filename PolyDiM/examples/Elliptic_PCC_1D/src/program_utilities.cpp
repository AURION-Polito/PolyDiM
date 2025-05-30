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

#include "program_utilities.hpp"

namespace Polydim
{
namespace examples
{
namespace Elliptic_PCC_1D
{
namespace program_utilities
{
// ***************************************************************************
std::unique_ptr<Polydim::examples::Elliptic_PCC_1D::test::I_Test> create_test(const Polydim::examples::Elliptic_PCC_1D::Program_configuration &config)
{
    switch (config.TestType())
    {
    case Polydim::examples::Elliptic_PCC_1D::test::Test_Types::Patch_Test:
        return std::make_unique<Polydim::examples::Elliptic_PCC_1D::test::Patch_Test>();
    case Polydim::examples::Elliptic_PCC_1D::test::Test_Types::Poisson_Polynomial_Problem:
        return std::make_unique<Polydim::examples::Elliptic_PCC_1D::test::Poisson_Polynomial_Problem>();
    default:
        throw std::runtime_error("Test type " + std::to_string((unsigned int)config.TestType()) + " not supported");
    }
}
// ***************************************************************************
void create_domain_mesh(const Polydim::examples::Elliptic_PCC_1D::Program_configuration &config,
                        const Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_1D &domain,
                        Gedim::MeshMatricesDAO &mesh)
{
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = config.GeometricTolerance1D();
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::MeshUtilities meshUtilities;

    switch (config.MeshGenerator())
    {
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_1D::Minimal:
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_1D::Equispaced: {
        Polydim::PDETools::Mesh::PDE_Mesh_Utilities::create_mesh_1D(geometryUtilities,
                                                                    meshUtilities,
                                                                    config.MeshGenerator(),
                                                                    domain,
                                                                    config.MeshMaxLength(),
                                                                    mesh);
    }
    break;
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_1D::CsvImporter: {
        Polydim::PDETools::Mesh::PDE_Mesh_Utilities::import_mesh_1D(geometryUtilities,
                                                                    meshUtilities,
                                                                    config.MeshGenerator(),
                                                                    config.MeshImportFilePath(),
                                                                    mesh);
    }
    break;
    default:
        throw std::runtime_error("MeshGenerator " + std::to_string((unsigned int)config.MeshGenerator()) + " not supported");
    }
}
// ***************************************************************************
Gedim::MeshUtilities::MeshGeometricData1D create_domain_mesh_geometric_properties(const Polydim::examples::Elliptic_PCC_1D::Program_configuration &config,
                                                                                  const Gedim::MeshMatricesDAO &mesh)
{
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = config.GeometricTolerance1D();
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::MeshUtilities meshUtilities;

    return Polydim::PDETools::Mesh::PDE_Mesh_Utilities::compute_mesh_1D_geometry_data(geometryUtilities, meshUtilities, mesh);
}
// ***************************************************************************
void export_solution(const Polydim::examples::Elliptic_PCC_1D::Program_configuration &config,
                     const Gedim::MeshMatricesDAO &mesh,
                     const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                     const Polydim::examples::Elliptic_PCC_1D::Assembler::Elliptic_PCC_1D_Problem_Data &assembler_data,
                     const Polydim::examples::Elliptic_PCC_1D::Assembler::PostProcess_Data &post_process_data,
                     const std::string &exportSolutionFolder,
                     const std::string &exportVtuFolder)
{
    const unsigned int METHOD_ID = static_cast<unsigned int>(config.MethodType());
    const unsigned int TEST_ID = static_cast<unsigned int>(config.TestType());

    {
        const char separator = ';';

        std::cout << "ProgramType" << separator;
        std::cout << "MethodType" << separator;
        std::cout << "MethodOrder" << separator;
        std::cout << "Cell1Ds" << separator;
        std::cout << "Dofs" << separator;
        std::cout << "Strongs" << separator;
        std::cout << "h" << separator;
        std::cout << "errorL2" << separator;
        std::cout << "errorH1" << separator;
        std::cout << "normL2" << separator;
        std::cout << "normH1" << separator;
        std::cout << "nnzA" << separator;
        std::cout << "residual" << std::endl;

        std::cout.precision(2);
        std::cout << std::scientific << TEST_ID << separator;
        std::cout << std::scientific << METHOD_ID << separator;
        std::cout << std::scientific << config.MethodOrder() << separator;
        std::cout << std::scientific << mesh.Cell1DTotalNumber() << separator;
        std::cout << std::scientific << dofs_data.NumberDOFs << separator;
        std::cout << std::scientific << dofs_data.NumberStrongs << separator;
        std::cout << std::scientific << post_process_data.mesh_size << separator;
        std::cout << std::scientific << post_process_data.error_L2 << separator;
        std::cout << std::scientific << post_process_data.error_H1 << separator;
        std::cout << std::scientific << post_process_data.norm_L2 << separator;
        std::cout << std::scientific << post_process_data.norm_H1 << separator;
        std::cout << std::scientific << assembler_data.globalMatrixA.NonZeros() << separator;
        std::cout << std::scientific << post_process_data.residual_norm << std::endl;
    }

    {
        const char separator = ';';
        const std::string errorFileName = exportSolutionFolder + "/Errors_" + std::to_string(TEST_ID) + "_" +
                                          std::to_string(METHOD_ID) + +"_" + std::to_string(config.MethodOrder()) + ".csv";
        const bool errorFileExists = Gedim::Output::FileExists(errorFileName);

        std::ofstream errorFile(errorFileName, std::ios_base::app | std::ios_base::out);
        if (!errorFileExists)
        {
            errorFile << "ProgramType" << separator;
            errorFile << "MethodType" << separator;
            errorFile << "MethodOrder" << separator;
            errorFile << "Cell1Ds" << separator;
            errorFile << "Dofs" << separator;
            errorFile << "Strongs" << separator;
            errorFile << "h" << separator;
            errorFile << "errorL2" << separator;
            errorFile << "errorH1" << separator;
            errorFile << "normL2" << separator;
            errorFile << "normH1" << separator;
            errorFile << "nnzA" << separator;
            errorFile << "residual" << std::endl;
        }

        errorFile.precision(16);
        errorFile << std::scientific << TEST_ID << separator;
        errorFile << std::scientific << METHOD_ID << separator;
        errorFile << std::scientific << config.MethodOrder() << separator;
        errorFile << std::scientific << mesh.Cell1DTotalNumber() << separator;
        errorFile << std::scientific << dofs_data.NumberDOFs << separator;
        errorFile << std::scientific << dofs_data.NumberStrongs << separator;
        errorFile << std::scientific << post_process_data.mesh_size << separator;
        errorFile << std::scientific << post_process_data.error_L2 << separator;
        errorFile << std::scientific << post_process_data.error_H1 << separator;
        errorFile << std::scientific << post_process_data.norm_L2 << separator;
        errorFile << std::scientific << post_process_data.norm_H1 << separator;
        errorFile << std::scientific << assembler_data.globalMatrixA.NonZeros() << separator;
        errorFile << std::scientific << post_process_data.residual_norm << std::endl;

        errorFile.close();
    }

    {
        {
            Gedim::VTKUtilities exporter;
            exporter.AddSegments(mesh.Cell0DsCoordinates(),
                                 mesh.Cell1DsExtremes(),
                                 {{"Numeric",
                                   Gedim::VTPProperty::Formats::Points,
                                   static_cast<unsigned int>(post_process_data.cell0Ds_numeric.size()),
                                   post_process_data.cell0Ds_numeric.data()},
                                  {"Exact",
                                   Gedim::VTPProperty::Formats::Points,
                                   static_cast<unsigned int>(post_process_data.cell0Ds_exact.size()),
                                   post_process_data.cell0Ds_exact.data()},
                                  {"ErrorL2",
                                   Gedim::VTPProperty::Formats::Cells,
                                   static_cast<unsigned int>(post_process_data.cell1Ds_error_L2.size()),
                                   post_process_data.cell1Ds_error_L2.data()},
                                  {"ErrorH1",
                                   Gedim::VTPProperty::Formats::Cells,
                                   static_cast<unsigned int>(post_process_data.cell1Ds_error_H1.size()),
                                   post_process_data.cell1Ds_error_H1.data()}});

            exporter.Export(exportVtuFolder + "/Solution_" + std::to_string(TEST_ID) + "_" + std::to_string(TEST_ID) +
                            +"_" + std::to_string(config.MethodOrder()) + ".vtu");
        }
    }
}
// ***************************************************************************
void export_dofs(const Polydim::examples::Elliptic_PCC_1D::Program_configuration &config,
                 const Gedim::MeshMatricesDAO &mesh,
                 const Gedim::MeshUtilities::MeshGeometricData1D &mesh_geometric_data,
                 const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                 const FEM::PCC::FEM_PCC_1D_ReferenceElement_Data &vem_reference_element_data,
                 const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                 const Polydim::examples::Elliptic_PCC_1D::Assembler::Elliptic_PCC_1D_Problem_Data &assembler_data,
                 const Polydim::examples::Elliptic_PCC_1D::Assembler::PostProcess_Data &post_process_data,
                 const std::string &exportVtuFolder)
{
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = config.GeometricTolerance1D();
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    std::list<Eigen::Vector3d> dofs_coordinate;
    std::list<double> solution_values;
    std::list<double> rhs_values;
    std::list<double> dof_global_index_values;
    std::list<double> dof_type_values;
    std::list<double> dof_cell_index_values;
    std::list<double> dof_dimension_values;
    std::list<double> dof_boundary_type_values;
    std::list<double> dof_boundary_marker_values;

    for (unsigned int c = 0; c < mesh.Cell0DTotalNumber(); ++c)
    {
        const auto &boundary_info = mesh_dofs_info.CellsBoundaryInfo.at(0).at(c);

        const auto &local_dofs = dofs_data.CellsDOFs[0].at(c);

        const unsigned int num_loc_dofs = local_dofs.size();

        if (num_loc_dofs == 0)
            continue;

        for (unsigned int loc_i = 0; loc_i < num_loc_dofs; ++loc_i)
        {
            const auto &local_dof = local_dofs.at(loc_i);

            dof_cell_index_values.push_back(c);
            dof_dimension_values.push_back(0);
            dof_boundary_type_values.push_back(static_cast<double>(boundary_info.Type));
            dof_boundary_marker_values.push_back(boundary_info.Marker);
            dofs_coordinate.push_back(mesh.Cell0DCoordinates(c));
            dof_type_values.push_back(static_cast<double>(local_dof.Type));
            dof_global_index_values.push_back(local_dof.Global_Index);

            switch (local_dof.Type)
            {
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
                solution_values.push_back(assembler_data.solutionDirichlet.GetValue(local_dof.Global_Index));
                rhs_values.push_back(std::nan(""));
                break;
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                solution_values.push_back(assembler_data.solution.GetValue(local_dof.Global_Index));
                rhs_values.push_back(assembler_data.rightHandSide.GetValue(local_dof.Global_Index));
                break;
            default:
                throw std::runtime_error("Unknown DOF Type");
            }
        }
    }

    for (unsigned int c = 0; c < mesh.Cell1DTotalNumber(); ++c)
    {
        const auto &boundary_info = mesh_dofs_info.CellsBoundaryInfo.at(1).at(c);

        const auto &local_dofs = dofs_data.CellsDOFs[1].at(c);

        const unsigned int num_loc_dofs = local_dofs.size();

        if (num_loc_dofs == 0)
            continue;

        const std::vector<double> local_edge_coordinates = geometryUtilities.EquispaceCoordinates(num_loc_dofs, 0.0, 1.0, false);
        const Eigen::Vector3d edge_origin = mesh.Cell1DOriginCoordinates(c);
        const Eigen::Vector3d edge_tangent = mesh.Cell1DEndCoordinates(c) - edge_origin;

        for (unsigned int loc_i = 0; loc_i < num_loc_dofs; ++loc_i)
        {
            const auto &local_dof = local_dofs.at(loc_i);

            dof_cell_index_values.push_back(c);
            dof_dimension_values.push_back(1);
            dof_boundary_type_values.push_back(static_cast<double>(boundary_info.Type));
            dof_boundary_marker_values.push_back(boundary_info.Marker);
            dofs_coordinate.push_back(edge_origin + local_edge_coordinates[loc_i] * edge_tangent);
            dof_type_values.push_back(static_cast<double>(local_dof.Type));
            dof_global_index_values.push_back(local_dof.Global_Index);

            switch (local_dof.Type)
            {
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
                solution_values.push_back(assembler_data.solutionDirichlet.GetValue(local_dof.Global_Index));
                rhs_values.push_back(std::nan(""));
                break;
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                solution_values.push_back(assembler_data.solution.GetValue(local_dof.Global_Index));
                rhs_values.push_back(assembler_data.rightHandSide.GetValue(local_dof.Global_Index));
                break;
            default:
                throw std::runtime_error("Unknown DOF Type");
            }
        }
    }

    {
        Eigen::MatrixXd coordinates(3, dofs_coordinate.size());
        unsigned int c = 0;
        for (const auto &dof_coordinate : dofs_coordinate)
            coordinates.col(c++) << dof_coordinate;
        const auto rhs_values_data = std::vector<double>(rhs_values.begin(), rhs_values.end());
        const auto solution_values_data = std::vector<double>(solution_values.begin(), solution_values.end());
        const auto dof_global_index_values_data =
            std::vector<double>(dof_global_index_values.begin(), dof_global_index_values.end());
        const auto dof_type_values_data = std::vector<double>(dof_type_values.begin(), dof_type_values.end());
        const auto dof_cell_index_values_data = std::vector<double>(dof_cell_index_values.begin(), dof_cell_index_values.end());
        const auto dof_dimension_values_data = std::vector<double>(dof_dimension_values.begin(), dof_dimension_values.end());
        const auto dof_boundary_type_values_data =
            std::vector<double>(dof_boundary_type_values.begin(), dof_boundary_type_values.end());
        const auto dof_boundary_marker_values_data =
            std::vector<double>(dof_boundary_marker_values.begin(), dof_boundary_marker_values.end());

        Gedim::VTKUtilities exporter;
        exporter.AddPoints(coordinates,
                           {{"cell_dimension",
                             Gedim::VTPProperty::Formats::Points,
                             static_cast<unsigned int>(dof_dimension_values_data.size()),
                             dof_dimension_values_data.data()},
                            {"cell_index",
                             Gedim::VTPProperty::Formats::Points,
                             static_cast<unsigned int>(dof_cell_index_values_data.size()),
                             dof_cell_index_values_data.data()},
                            {"boundary_type",
                             Gedim::VTPProperty::Formats::Points,
                             static_cast<unsigned int>(dof_boundary_type_values_data.size()),
                             dof_boundary_type_values_data.data()},
                            {"boundary_marker",
                             Gedim::VTPProperty::Formats::Points,
                             static_cast<unsigned int>(dof_boundary_marker_values_data.size()),
                             dof_boundary_marker_values_data.data()},
                            {"dof_global_index",
                             Gedim::VTPProperty::Formats::Points,
                             static_cast<unsigned int>(dof_global_index_values_data.size()),
                             dof_global_index_values_data.data()},
                            {"dof_type",
                             Gedim::VTPProperty::Formats::Points,
                             static_cast<unsigned int>(dof_type_values_data.size()),
                             dof_type_values_data.data()},
                            {"rhs",
                             Gedim::VTPProperty::Formats::Points,
                             static_cast<unsigned int>(rhs_values_data.size()),
                             rhs_values_data.data()},
                            {"solution",
                             Gedim::VTPProperty::Formats::Points,
                             static_cast<unsigned int>(solution_values_data.size()),
                             solution_values_data.data()}});

        const unsigned int METHOD_ID = static_cast<unsigned int>(config.MethodType());
        const unsigned int TEST_ID = static_cast<unsigned int>(config.TestType());
        exporter.Export(exportVtuFolder + "/dofs_" + std::to_string(TEST_ID) + "_" + std::to_string(METHOD_ID) + +"_" +
                        std::to_string(config.MethodOrder()) + ".vtu");
    }
}
// ***************************************************************************
} // namespace program_utilities
} // namespace Elliptic_PCC_1D
} // namespace examples
} // namespace Polydim
