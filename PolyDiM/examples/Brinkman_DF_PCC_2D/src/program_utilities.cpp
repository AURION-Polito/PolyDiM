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

#include "VTKUtilities.hpp"

namespace Polydim
{
namespace examples
{
namespace Brinkman_DF_PCC_2D
{

unsigned int Polydim::examples::Brinkman_DF_PCC_2D::test::Patch_Test::order;

namespace program_utilities
{
// ***************************************************************************
std::unique_ptr<Polydim::examples::Brinkman_DF_PCC_2D::test::I_Test> create_test(const Polydim::examples::Brinkman_DF_PCC_2D::Program_configuration &config)
{
    switch (config.TestType())
    {
    case Polydim::examples::Brinkman_DF_PCC_2D::test::Test_Types::Patch_Test:
        Polydim::examples::Brinkman_DF_PCC_2D::test::Patch_Test::order = config.MethodOrder();
        return std::make_unique<Polydim::examples::Brinkman_DF_PCC_2D::test::Patch_Test>();
    case Polydim::examples::Brinkman_DF_PCC_2D::test::Test_Types::StokesSinSin:
        return std::make_unique<Polydim::examples::Brinkman_DF_PCC_2D::test::StokesSinSin>();
    case Polydim::examples::Brinkman_DF_PCC_2D::test::Test_Types::Stokes_ZeroVelocity_1:
        return std::make_unique<Polydim::examples::Brinkman_DF_PCC_2D::test::Stokes_ZeroVelocity_1>();
    case Polydim::examples::Brinkman_DF_PCC_2D::test::Test_Types::Stokes_ZeroVelocity_2:
        return std::make_unique<Polydim::examples::Brinkman_DF_PCC_2D::test::Stokes_ZeroVelocity_2>();
    case Polydim::examples::Brinkman_DF_PCC_2D::test::Test_Types::Darcy:
        return std::make_unique<Polydim::examples::Brinkman_DF_PCC_2D::test::Darcy>();
    case Polydim::examples::Brinkman_DF_PCC_2D::test::Test_Types::Brinkman:
        return std::make_unique<Polydim::examples::Brinkman_DF_PCC_2D::test::Brinkman>();
    case Polydim::examples::Brinkman_DF_PCC_2D::test::Test_Types::DarcyStokes_1:
        return std::make_unique<Polydim::examples::Brinkman_DF_PCC_2D::test::DarcyStokes_1>();
    case Polydim::examples::Brinkman_DF_PCC_2D::test::Test_Types::DarcyStokes_2:
        return std::make_unique<Polydim::examples::Brinkman_DF_PCC_2D::test::DarcyStokes_2>();
    default:
        throw std::runtime_error("Test type " + std::to_string((unsigned int)config.TestType()) + " not supported");
    }
}
// ***************************************************************************
void create_domain_mesh(const Polydim::examples::Brinkman_DF_PCC_2D::Program_configuration &config,
                        const Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D &domain,
                        Gedim::MeshMatricesDAO &mesh)
{
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = config.GeometricTolerance1D();
    geometryUtilitiesConfig.Tolerance2D = config.GeometricTolerance2D();
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::MeshUtilities meshUtilities;

    switch (config.MeshGenerator())
    {
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::Triangular:
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::Minimal:
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::Polygonal:
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::Squared: {
        Polydim::PDETools::Mesh::PDE_Mesh_Utilities::create_mesh_2D(geometryUtilities,
                                                                    meshUtilities,
                                                                    config.MeshGenerator(),
                                                                    domain,
                                                                    config.MeshMaxArea(),
                                                                    mesh);
    }
    break;
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::OFFImporter:
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::CsvImporter:
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::TriangularSimpleImporter: {
        Polydim::PDETools::Mesh::PDE_Mesh_Utilities::import_mesh_2D(geometryUtilities,
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
Gedim::MeshUtilities::MeshGeometricData2D create_domain_mesh_geometric_properties(const Polydim::examples::Brinkman_DF_PCC_2D::Program_configuration &config,
                                                                                  Gedim::MeshMatricesDAO &mesh)
{
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = config.GeometricTolerance1D();
    geometryUtilitiesConfig.Tolerance2D = config.GeometricTolerance2D();
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::MeshUtilities meshUtilities;
    meshUtilities.ComputeCell1DCell2DNeighbours(mesh);
    return Polydim::PDETools::Mesh::PDE_Mesh_Utilities::compute_mesh_2D_geometry_data(geometryUtilities, meshUtilities, mesh);
}
// ***************************************************************************
void export_solution(const Polydim::examples::Brinkman_DF_PCC_2D::Program_configuration &config,
                     const Gedim::MeshMatricesDAO &mesh,
                     const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
                     const Polydim::PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                     const Polydim::examples::Brinkman_DF_PCC_2D::Assembler::Stokes_DF_PCC_2D_Problem_Data &assembler_data,
                     const Polydim::examples::Brinkman_DF_PCC_2D::Assembler::PostProcess_Data &post_process_data,
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
        std::cout << "Cell2Ds" << separator;
        std::cout << "Dofs" << separator;
        std::cout << "Strongs" << separator;
        std::cout << "h" << separator;
        std::cout << "errorH1Velocity" << separator;
        std::cout << "errorL2Pressure" << separator;
        std::cout << "normH1Velocity" << separator;
        std::cout << "normL2Pressure" << separator;
        std::cout << "nnzA" << separator;
        std::cout << "residual" << std::endl;

        std::cout.precision(2);
        std::cout << std::scientific << TEST_ID << separator;
        std::cout << std::scientific << METHOD_ID << separator;
        std::cout << std::scientific << config.MethodOrder() << separator;
        std::cout << std::scientific << mesh.Cell2DTotalNumber() << separator;
        std::cout << std::scientific << count_dofs.num_total_dofs << separator;
        std::cout << std::scientific << count_dofs.num_total_strong << separator;
        std::cout << std::scientific << post_process_data.mesh_size << separator;
        std::cout << std::scientific << post_process_data.error_H1_velocity << separator;
        std::cout << std::scientific << post_process_data.error_L2_pressure << separator;
        std::cout << std::scientific << post_process_data.norm_H1_velocity << separator;
        std::cout << std::scientific << post_process_data.norm_L2_pressure << separator;
        std::cout << std::scientific << assembler_data.globalMatrixA.NonZeros() << separator;
        std::cout << std::scientific << post_process_data.residual_norm << std::endl;
    }

    {
        const char separator = ';';
        const std::string errorFileName = exportSolutionFolder + "/Errors_" + std::to_string(TEST_ID) + "_" +
                                          std::to_string(METHOD_ID) + "_" + std::to_string(config.MethodOrder()) + ".csv";
        const bool errorFileExists = Gedim::Output::FileExists(errorFileName);

        std::ofstream errorFile(errorFileName, std::ios_base::app | std::ios_base::out);

        if (!errorFileExists)
        {
            errorFile << "ProgramType" << separator;
            errorFile << "MethodType" << separator;
            errorFile << "MethodOrder" << separator;
            errorFile << "Cell2Ds" << separator;
            errorFile << "Dofs" << separator;
            errorFile << "Strongs" << separator;
            errorFile << "h" << separator;
            errorFile << "errorH1Velocity" << separator;
            errorFile << "errorL2Pressure" << separator;
            errorFile << "normH1Velocity" << separator;
            errorFile << "normL2Pressure" << separator;
            errorFile << "nnzA" << separator;
            errorFile << "residual" << std::endl;
        }

        errorFile.precision(16);
        errorFile << std::scientific << TEST_ID << separator;
        errorFile << std::scientific << METHOD_ID << separator;
        errorFile << std::scientific << config.MethodOrder() << separator;
        errorFile << std::scientific << mesh.Cell2DTotalNumber() << separator;
        errorFile << std::scientific << count_dofs.num_total_dofs << separator;
        errorFile << std::scientific << count_dofs.num_total_strong << separator;
        errorFile << std::scientific << post_process_data.mesh_size << separator;
        errorFile << std::scientific << post_process_data.error_H1_velocity << separator;
        errorFile << std::scientific << post_process_data.error_L2_pressure << separator;
        errorFile << std::scientific << post_process_data.norm_H1_velocity << separator;
        errorFile << std::scientific << post_process_data.norm_L2_pressure << separator;
        errorFile << std::scientific << assembler_data.globalMatrixA.NonZeros() << separator;
        errorFile << std::scientific << post_process_data.residual_norm << std::endl;

        errorFile.close();
    }

    {
        const char separator = ';';
        const std::string fluxFileName = exportSolutionFolder + "/Flux_" + std::to_string(TEST_ID) + "_" +
                                         std::to_string(METHOD_ID) + "_" + std::to_string(config.MethodOrder()) + ".csv";

        std::ofstream fluxFile(fluxFileName, std::ios_base::out);

        fluxFile << "EdgeMarker" << separator;
        fluxFile << "Flux" << std::endl;

        fluxFile.precision(16);
        double sum = 0.0;
        for (const auto &f : post_process_data.flux)
        {
            fluxFile << std::scientific << f.first << separator << f.second << std::endl;
            sum += f.second;
        }

        std::cout.precision(16);
        std::cout << std::scientific << "sum flux: " << sum << std::endl;

        fluxFile.close();
    }

    if (!std::isnan(post_process_data.cell0Ds_exact_velocity[0](0)))
    {
        {
            Gedim::VTKUtilities exporter;
            exporter.AddPolygons(mesh.Cell0DsCoordinates(),
                                 mesh.Cell2DsVertices(),
                                 {{"Numeric Velocity - X",
                                   Gedim::VTPProperty::Formats::Points,
                                   static_cast<unsigned int>(post_process_data.cell0Ds_numeric_velocity[0].size()),
                                   post_process_data.cell0Ds_numeric_velocity[0].data()},
                                  {"Numeric Velocity - Y",
                                   Gedim::VTPProperty::Formats::Points,
                                   static_cast<unsigned int>(post_process_data.cell0Ds_numeric_velocity[1].size()),
                                   post_process_data.cell0Ds_numeric_velocity[1].data()},
                                  {"Exact Velocity - X",
                                   Gedim::VTPProperty::Formats::Points,
                                   static_cast<unsigned int>(post_process_data.cell0Ds_exact_velocity[0].size()),
                                   post_process_data.cell0Ds_exact_velocity[0].data()},
                                  {"Exact Velocity - Y",
                                   Gedim::VTPProperty::Formats::Points,
                                   static_cast<unsigned int>(post_process_data.cell0Ds_exact_velocity[1].size()),
                                   post_process_data.cell0Ds_exact_velocity[1].data()},
                                  {"ErrorH1Velocity",
                                   Gedim::VTPProperty::Formats::Cells,
                                   static_cast<unsigned int>(post_process_data.cell2Ds_error_H1_velocity.size()),
                                   post_process_data.cell2Ds_error_H1_velocity.data()},
                                  {"Permeability",
                                   Gedim::VTPProperty::Formats::Cells,
                                   static_cast<unsigned int>(post_process_data.inverse_diffusion_coeff_values.size()),
                                   post_process_data.inverse_diffusion_coeff_values.data()},
                                  {"Viscosity",
                                   Gedim::VTPProperty::Formats::Cells,
                                   static_cast<unsigned int>(post_process_data.viscosity_values.size()),
                                   post_process_data.viscosity_values.data()}});

            exporter.Export(exportVtuFolder + "/Velocity_" + std::to_string(TEST_ID) + "_" + std::to_string(METHOD_ID) +
                            +"_" + std::to_string(config.MethodOrder()) + ".vtu");
        }

        {
            Gedim::VTKUtilities exporter;
            exporter.AddPolygons(post_process_data.repeated_vertices_coordinates,
                                 post_process_data.repeated_connectivity,
                                 {{"Numeric Pressure",
                                   Gedim::VTPProperty::Formats::Points,
                                   static_cast<unsigned int>(post_process_data.cell0Ds_numeric_pressure.size()),
                                   post_process_data.cell0Ds_numeric_pressure.data()},
                                  {"Exact Pressure",
                                   Gedim::VTPProperty::Formats::Points,
                                   static_cast<unsigned int>(post_process_data.cell0Ds_exact_pressure.size()),
                                   post_process_data.cell0Ds_exact_pressure.data()},
                                  {"ErrorL2Pressure",
                                   Gedim::VTPProperty::Formats::Cells,
                                   static_cast<unsigned int>(post_process_data.cell2Ds_error_L2_pressure.size()),
                                   post_process_data.cell2Ds_error_L2_pressure.data()}});

            exporter.Export(exportVtuFolder + "/Pressure_" + std::to_string(TEST_ID) + "_" + std::to_string(METHOD_ID) +
                            +"_" + std::to_string(config.MethodOrder()) + ".vtu");
        }
    }
    else
    {
        {
            Gedim::VTKUtilities exporter;
            exporter.AddPolygons(mesh.Cell0DsCoordinates(),
                                 mesh.Cell2DsVertices(),
                                 {{"Numeric Velocity - X",
                                   Gedim::VTPProperty::Formats::Points,
                                   static_cast<unsigned int>(post_process_data.cell0Ds_numeric_velocity[0].size()),
                                   post_process_data.cell0Ds_numeric_velocity[0].data()},
                                  {"Numeric Velocity - Y",
                                   Gedim::VTPProperty::Formats::Points,
                                   static_cast<unsigned int>(post_process_data.cell0Ds_numeric_velocity[1].size()),
                                   post_process_data.cell0Ds_numeric_velocity[1].data()},
                                  {"Permeability",
                                   Gedim::VTPProperty::Formats::Cells,
                                   static_cast<unsigned int>(post_process_data.inverse_diffusion_coeff_values.size()),
                                   post_process_data.inverse_diffusion_coeff_values.data()},
                                  {"Viscosity",
                                   Gedim::VTPProperty::Formats::Cells,
                                   static_cast<unsigned int>(post_process_data.viscosity_values.size()),
                                   post_process_data.viscosity_values.data()}});

            exporter.Export(exportVtuFolder + "/Velocity_" + std::to_string(TEST_ID) + "_" + std::to_string(METHOD_ID) +
                            +"_" + std::to_string(config.MethodOrder()) + ".vtu");
        }

        {
            Gedim::VTKUtilities exporter;
            exporter.AddPolygons(post_process_data.repeated_vertices_coordinates,
                                 post_process_data.repeated_connectivity,
                                 {{"Numeric Pressure",
                                   Gedim::VTPProperty::Formats::Points,
                                   static_cast<unsigned int>(post_process_data.cell0Ds_numeric_pressure.size()),
                                   post_process_data.cell0Ds_numeric_pressure.data()}});

            exporter.Export(exportVtuFolder + "/Pressure_" + std::to_string(TEST_ID) + "_" + std::to_string(METHOD_ID) +
                            +"_" + std::to_string(config.MethodOrder()) + ".vtu");
        }
    }
}
// ***************************************************************************
void export_velocity_dofs(const Polydim::examples::Brinkman_DF_PCC_2D::Program_configuration &config,
                          const Gedim::MeshMatricesDAO &mesh,
                          const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                          const std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> &mesh_dofs_info,
                          const Polydim::PDETools::LocalSpace_DF_PCC_2D::ReferenceElement_Data &reference_element_data,
                          const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
                          const Polydim::PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                          const Polydim::examples::Brinkman_DF_PCC_2D::Assembler::Stokes_DF_PCC_2D_Problem_Data &assembler_data,
                          const std::string &exportVtuFolder)
{
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = config.GeometricTolerance1D();
    geometryUtilitiesConfig.Tolerance2D = config.GeometricTolerance2D();
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    const unsigned int METHOD_ID = static_cast<unsigned int>(config.MethodType());
    const unsigned int TEST_ID = static_cast<unsigned int>(config.TestType());

    const std::string file_path = exportVtuFolder + "/dofs_" + std::to_string(TEST_ID) + "_" +
                                  std::to_string(METHOD_ID) + +"_" + std::to_string(config.MethodOrder()) + ".vtu";

    Polydim::PDETools::LocalSpace_DF_PCC_2D::export_velocity_dofs(geometryUtilities,
                                                                  mesh,
                                                                  mesh_geometric_data,
                                                                  mesh_dofs_info,
                                                                  dofs_data,
                                                                  count_dofs,
                                                                  assembler_data.rightHandSide,
                                                                  assembler_data.solution,
                                                                  assembler_data.solutionDirichlet,
                                                                  file_path);
}
// ***************************************************************************
void export_discrepancy_errors(const Polydim::examples::Brinkman_DF_PCC_2D::Program_configuration &config,
                               const Gedim::MeshMatricesDAO &mesh,
                               const Polydim::examples::Brinkman_DF_PCC_2D::Assembler::DiscrepancyErrors_Data &discrepancy_errors_data,
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
        std::cout << "Cell2Ds" << separator;
        std::cout << "VelocityDofsRatio" << separator;
        std::cout << "PressureDofsRatio" << separator;
        std::cout << "discrepancyErrorH1Velocity" << separator;
        std::cout << "discrepancyErrorL2Pressure" << separator;
        std::cout << "normH1FULLVelocity" << separator;
        std::cout << "normL2FULLProjectedPressure" << separator;
        std::cout << "reducedResidual" << separator;
        std::cout << "fullResidual" << std::endl;

        std::cout.precision(2);
        std::cout << std::scientific << TEST_ID << separator;
        std::cout << std::scientific << METHOD_ID << separator;
        std::cout << std::scientific << config.MethodOrder() << separator;
        std::cout << std::scientific << mesh.Cell2DTotalNumber() << separator;
        std::cout << std::scientific << discrepancy_errors_data.velocity_dofs_ratio << separator;
        std::cout << std::scientific << discrepancy_errors_data.pressure_dofs_ratio << separator;
        std::cout << std::scientific << discrepancy_errors_data.discrepancy_error_H1_velocity << separator;
        std::cout << std::scientific << discrepancy_errors_data.discrepancy_error_L2_pressure << separator;
        std::cout << std::scientific << discrepancy_errors_data.full_norm_H1_velocity << separator;
        std::cout << std::scientific << discrepancy_errors_data.full_norm_L2_pressure << separator;
        std::cout << std::scientific << discrepancy_errors_data.reduced_residual_norm << separator;
        std::cout << std::scientific << discrepancy_errors_data.residual_norm << std::endl;
    }

    {
        const char separator = ';';
        const std::string errorFileName = exportSolutionFolder + "/DiscrepancyErrors_" + std::to_string(TEST_ID) + "_" +
                                          std::to_string(METHOD_ID) + +"_" + std::to_string(config.MethodOrder()) + ".csv";
        const bool errorFileExists = Gedim::Output::FileExists(errorFileName);

        std::ofstream errorFile(errorFileName, std::ios_base::app | std::ios_base::out);

        if (!errorFileExists)
        {
            errorFile << "ProgramType" << separator;
            errorFile << "MethodType" << separator;
            errorFile << "MethodOrder" << separator;
            errorFile << "Cell2Ds" << separator;
            errorFile << "VelocityDofsRatio" << separator;
            errorFile << "PressureDofsRatio" << separator;
            errorFile << "discrepancyErrorH1Velocity" << separator;
            errorFile << "discrepancyErrorL2Pressure" << separator;
            errorFile << "normH1FULLVelocity" << separator;
            errorFile << "normL2FULLProjectedPressure" << separator;
            errorFile << "reducedResidual" << separator;
            errorFile << "fullResidual" << std::endl;
        }

        errorFile.precision(16);
        errorFile << std::scientific << TEST_ID << separator;
        errorFile << std::scientific << METHOD_ID << separator;
        errorFile << std::scientific << config.MethodOrder() << separator;
        errorFile << std::scientific << mesh.Cell2DTotalNumber() << separator;
        errorFile << std::scientific << discrepancy_errors_data.velocity_dofs_ratio << separator;
        errorFile << std::scientific << discrepancy_errors_data.pressure_dofs_ratio << separator;
        errorFile << std::scientific << discrepancy_errors_data.discrepancy_error_H1_velocity << separator;
        errorFile << std::scientific << discrepancy_errors_data.discrepancy_error_L2_pressure << separator;
        errorFile << std::scientific << discrepancy_errors_data.full_norm_H1_velocity << separator;
        errorFile << std::scientific << discrepancy_errors_data.full_norm_L2_pressure << separator;
        errorFile << std::scientific << discrepancy_errors_data.reduced_residual_norm << separator;
        errorFile << std::scientific << discrepancy_errors_data.residual_norm << std::endl;

        errorFile.close();
    }

    {
        {
            Gedim::VTKUtilities exporter;
            exporter.AddPolygons(
                mesh.Cell0DsCoordinates(),
                mesh.Cell2DsVertices(),
                {{"ErrorL2Pressure",
                  Gedim::VTPProperty::Formats::Cells,
                  static_cast<unsigned int>(discrepancy_errors_data.cell2Ds_discrepancy_error_L2_pressure.size()),
                  discrepancy_errors_data.cell2Ds_discrepancy_error_L2_pressure.data()},
                 {"ErrorH1Velocity",
                  Gedim::VTPProperty::Formats::Cells,
                  static_cast<unsigned int>(discrepancy_errors_data.cell2Ds_discrepancy_error_H1_velocity.size()),
                  discrepancy_errors_data.cell2Ds_discrepancy_error_H1_velocity.data()}});

            exporter.Export(exportVtuFolder + "/DiscrepancyErrors_" + std::to_string(TEST_ID) + "_" +
                            std::to_string(METHOD_ID) + +"_" + std::to_string(config.MethodOrder()) + ".vtu");
        }
    }
}
// ***************************************************************************
void export_performance(const Polydim::examples::Brinkman_DF_PCC_2D::Program_configuration &config,
                        const Assembler::Performance_Data &performance_data,
                        const std::string &exportFolder)
{

    const char separator = ',';
    std::ofstream exporter;
    const unsigned int Method_ID = static_cast<unsigned int>(config.MethodType());
    const unsigned int TEST_ID = static_cast<unsigned int>(config.TestType());
    exporter.open(exportFolder + "/Cell2Ds_MethodPerformance_" + std::to_string(TEST_ID) + "_" +
                  std::to_string(Method_ID) + "_" + std::to_string(config.MethodOrder()) + ".csv");
    exporter.precision(16);

    if (exporter.fail())
        throw std::runtime_error("Error on mesh cell2Ds file");

    switch (config.MethodType())
    {
    case PDETools::LocalSpace_DF_PCC_2D::MethodTypes::VEM_DF_PCC_FULL:
    case PDETools::LocalSpace_DF_PCC_2D::MethodTypes::VEM_DF_PCC_REDUCED: {

        exporter << "Cell2D_Index" << separator;
        exporter << "NumQuadPoints_Boundary" << separator;
        exporter << "NumQuadPoints_Internal" << separator;
        exporter << "max_PiNabla_Cond" << separator;
        exporter << "max_Pi0k_Cond" << separator;
        exporter << "max_PiNabla_Error" << separator;
        exporter << "max_Pi0k_Error" << separator;
        exporter << "max_GBD_Error" << separator;
        exporter << "max_HCD_Error" << separator;
        exporter << "Stab_Error" << std::endl;

        for (unsigned int v = 0; v < performance_data.Cell2DsPerformance.size(); v++)
        {
            const auto &cell2DPerformance = performance_data.Cell2DsPerformance[v].performance_data.vem_analysis_data;

            exporter << std::scientific << v << separator;
            exporter << std::scientific
                     << performance_data.Cell2DsPerformance[v].performance_data.NumBoundaryQuadraturePoints << separator;
            exporter << std::scientific
                     << performance_data.Cell2DsPerformance[v].performance_data.NumInternalQuadraturePoints << separator;
            exporter << std::scientific << cell2DPerformance.maxPiNablaConditioning << separator;
            exporter << std::scientific << cell2DPerformance.maxPi0kConditioning << separator;
            exporter << std::scientific << cell2DPerformance.maxErrorPiNabla << separator;
            exporter << std::scientific << cell2DPerformance.maxErrorPi0k << separator;
            exporter << std::scientific << cell2DPerformance.maxErrorGBD << separator;
            exporter << std::scientific << cell2DPerformance.maxErrorHCD << separator;
            exporter << std::scientific << cell2DPerformance.ErrorStabilization << std::endl;
        }
    }
    break;
    case PDETools::LocalSpace_DF_PCC_2D::MethodTypes::TAYLOR_HOOD:
        break;
    default:
        throw std::runtime_error("not valid method type");
    }

    exporter.close();
}
// ***************************************************************************
} // namespace program_utilities
} // namespace Brinkman_DF_PCC_2D
} // namespace examples
} // namespace Polydim
