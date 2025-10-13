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
namespace Elliptic_PCC_BulkFace_2D
{
namespace program_utilities
{
// ***************************************************************************
std::unique_ptr<Polydim::examples::Elliptic_PCC_BulkFace_2D::test::I_Test> create_test(
    const Polydim::examples::Elliptic_PCC_BulkFace_2D::Program_configuration &config)
{
    switch (config.TestType())
    {
    case Polydim::examples::Elliptic_PCC_BulkFace_2D::test::Test_Types::Patch_Test:
        return std::make_unique<Polydim::examples::Elliptic_PCC_BulkFace_2D::test::Patch_Test>();
    // case Polydim::examples::Elliptic_PCC_BulkFace_2D::test::Test_Types::Elliptic_Problem:
    //     return std::make_unique<Polydim::examples::Elliptic_PCC_BulkFace_2D::test::Elliptic_Problem>();
    default:
        throw std::runtime_error("Test type " + std::to_string((unsigned int)config.TestType()) + " not supported");
    }
}
// ***************************************************************************
void create_domain_mesh_2D(const Polydim::examples::Elliptic_PCC_BulkFace_2D::Program_configuration &config,
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
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::Squared:
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::RandomDistorted: {
        Polydim::PDETools::Mesh::PDE_Mesh_Utilities::create_mesh_2D(geometryUtilities,
                                                                    meshUtilities,
                                                                    config.MeshGenerator(),
                                                                    domain,
                                                                    config.MeshMaxArea(),
                                                                    mesh);
    }
    break;
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::CsvImporter:
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::OFFImporter: {
        Polydim::PDETools::Mesh::PDE_Mesh_Utilities::import_mesh_2D(meshUtilities, config.MeshGenerator(), config.MeshImportFilePath(), mesh);
    }
    break;
    default:
        throw std::runtime_error("MeshGenerator " + std::to_string((unsigned int)config.MeshGenerator()) + " not supported");
    }
}
// ***************************************************************************
std::vector<double> create_time_steps(const Polydim::examples::Elliptic_PCC_BulkFace_2D::Program_configuration &config,
                                      const std::array<double, 2> &time_domain)
{
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = config.GeometricTolerance1D();
    geometryUtilitiesConfig.Tolerance2D = config.GeometricTolerance2D();
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::Output::Assert(geometryUtilities.IsValuePositive(config.TimeStep(), geometryUtilities.Tolerance1D()));
    std::vector<double> times = geometryUtilities.EquispaceCoordinates(config.TimeStep(), true);
    for (unsigned int t = 0; t < times.size(); t++)
        times[t] = (time_domain.at(1) - time_domain.at(0)) * times[t] + time_domain.at(0);

    return times;
}
// ***************************************************************************
void export_solution(const Polydim::examples::Elliptic_PCC_BulkFace_2D::Program_configuration &config,
                     const Gedim::MeshMatricesDAO &mesh_2D,
                     const Gedim::MeshMatricesDAO &mesh_1D,
                     const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
                     const Polydim::examples::Elliptic_PCC_BulkFace_2D::Assembler::Elliptic_PCC_BF_2D_Problem_Data &assembler_data,
                     const Polydim::examples::Elliptic_PCC_BulkFace_2D::Assembler::PostProcess_Data &post_process_data,
                     const std::string &exportSolutionFolder,
                     const std::string &exportVtuFolder)
{

    const char separator = ';';

    export_solution_2D(config, mesh_2D, dofs_data[0], post_process_data.post_process_data_2D, exportSolutionFolder, exportVtuFolder);

    export_solution_1D(config, mesh_1D, dofs_data[1], post_process_data.post_process_data_1D, exportSolutionFolder, exportVtuFolder);

    std::cout << "nnzA" << separator;
    std::cout << "residual" << std::endl;

    std::cout << std::scientific << assembler_data.globalMatrixA.NonZeros() << separator;
    std::cout << std::scientific << post_process_data.residual_norm << std::endl;
}
// ***************************************************************************
void export_solution_1D(const Polydim::examples::Elliptic_PCC_BulkFace_2D::Program_configuration &config,
                        const Gedim::MeshMatricesDAO &mesh,
                        const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                        const Polydim::examples::Elliptic_PCC_BulkFace_2D::Assembler::PostProcess_Data_1D &post_process_data,
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
        std::cout << "normH1" << std::endl;

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
        std::cout << std::scientific << post_process_data.norm_H1 << std::endl;
    }

    {
        const char separator = ';';
        const std::string errorFileName = exportSolutionFolder + "/Errors_1D_" + std::to_string(TEST_ID) + "_" +
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
            errorFile << "normH1" << std::endl;
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
        errorFile << std::scientific << post_process_data.norm_H1 << std::endl;

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

            exporter.Export(exportVtuFolder + "/Solution_1D_" + std::to_string(TEST_ID) + "_" +
                            std::to_string(TEST_ID) + +"_" + std::to_string(config.MethodOrder()) + ".vtu");
        }
    }
}
// ***************************************************************************
void export_solution_2D(const Polydim::examples::Elliptic_PCC_BulkFace_2D::Program_configuration &config,
                        const Gedim::MeshMatricesDAO &mesh,
                        const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                        const Polydim::examples::Elliptic_PCC_BulkFace_2D::Assembler::PostProcess_Data_2D &post_process_data,
                        const std::string &exportSolutionFolder,
                        const std::string &exportVtuFolder)
{
    const unsigned int Method_ID = static_cast<unsigned int>(config.MethodType());
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
        std::cout << "errorL2" << separator;
        std::cout << "errorH1" << separator;
        std::cout << "normL2" << separator;
        std::cout << "normH1" << std::endl;

        std::cout.precision(2);
        std::cout << std::scientific << TEST_ID << separator;
        std::cout << std::scientific << Method_ID << separator;
        std::cout << std::scientific << config.MethodOrder() << separator;
        std::cout << std::scientific << mesh.Cell2DTotalNumber() << separator;
        std::cout << std::scientific << dofs_data.NumberDOFs << separator;
        std::cout << std::scientific << dofs_data.NumberStrongs << separator;
        std::cout << std::scientific << post_process_data.mesh_size << separator;
        std::cout << std::scientific << post_process_data.error_L2 << separator;
        std::cout << std::scientific << post_process_data.error_H1 << separator;
        std::cout << std::scientific << post_process_data.norm_L2 << separator;
        std::cout << std::scientific << post_process_data.norm_H1 << std::endl;
    }

    {
        const char separator = ';';
        const std::string errorFileName = exportSolutionFolder + "/Errors_2D_" + std::to_string(TEST_ID) + "_" +
                                          std::to_string(Method_ID) + +"_" + std::to_string(config.MethodOrder()) + ".csv";
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
            errorFile << "errorL2" << separator;
            errorFile << "errorH1" << separator;
            errorFile << "normL2" << separator;
            errorFile << "normH1" << std::endl;
        }

        errorFile.precision(16);
        errorFile << std::scientific << TEST_ID << separator;
        errorFile << std::scientific << Method_ID << separator;
        errorFile << std::scientific << config.MethodOrder() << separator;
        errorFile << std::scientific << mesh.Cell2DTotalNumber() << separator;
        errorFile << std::scientific << dofs_data.NumberDOFs << separator;
        errorFile << std::scientific << dofs_data.NumberStrongs << separator;
        errorFile << std::scientific << post_process_data.mesh_size << separator;
        errorFile << std::scientific << post_process_data.error_L2 << separator;
        errorFile << std::scientific << post_process_data.error_H1 << separator;
        errorFile << std::scientific << post_process_data.norm_L2 << separator;
        errorFile << std::scientific << post_process_data.norm_H1 << std::endl;

        errorFile.close();
    }

    {
        {
            Gedim::VTKUtilities exporter;
            exporter.AddPolygons(mesh.Cell0DsCoordinates(),
                                 mesh.Cell2DsVertices(),
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
                                   static_cast<unsigned int>(post_process_data.cell2Ds_error_L2.size()),
                                   post_process_data.cell2Ds_error_L2.data()},
                                  {"ErrorH1",
                                   Gedim::VTPProperty::Formats::Cells,
                                   static_cast<unsigned int>(post_process_data.cell2Ds_error_H1.size()),
                                   post_process_data.cell2Ds_error_H1.data()}});

            exporter.Export(exportVtuFolder + "/Solution_2D_" + std::to_string(TEST_ID) + "_" +
                            std::to_string(Method_ID) + +"_" + std::to_string(config.MethodOrder()) + ".vtu");
        }
    }
}
// ***************************************************************************
void export_performance_2D(const Polydim::examples::Elliptic_PCC_BulkFace_2D::Program_configuration &config,
                           const Assembler::Performance_Data_2D &performance_data,
                           const std::string &exportFolder)
{
    {
        const char separator = ',';
        std::ofstream exporter;
        const unsigned int Method_ID = static_cast<unsigned int>(config.MethodType());
        const unsigned int TEST_ID = static_cast<unsigned int>(config.TestType());
        exporter.open(exportFolder + "/Cell2Ds_MethodPerformance_" + std::to_string(TEST_ID) + "_" +
                      std::to_string(Method_ID) + +"_" + std::to_string(config.MethodOrder()) + ".csv");
        exporter.precision(16);

        if (exporter.fail())
            throw std::runtime_error("Error on mesh cell2Ds file");

        exporter << "Cell2D_Index" << separator;
        exporter << "NumQuadPoints_Boundary" << separator;
        exporter << "NumQuadPoints_Internal" << separator;
        exporter << "PiNabla_Cond" << separator;
        exporter << "Pi0k_Cond" << separator;
        exporter << "Pi0km1_Cond" << separator;
        exporter << "PiNabla_Error" << separator;
        exporter << "Pi0k_Error" << separator;
        exporter << "Pi0km1_Error" << separator;
        exporter << "HCD_Error" << separator;
        exporter << "GBD_Error" << separator;
        exporter << "Stab_Error" << std::endl;

        for (unsigned int v = 0; v < performance_data.Cell2DsPerformance.size(); v++)
        {
            const auto &cell2D_performance = performance_data.Cell2DsPerformance[v].VEM_Performance_Data;

            exporter << std::scientific << v << separator;
            exporter << std::scientific << cell2D_performance.NumBoundaryQuadraturePoints << separator;
            exporter << std::scientific << cell2D_performance.NumInternalQuadraturePoints << separator;
            exporter << std::scientific << cell2D_performance.Analysis.PiNablaConditioning << separator;
            exporter << std::scientific << cell2D_performance.Analysis.Pi0kConditioning << separator;
            exporter << std::scientific << cell2D_performance.Analysis.Pi0km1Conditioning << separator;
            exporter << std::scientific << cell2D_performance.Analysis.ErrorPiNabla << separator;
            exporter << std::scientific << cell2D_performance.Analysis.ErrorPi0k << separator;
            exporter << std::scientific << cell2D_performance.Analysis.ErrorPi0km1 << separator;
            exporter << std::scientific << cell2D_performance.Analysis.ErrorHCD << separator;
            exporter << std::scientific << cell2D_performance.Analysis.ErrorGBD << separator;
            exporter << std::scientific << cell2D_performance.Analysis.ErrorStabilization << std::endl;
        }

        exporter.close();
    }
}
// ***************************************************************************
} // namespace program_utilities
} // namespace Elliptic_PCC_BulkFace_2D
} // namespace examples
} // namespace Polydim
