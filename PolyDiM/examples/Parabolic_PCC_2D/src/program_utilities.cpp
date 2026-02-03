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
#include <numbers>

namespace Polydim
{
namespace examples
{
namespace Parabolic_PCC_2D
{
namespace program_utilities
{
// ***************************************************************************
std::unique_ptr<Polydim::examples::Parabolic_PCC_2D::test::I_Test> create_test(const Polydim::examples::Parabolic_PCC_2D::Program_configuration &config)
{
    switch (config.TestType())
    {
    case Polydim::examples::Parabolic_PCC_2D::test::Test_Types::Patch_Test: {
        Polydim::examples::Parabolic_PCC_2D::test::Patch_Test::space_order = config.MethodOrder();
        Polydim::examples::Parabolic_PCC_2D::test::Patch_Test::time_order = config.Theta() == 0.5 ? 2 : 1;
        return std::make_unique<Polydim::examples::Parabolic_PCC_2D::test::Patch_Test>();
    }
    case Polydim::examples::Parabolic_PCC_2D::test::Test_Types::Space_Test: {
        Polydim::examples::Parabolic_PCC_2D::test::Patch_Test::space_order = config.MethodOrder() + 1;
        Polydim::examples::Parabolic_PCC_2D::test::Patch_Test::time_order = 0; // config.Theta() == 0.5 ? 2 : 1;
        return std::make_unique<Polydim::examples::Parabolic_PCC_2D::test::Patch_Test>();
    }
    case Polydim::examples::Parabolic_PCC_2D::test::Test_Types::Time_Test: {
        Polydim::examples::Parabolic_PCC_2D::test::Patch_Test::space_order = config.MethodOrder();
        Polydim::examples::Parabolic_PCC_2D::test::Patch_Test::time_order = config.Theta() == 0.5 ? 2 + 1 : 1 + 1;
        return std::make_unique<Polydim::examples::Parabolic_PCC_2D::test::Patch_Test>();
    }
    case Polydim::examples::Parabolic_PCC_2D::test::Test_Types::Parabolic_Problem:
        return std::make_unique<Polydim::examples::Parabolic_PCC_2D::test::Parabolic_Problem>();
    default:
        throw std::runtime_error("Test type " + std::to_string((unsigned int)config.TestType()) + " not supported");
    }
}
// ***************************************************************************
void create_domain_mesh(const Polydim::examples::Parabolic_PCC_2D::Program_configuration &config,
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
Gedim::MeshUtilities::MeshGeometricData2D create_domain_mesh_geometric_properties(const Polydim::examples::Parabolic_PCC_2D::Program_configuration &config,
                                                                                  const Gedim::MeshMatricesDAO &mesh)
{
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = config.GeometricTolerance1D();
    geometryUtilitiesConfig.Tolerance2D = config.GeometricTolerance2D();
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::MeshUtilities meshUtilities;

    return Polydim::PDETools::Mesh::PDE_Mesh_Utilities::compute_mesh_2D_geometry_data(geometryUtilities, meshUtilities, mesh);
}
// ***************************************************************************
void export_solution(const Polydim::examples::Parabolic_PCC_2D::Program_configuration &config,
                     const Gedim::MeshMatricesDAO &mesh,
                     const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                     const Gedim::Eigen_SparseArray<> &A,
                     const Polydim::examples::Parabolic_PCC_2D::Assembler::PostProcess_Data &post_process_data,
                     const unsigned int time_index,
                     const double &time_value,
                     const std::string &exportSolutionFolder,
                     const std::string &exportVtuFolder)
{
    const unsigned int Method_ID = static_cast<unsigned int>(config.MethodType());
    const unsigned int TEST_ID = static_cast<unsigned int>(config.TestType());
    const unsigned int space_order = config.MethodOrder();
    const unsigned int time_order = config.Theta() == 0.5 ? 2 : 1;

    {
        const char separator = ';';

        std::cout << "ProgramType" << separator;
        std::cout << "MethodType" << separator;
        std::cout << "space_method_order" << separator;
        std::cout << "time_method_order" << separator;
        std::cout << "Cell2Ds" << separator;
        std::cout << "Dofs" << separator;
        std::cout << "Strongs" << separator;
        std::cout << "h" << separator;
        std::cout << "time_index" << separator;
        std::cout << "time_value" << separator;
        std::cout << "errorL2" << separator;
        std::cout << "errorH1" << separator;
        std::cout << "normL2" << separator;
        std::cout << "normH1" << separator;
        std::cout << "nnzA" << separator;
        std::cout << "residual" << std::endl;

        std::cout.precision(2);
        std::cout << std::scientific << TEST_ID << separator;
        std::cout << std::scientific << Method_ID << separator;
        std::cout << std::scientific << space_order << separator;
        std::cout << std::scientific << time_order << separator;
        std::cout << std::scientific << mesh.Cell2DTotalNumber() << separator;
        std::cout << std::scientific << dofs_data.NumberDOFs << separator;
        std::cout << std::scientific << dofs_data.NumberStrongs << separator;
        std::cout << std::scientific << post_process_data.mesh_size << separator;
        std::cout << std::scientific << time_index << separator;
        std::cout << std::scientific << time_value << separator;
        std::cout << std::scientific << post_process_data.error_L2 << separator;
        std::cout << std::scientific << post_process_data.error_H1 << separator;
        std::cout << std::scientific << post_process_data.norm_L2 << separator;
        std::cout << std::scientific << post_process_data.norm_H1 << separator;
        std::cout << std::scientific << A.NonZeros() << separator;
        std::cout << std::scientific << post_process_data.residual_norm << std::endl;
    }

    {
        const char separator = ';';
        const std::string errorFileName = exportSolutionFolder + "/Errors" + "_" + std::to_string(TEST_ID) + "_" +
                                          std::to_string(Method_ID) + "_" + std::to_string(space_order) + "_" +
                                          std::to_string(time_order) + ".csv";
        const bool errorFileExists = Gedim::Output::FileExists(errorFileName);

        std::ofstream errorFile(errorFileName, std::ios_base::app | std::ios_base::out);
        if (!errorFileExists)
        {
            errorFile << "ProgramType" << separator;
            errorFile << "MethodType" << separator;
            errorFile << "space_method_order" << separator;
            errorFile << "time_method_order" << separator;
            errorFile << "Cell2Ds" << separator;
            errorFile << "Dofs" << separator;
            errorFile << "Strongs" << separator;
            errorFile << "h" << separator;
            errorFile << "time_index" << separator;
            errorFile << "time_value" << separator;
            errorFile << "errorL2" << separator;
            errorFile << "errorH1" << separator;
            errorFile << "normL2" << separator;
            errorFile << "normH1" << separator;
            errorFile << "nnzA" << separator;
            errorFile << "residual" << std::endl;
        }

        errorFile.precision(16);
        errorFile << std::scientific << TEST_ID << separator;
        errorFile << std::scientific << Method_ID << separator;
        errorFile << std::scientific << space_order << separator;
        errorFile << std::scientific << time_order << separator;
        errorFile << std::scientific << mesh.Cell2DTotalNumber() << separator;
        errorFile << std::scientific << dofs_data.NumberDOFs << separator;
        errorFile << std::scientific << dofs_data.NumberStrongs << separator;
        errorFile << std::scientific << post_process_data.mesh_size << separator;
        errorFile << std::scientific << time_index << separator;
        errorFile << std::scientific << time_value << separator;
        errorFile << std::scientific << post_process_data.error_L2 << separator;
        errorFile << std::scientific << post_process_data.error_H1 << separator;
        errorFile << std::scientific << post_process_data.norm_L2 << separator;
        errorFile << std::scientific << post_process_data.norm_H1 << separator;
        errorFile << std::scientific << A.NonZeros() << separator;
        errorFile << std::scientific << post_process_data.residual_norm << std::endl;

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

            exporter.Export(exportVtuFolder + "/Solution" + "_" + std::to_string(TEST_ID) + "_" +
                            std::to_string(Method_ID) + "_" + std::to_string(space_order) + "_" +
                            std::to_string(time_order) + "_" + std::to_string(time_index) + ".vtu");
        }
    }
}
// ***************************************************************************
void export_dofs(const Polydim::examples::Parabolic_PCC_2D::Program_configuration &config,
                 const Gedim::MeshMatricesDAO &mesh,
                 const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                 const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                 const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                 const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                 const Polydim::examples::Parabolic_PCC_2D::Assembler::Parabolic_PCC_2D_Problem_Data &assembler_data,
                 const Polydim::examples::Parabolic_PCC_2D::Assembler::PostProcess_Data &post_process_data,
                 const unsigned int time_index,
                 const double &time_value,
                 const std::string &exportVtuFolder)
{
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = config.GeometricTolerance1D();
    geometryUtilitiesConfig.Tolerance2D = config.GeometricTolerance2D();
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    const unsigned int space_order = config.MethodOrder();
    const unsigned int time_order = config.Theta() == 0.5 ? 2 : 1;

    const unsigned int Method_ID = static_cast<unsigned int>(config.MethodType());
    const unsigned int TEST_ID = static_cast<unsigned int>(config.TestType());

    const std::string file_path = exportVtuFolder + "/dofs" + "_" + std::to_string(TEST_ID) + "_" +
                                  std::to_string(Method_ID) + "_" + std::to_string(space_order) + "_" +
                                  std::to_string(time_order) + "_" + std::to_string(time_index) + ".vtu";

    PDETools::LocalSpace_PCC_2D::export_dofs(geometryUtilities,
                                             mesh,
                                             mesh_geometric_data,
                                             mesh_dofs_info,
                                             dofs_data,
                                             assembler_data.rightHandSide,
                                             assembler_data.solution,
                                             assembler_data.solutionDirichlet,
                                             file_path);
}
// ***************************************************************************
void export_performance(const Polydim::examples::Parabolic_PCC_2D::Program_configuration &config,
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

    case PDETools::LocalSpace_PCC_2D::MethodTypes::FEM_PCC:
        break;
    case PDETools::LocalSpace_PCC_2D::MethodTypes::VEM_PCC:
    case PDETools::LocalSpace_PCC_2D::MethodTypes::VEM_PCC_Inertia:
    case PDETools::LocalSpace_PCC_2D::MethodTypes::VEM_PCC_Ortho: {
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
            const auto &cell2D_performance = performance_data.Cell2DsPerformance[v].performance_data;

            exporter << std::scientific << v << separator;
            exporter << std::scientific << cell2D_performance.NumBoundaryQuadraturePoints << separator;
            exporter << std::scientific << cell2D_performance.NumInternalQuadraturePoints << separator;
            exporter << std::scientific << cell2D_performance.vem_analysis_data.PiNablaConditioning << separator;
            exporter << std::scientific << cell2D_performance.vem_analysis_data.Pi0kConditioning << separator;
            exporter << std::scientific << cell2D_performance.vem_analysis_data.Pi0km1Conditioning << separator;
            exporter << std::scientific << cell2D_performance.vem_analysis_data.ErrorPiNabla << separator;
            exporter << std::scientific << cell2D_performance.vem_analysis_data.ErrorPi0k << separator;
            exporter << std::scientific << cell2D_performance.vem_analysis_data.ErrorPi0km1 << separator;
            exporter << std::scientific << cell2D_performance.vem_analysis_data.ErrorHCD << separator;
            exporter << std::scientific << cell2D_performance.vem_analysis_data.ErrorGBD << separator;
            exporter << std::scientific << cell2D_performance.vem_analysis_data.ErrorStabilization << std::endl;
        }
    }
    break;
    case PDETools::LocalSpace_PCC_2D::MethodTypes::ZFEM_PCC:
        break;
    }

    exporter.close();
}
// ***************************************************************************
std::vector<double> create_time_steps(const Polydim::examples::Parabolic_PCC_2D::Program_configuration &config,
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
} // namespace program_utilities
} // namespace Parabolic_PCC_2D
} // namespace examples
} // namespace Polydim
