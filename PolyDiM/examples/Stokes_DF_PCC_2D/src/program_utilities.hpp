#ifndef __program_utilities_H
#define __program_utilities_H

#include "DOFsManager.hpp"
#include "VTKUtilities.hpp"
#include "assembler.hpp"
#include "program_configuration.hpp"
#include "test_definition.hpp"

#include <typeindex>
#include <unordered_map>

namespace Polydim
{
namespace examples
{
namespace Stokes_DF_PCC_2D
{
namespace program_utilities
{
// ***************************************************************************
unique_ptr<Polydim::examples::Stokes_DF_PCC_2D::test::I_Test> create_test(const Polydim::examples::Stokes_DF_PCC_2D::Program_configuration &config)
{
    switch (config.TestType())
    {
    case Polydim::examples::Stokes_DF_PCC_2D::test::Test_Types::Patch_Test:
        return make_unique<Polydim::examples::Stokes_DF_PCC_2D::test::Patch_Test>();
    default:
        throw runtime_error("Test type " + to_string((unsigned int)config.TestType()) + " not supported");
    }
}
// ***************************************************************************
void create_domain_mesh(const Polydim::examples::Stokes_DF_PCC_2D::Program_configuration &config,
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
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::Polygonal: {
        Polydim::PDETools::Mesh::PDE_Mesh_Utilities::create_mesh_2D(geometryUtilities,
                                                                    meshUtilities,
                                                                    config.MeshGenerator(),
                                                                    domain,
                                                                    config.MeshMaxArea(),
                                                                    mesh);
    }
    break;
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::OFFImporter: {
        Polydim::PDETools::Mesh::PDE_Mesh_Utilities::import_mesh_2D(geometryUtilities,
                                                                    meshUtilities,
                                                                    config.MeshGenerator(),
                                                                    config.MeshImportFilePath(),
                                                                    mesh);
    }
    break;
    default:
        throw runtime_error("MeshGenerator " + to_string((unsigned int)config.MeshGenerator()) + " not supported");
    }
}
// ***************************************************************************
Gedim::MeshUtilities::MeshGeometricData2D create_domain_mesh_geometric_properties(const Polydim::examples::Stokes_DF_PCC_2D::Program_configuration &config,
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
void export_solution(const Polydim::examples::Stokes_DF_PCC_2D::Program_configuration &config,
                     const Gedim::MeshMatricesDAO &mesh,
                     const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
                     const Polydim::examples::Stokes_DF_PCC_2D::Assembler::Stokes_DF_PCC_2D_Problem_Data &assembler_data,
                     const Polydim::examples::Stokes_DF_PCC_2D::Assembler::PostProcess_Data &post_process_data,
                     const std::string &exportSolutionFolder,
                     const std::string &exportVtuFolder)
{
    const unsigned int VEM_ID = static_cast<unsigned int>(config.VemType());
    const unsigned int TEST_ID = static_cast<unsigned int>(config.TestType());

    {
        const char separator = ';';
        std::cout << "ProgramType" << separator;
        std::cout << "VemType" << separator;
        cout << "VemOrder" << separator;
        cout << "Cell2Ds" << separator;
        cout << "Dofs" << separator;
        cout << "Strongs" << separator;
        cout << "h" << separator;
        cout << "errorH1Velocity" << separator;
        cout << "errorL2Pressure" << separator;
        cout << "normH1Velocity" << separator;
        cout << "normL2Pressure" << separator;
        cout << "nnzA" << separator;
        cout << "residual" << endl;

        cout.precision(2);
        std::cout << scientific << TEST_ID << separator;
        std::cout << scientific << VEM_ID << separator;
        cout << scientific << config.VemOrder() << separator;
        cout << scientific << mesh.Cell2DTotalNumber() << separator;
        std::cout << scientific
                  << dofs_data[0].NumberDOFs + dofs_data[1].NumberDOFs + dofs_data[2].NumberDOFs + dofs_data[3].NumberDOFs
                  << separator;
        std::cout << scientific
                  << dofs_data[0].NumberStrongs + dofs_data[1].NumberStrongs + dofs_data[2].NumberStrongs + dofs_data[3].NumberStrongs
                  << separator;
        cout << scientific << post_process_data.mesh_size << separator;
        cout << scientific << post_process_data.error_H1_velocity << separator;
        cout << scientific << post_process_data.error_L2_pressure << separator;
        cout << scientific << post_process_data.norm_H1_velocity << separator;
        cout << scientific << post_process_data.norm_L2_pressure << separator;
        cout << scientific << assembler_data.globalMatrixA.NonZeros() << separator;
        cout << scientific << post_process_data.residual_norm << endl;
    }

    {
        const char separator = ';';
        const string errorFileName = exportSolutionFolder + "/Errors.csv";
        const bool errorFileExists = Gedim::Output::FileExists(errorFileName);

        std::ofstream errorFile(errorFileName, std::ios_base::app | std::ios_base::out);

        if (!errorFileExists)
        {
            errorFile << "ProgramType" << separator;
            errorFile << "VemType" << separator;
            errorFile << "VemOrder" << separator;
            errorFile << "Cell2Ds" << separator;
            errorFile << "Dofs" << separator;
            errorFile << "Strongs" << separator;
            errorFile << "h" << separator;
            errorFile << "errorH1Velocity" << separator;
            errorFile << "errorL2Pressure" << separator;
            errorFile << "normH1Velocity" << separator;
            errorFile << "normL2Pressure" << separator;
            errorFile << "nnzA" << separator;
            errorFile << "residual" << endl;
        }

        errorFile.precision(16);
        errorFile << scientific << TEST_ID << separator;
        errorFile << scientific << VEM_ID << separator;
        errorFile << scientific << config.VemOrder() << separator;
        errorFile << scientific << mesh.Cell2DTotalNumber() << separator;
        errorFile << scientific
                  << dofs_data[0].NumberDOFs + dofs_data[1].NumberDOFs + dofs_data[2].NumberDOFs + dofs_data[3].NumberDOFs
                  << separator;
        errorFile << scientific
                  << dofs_data[0].NumberStrongs + dofs_data[1].NumberStrongs + dofs_data[2].NumberStrongs + dofs_data[3].NumberStrongs
                  << separator;
        errorFile << scientific << post_process_data.mesh_size << separator;
        errorFile << scientific << post_process_data.error_H1_velocity << separator;
        errorFile << scientific << post_process_data.error_L2_pressure << separator;
        errorFile << scientific << post_process_data.norm_H1_velocity << separator;
        errorFile << scientific << post_process_data.norm_L2_pressure << separator;
        errorFile << scientific << assembler_data.globalMatrixA.NonZeros() << separator;
        errorFile << scientific << post_process_data.residual_norm << endl;

        errorFile.close();

        errorFile.close();
    }

    {
        //        {
        //            Gedim::VTKUtilities exporter;
        //            exporter.AddPolygons(mesh.Cell0DsCoordinates(),
        //                                 mesh.Cell2DsVertices(),
        //                                 {
        //                                     {
        //                                         "Numeric",
        //                                         Gedim::VTPProperty::Formats::Points,
        //                                         static_cast<unsigned
        //                                         int>(post_process_data.cell0Ds_numeric.size()),
        //                                         post_process_data.cell0Ds_numeric.data()
        //                                     },
        //                                     {
        //                                         "Exact",
        //                                         Gedim::VTPProperty::Formats::Points,
        //                                         static_cast<unsigned
        //                                         int>(post_process_data.cell0Ds_exact.size()),
        //                                         post_process_data.cell0Ds_exact.data()
        //                                     },
        //                                     {
        //                                         "ErrorL2",
        //                                         Gedim::VTPProperty::Formats::Cells,
        //                                         static_cast<unsigned
        //                                         int>(post_process_data.cell2Ds_error_L2.size()),
        //                                         post_process_data.cell2Ds_error_L2.data()
        //                                     },
        //                                     {
        //                                         "ErrorH1",
        //                                         Gedim::VTPProperty::Formats::Cells,
        //                                         static_cast<unsigned
        //                                         int>(post_process_data.cell2Ds_error_H1.size()),
        //                                         post_process_data.cell2Ds_error_H1.data()
        //                                     }
        //                                 });

        //            exporter.Export(exportVtuFolder + "/Solution_" +
        //            to_string(TEST_ID) + "_" + to_string(VEM_ID) + +
        //            "_" + to_string(config.VemOrder()) + ".vtu");
        //        }
    }
}
// ***************************************************************************
} // namespace program_utilities
} // namespace Stokes_DF_PCC_2D
} // namespace examples
} // namespace Polydim

#endif
