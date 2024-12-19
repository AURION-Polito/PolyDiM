#ifndef __program_utilities_H
#define __program_utilities_H

#include "assembler.hpp"
#include "program_configuration.hpp"
#include "DOFsManager.hpp"
#include "VTKUtilities.hpp"
#include "VEM_MCC_2D_Velocity_LocalSpace.hpp"
#include "VEM_MCC_2D_Partial_Velocity_LocalSpace.hpp"
#include "VEM_MCC_2D_Ortho_Velocity_LocalSpace.hpp"
#include "test_definition.hpp"

#include <unordered_map>
#include <typeindex>

namespace Polydim
{
namespace examples
{
namespace Elliptic_MCC_2D
{
namespace program_utilities
{
// ***************************************************************************
// Centralized mapping function for IDs
unsigned int VemType(const std::type_index& type)
{
    static const std::unordered_map<std::type_index, unsigned int> typeToID = {
        {typeid(Polydim::VEM::MCC::VEM_MCC_2D_Velocity_LocalSpace), 1},
        {typeid(Polydim::VEM::MCC::VEM_MCC_2D_Partial_Velocity_LocalSpace), 2},
        {typeid(Polydim::VEM::MCC::VEM_MCC_2D_Ortho_Velocity_LocalSpace), 3}
    };

    auto it = typeToID.find(type);
    if (it != typeToID.end()) {
        return it->second;
    } else {
        throw std::runtime_error("Class type not recognized.");
    }
}
// Helper template to get the ID for a specific class
template <typename T>
unsigned int VemType() {
    return VemType(typeid(T));
}
// ***************************************************************************
void create_domain_mesh(const Polydim::examples::Elliptic_MCC_2D::Program_configuration& config,
                        const Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D& domain,
                        Gedim::MeshMatricesDAO& mesh)
{
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = 1.0e-8;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::MeshUtilities meshUtilities;

    switch (config.MeshGenerator())
    {
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::Triangular:
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::Minimal:
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::Polygonal:
    {
        Polydim::PDETools::Mesh::PDE_Mesh_Utilities::create_mesh_2D(geometryUtilities,
                                                                    meshUtilities,
                                                                    config.MeshGenerator(),
                                                                    domain,
                                                                    config.MeshMaxArea(),
                                                                    mesh);
    }
    break;
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::OFFImporter:
    {
        Polydim::PDETools::Mesh::PDE_Mesh_Utilities::import_mesh_2D(geometryUtilities,
                                                                    meshUtilities,
                                                                    config.MeshGenerator(),
                                                                    config.MeshImportFilePath(),
                                                                    mesh);
    }
    break;
    default:
        throw runtime_error("MeshGenerator " +
                            to_string((unsigned int)config.MeshGenerator()) +
                            " not supported");
    }
}
// ***************************************************************************
Gedim::MeshUtilities::MeshGeometricData2D create_domain_mesh_geometric_properties(const Polydim::examples::Elliptic_MCC_2D::Program_configuration& config,
                                                                                  const Gedim::MeshMatricesDAO& mesh)
{
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = config.GeometricTolerance1D();
    geometryUtilitiesConfig.Tolerance2D = config.GeometricTolerance2D();
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::MeshUtilities meshUtilities;

    return Polydim::PDETools::Mesh::PDE_Mesh_Utilities::compute_mesh_2D_geometry_data(geometryUtilities,
                                                                                      meshUtilities,
                                                                                      mesh);
}
// ***************************************************************************
template<typename VEM_LOCAL_SPACE_TYPE, typename TEST_TYPE>
void export_solution(const Polydim::examples::Elliptic_MCC_2D::Program_configuration& config,
                     const Gedim::MeshMatricesDAO& mesh,
                     const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData>& dofs_data,
                     const typename Polydim::examples::Elliptic_MCC_2D::Assembler<VEM_LOCAL_SPACE_TYPE>::Elliptic_MCC_2D_Problem_Data& assembler_data,
                     const typename Polydim::examples::Elliptic_MCC_2D::Assembler<VEM_LOCAL_SPACE_TYPE>::PostProcess_Data& post_process_data,
                     const std::string& exportSolutionFolder,
                     const std::string& exportVtuFolder)
{

    const unsigned int VEM_ID = Polydim::examples::Elliptic_MCC_2D::program_utilities::VemType<VEM_LOCAL_SPACE_TYPE>();
    const unsigned int TEST_ID = Polydim::examples::Elliptic_MCC_2D::test::TestType<TEST_TYPE>();

    {
        const char separator = ';';
        std::cout << "ProgramType" << separator;
        std::cout << "VemType" << separator;
        std::cout << "VemOrder" << separator;
        std::cout << "Cell2Ds" <<  separator;
        std::cout << "Dofs" <<  separator;
        std::cout << "Strongs" <<  separator;
        std::cout << "h" <<  separator;
        std::cout << "errorL2Velocity" <<  separator;
        std::cout << "errorL2Pressure" << separator;
        std::cout << "normL2Velocity" <<  separator;
        std::cout << "normL2Pressure" << separator;
        std::cout << "nnzA" << separator;
        std::cout << "residual" << endl;

        std::cout.precision(2);
        std::cout << scientific<< TEST_ID << separator;
        std::cout << scientific<< VEM_ID << separator;
        std::cout << scientific << mesh.Cell2DTotalNumber()<< separator;
        std::cout << scientific << dofs_data[0].NumberDOFs + dofs_data[1].NumberDOFs << separator;
        std::cout << scientific << dofs_data[0].NumberStrongs + dofs_data[1].NumberStrongs << separator;
        std::cout << scientific << post_process_data.mesh_size << separator;
        std::cout << scientific << post_process_data.error_L2_velocity << separator;
        std::cout << scientific << post_process_data.error_L2_pressure << separator;
        std::cout << scientific << post_process_data.norm_L2_velocity << separator;
        std::cout << scientific << post_process_data.norm_L2_pressure << separator;
        std::cout << scientific << assembler_data.globalMatrixA.NonZeros() << separator;
        std::cout << scientific << post_process_data.residual_norm << endl;

    }

    {    
        const char separator = ';';
        const string errorFileName = exportSolutionFolder +
                                     "/Errors.csv";
        const bool errorFileExists = Gedim::Output::FileExists(errorFileName);

        std::ofstream errorFile(errorFileName,
                                std::ios_base::app | std::ios_base::out);
        if (!errorFileExists)
        {
            errorFile << "ProgramType" << separator;
            errorFile << "VemType" << separator;
            errorFile << "VemOrder" << separator;
            errorFile << "Cell2Ds" <<  separator;
            errorFile << "Dofs" <<  separator;
            errorFile << "Strongs" <<  separator;
            errorFile << "h" <<  separator;
            errorFile << "errorL2Velocity" <<  separator;
            errorFile << "errorL2Pressure" << separator;
            errorFile << "normL2Velocity" <<  separator;
            errorFile << "normL2Pressure" << separator;
            errorFile << "nnzA" << separator;
            errorFile << "residual" << endl;
        }

        errorFile.precision(16);
        errorFile << scientific<< TEST_ID << separator;
        errorFile << scientific<< VEM_ID << separator;
        errorFile << scientific << config.VemOrder()<< separator;
        errorFile << scientific << mesh.Cell2DTotalNumber()<< separator;
        errorFile << scientific << dofs_data[0].NumberDOFs + dofs_data[1].NumberDOFs << separator;
        errorFile << scientific << dofs_data[0].NumberStrongs + dofs_data[1].NumberStrongs << separator;
        errorFile << scientific << post_process_data.mesh_size << separator;
        errorFile << scientific << post_process_data.error_L2_velocity << separator;
        errorFile << scientific << post_process_data.error_L2_pressure << separator;
        errorFile << scientific << post_process_data.norm_L2_velocity << separator;
        errorFile << scientific << post_process_data.norm_L2_pressure << separator;
        errorFile << scientific << assembler_data.globalMatrixA.NonZeros() << separator;
        errorFile << scientific << post_process_data.residual_norm << endl;

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
//                                         static_cast<unsigned int>(post_process_data.cell0Ds_numeric.size()),
//                                         post_process_data.cell0Ds_numeric.data()
//                                     },
//                                     {
//                                         "Exact",
//                                         Gedim::VTPProperty::Formats::Points,
//                                         static_cast<unsigned int>(post_process_data.cell0Ds_exact.size()),
//                                         post_process_data.cell0Ds_exact.data()
//                                     },
//                                     {
//                                         "ErrorL2",
//                                         Gedim::VTPProperty::Formats::Cells,
//                                         static_cast<unsigned int>(post_process_data.cell2Ds_error_L2.size()),
//                                         post_process_data.cell2Ds_error_L2.data()
//                                     },
//                                     {
//                                         "ErrorH1",
//                                         Gedim::VTPProperty::Formats::Cells,
//                                         static_cast<unsigned int>(post_process_data.cell2Ds_error_H1.size()),
//                                         post_process_data.cell2Ds_error_H1.data()
//                                     }
//                                 });

//            exporter.Export(exportVtuFolder + "/Solution_" + to_string(TEST_ID) + "_" + to_string(VEM_ID) + + "_" + to_string(config.VemOrder()) + ".vtu");
//        }
    }
}
// ***************************************************************************
}
}
}
}

#endif
