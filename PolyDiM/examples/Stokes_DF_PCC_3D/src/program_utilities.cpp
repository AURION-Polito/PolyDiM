#include "program_utilities.hpp"

namespace Polydim
{
namespace examples
{
namespace Stokes_DF_PCC_3D
{
namespace program_utilities
{
// ***************************************************************************
unique_ptr<Polydim::examples::Stokes_DF_PCC_3D::test::I_Test> create_test(const Polydim::examples::Stokes_DF_PCC_3D::Program_configuration &config)
{
    switch (config.TestType())
    {
    case Polydim::examples::Stokes_DF_PCC_3D::test::Test_Types::Patch_Test:
        return make_unique<Polydim::examples::Stokes_DF_PCC_3D::test::Patch_Test>();
    default:
        throw runtime_error("Test type " + to_string((unsigned int)config.TestType()) + " not supported");
    }
}
// ***************************************************************************
void create_domain_mesh(const Polydim::examples::Stokes_DF_PCC_3D::Program_configuration &config,
                        const Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_3D &domain,
                        Gedim::MeshMatricesDAO &mesh)
{
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = config.GeometricTolerance1D();
    geometryUtilitiesConfig.Tolerance2D = config.GeometricTolerance2D();
    geometryUtilitiesConfig.Tolerance3D = config.GeometricTolerance3D();
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::MeshUtilities meshUtilities;

    switch (config.MeshGenerator())
    {
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_3D::Tetrahedral:
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_3D::Minimal:
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_3D::Polyhedral: {
        Polydim::PDETools::Mesh::PDE_Mesh_Utilities::create_mesh_3D(geometryUtilities,
                                                                    meshUtilities,
                                                                    config.MeshGenerator(),
                                                                    domain,
                                                                    config.MeshMaxArea(),
                                                                    mesh);
    }
    break;
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_3D::OVMImporter:
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_3D::VtkImporter:
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_3D::CsvImporter: {
        Polydim::PDETools::Mesh::PDE_Mesh_Utilities::import_mesh_3D(geometryUtilities,
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
Gedim::MeshUtilities::MeshGeometricData3D create_domain_mesh_geometric_properties(const Polydim::examples::Stokes_DF_PCC_3D::Program_configuration &config,
                                                                                  const Gedim::MeshMatricesDAO &mesh)
{
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = config.GeometricTolerance1D();
    geometryUtilitiesConfig.Tolerance2D = config.GeometricTolerance2D();
    geometryUtilitiesConfig.Tolerance3D = config.GeometricTolerance3D();
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::MeshUtilities meshUtilities;

    return Polydim::PDETools::Mesh::PDE_Mesh_Utilities::compute_mesh_3D_geometry_data(geometryUtilities, meshUtilities, mesh);
}
// ***************************************************************************
void export_solution(const Polydim::examples::Stokes_DF_PCC_3D::Program_configuration &config,
                     const Gedim::MeshMatricesDAO &mesh,
                     const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
                     const Polydim::examples::Stokes_DF_PCC_3D::Assembler::Stokes_DF_PCC_3D_Problem_Data &assembler_data,
                     const Polydim::examples::Stokes_DF_PCC_3D::Assembler::PostProcess_Data &post_process_data,
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
        cout << "Cell3Ds" << separator;
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
        cout << scientific << mesh.Cell3DTotalNumber() << separator;
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
            errorFile << "Cell3Ds" << separator;
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
        errorFile << scientific << mesh.Cell3DTotalNumber() << separator;
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
    }

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
                                  {"Numeric Velocity - Z",
                                   Gedim::VTPProperty::Formats::Points,
                                   static_cast<unsigned int>(post_process_data.cell0Ds_numeric_velocity[2].size()),
                                   post_process_data.cell0Ds_numeric_velocity[2].data()},
                                  {"Exact Velocity - X",
                                   Gedim::VTPProperty::Formats::Points,
                                   static_cast<unsigned int>(post_process_data.cell0Ds_exact_velocity[0].size()),
                                   post_process_data.cell0Ds_exact_velocity[0].data()},
                                  {"Exact Velocity - Y",
                                   Gedim::VTPProperty::Formats::Points,
                                   static_cast<unsigned int>(post_process_data.cell0Ds_exact_velocity[1].size()),
                                   post_process_data.cell0Ds_exact_velocity[1].data()},
                                  {"Exact Velocity - Z",
                                   Gedim::VTPProperty::Formats::Points,
                                   static_cast<unsigned int>(post_process_data.cell0Ds_exact_velocity[2].size()),
                                   post_process_data.cell0Ds_exact_velocity[2].data()},
                                  {"ErrorL2Pressure",
                                   Gedim::VTPProperty::Formats::Cells,
                                   static_cast<unsigned int>(post_process_data.cell2Ds_error_L2_pressure.size()),
                                   post_process_data.cell2Ds_error_L2_pressure.data()},
                                  {"ErrorH1Velocity",
                                   Gedim::VTPProperty::Formats::Cells,
                                   static_cast<unsigned int>(post_process_data.cell2Ds_error_H1_velocity.size()),
                                   post_process_data.cell2Ds_error_H1_velocity.data()}});

            exporter.Export(exportVtuFolder + "/Solution_" + to_string(TEST_ID) + "_" + to_string(VEM_ID) + +"_" +
                            to_string(config.VemOrder()) + ".vtu");
        }
    }
}
// ***************************************************************************
void export_velocity_dofs(const Polydim::examples::Stokes_DF_PCC_3D::Program_configuration &config,
                          const Gedim::MeshMatricesDAO &mesh,
                          const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                          const std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> &mesh_dofs_info,
                          const VEM::DF_PCC::VEM_DF_PCC_3D_Velocity_ReferenceElement_Data &vem_velocity_reference_element_data,
                          const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
                          const Polydim::examples::Stokes_DF_PCC_3D::Assembler::Stokes_DF_PCC_3D_Problem_Data &assembler_data,
                          const Polydim::examples::Stokes_DF_PCC_3D::Assembler::PostProcess_Data &post_process_data,
                          const std::string &exportVtuFolder)
{
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = config.GeometricTolerance1D();
    geometryUtilitiesConfig.Tolerance2D = config.GeometricTolerance2D();
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    //    std::list<Eigen::Vector3d> dofs_coordinate;
    //    std::list<double> solution_values;
    //    std::list<double> rhs_values;
    //    std::list<double> dof_global_index_values;
    //    std::list<double> dof_type_values;
    //    std::list<double> dof_cell_index_values;
    //    std::list<double> dof_dimension_values;
    //    std::list<double> dof_boundary_type_values;
    //    std::list<double> dof_boundary_marker_values;

    //    for (unsigned int c = 0; c < mesh.Cell0DTotalNumber(); ++c)
    //    {
    //        const auto &boundary_info = mesh_dofs_info.CellsBoundaryInfo.at(0).at(c);

    //        const auto &local_dofs = dofs_data.CellsDOFs[0].at(c);

    //        for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
    //        {
    //            const auto &local_dof = local_dofs.at(loc_i);

    //            dof_cell_index_values.push_back(c);
    //            dof_dimension_values.push_back(0);
    //            dof_boundary_type_values.push_back(static_cast<double>(boundary_info.Type));
    //            dof_boundary_marker_values.push_back(boundary_info.Marker);
    //            dofs_coordinate.push_back(mesh.Cell0DCoordinates(c));
    //            dof_type_values.push_back(static_cast<double>(local_dof.Type));
    //            dof_global_index_values.push_back(local_dof.Global_Index);

    //            switch (local_dof.Type)
    //            {
    //            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
    //                solution_values.push_back(assembler_data.solutionDirichlet.GetValue(local_dof.Global_Index));
    //                rhs_values.push_back(std::nan(""));
    //                break;
    //            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
    //                solution_values.push_back(assembler_data.solution.GetValue(local_dof.Global_Index));
    //                rhs_values.push_back(assembler_data.rightHandSide.GetValue(local_dof.Global_Index));
    //                break;
    //            default:
    //                throw std::runtime_error("Unknown DOF Type");
    //            }
    //        }
    //    }

    //    for (unsigned int c = 0; c < mesh.Cell1DTotalNumber(); ++c)
    //    {
    //        const auto &boundary_info = mesh_dofs_info.CellsBoundaryInfo.at(1).at(c);

    //        const auto &local_dofs = dofs_data.CellsDOFs[1].at(c);

    //        const unsigned int num_loc_dofs = local_dofs.size();
    //        const Eigen::VectorXd local_edge_coordinates =
    //            vem_reference_element_data.Quadrature.ReferenceEdgeDOFsInternalPoints.row(0);
    //        const Eigen::Vector3d edge_origin = mesh.Cell1DOriginCoordinates(c);
    //        const Eigen::Vector3d edge_tangent = mesh.Cell1DEndCoordinates(c) - edge_origin;

    //        for (unsigned int loc_i = 0; loc_i < num_loc_dofs; ++loc_i)
    //        {
    //            const auto &local_dof = local_dofs.at(loc_i);

    //            dof_cell_index_values.push_back(c);
    //            dof_dimension_values.push_back(1);
    //            dof_boundary_type_values.push_back(static_cast<double>(boundary_info.Type));
    //            dof_boundary_marker_values.push_back(boundary_info.Marker);
    //            dofs_coordinate.push_back(edge_origin + local_edge_coordinates(loc_i) * edge_tangent);
    //            dof_type_values.push_back(static_cast<double>(local_dof.Type));
    //            dof_global_index_values.push_back(local_dof.Global_Index);

    //            switch (local_dof.Type)
    //            {
    //            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
    //                solution_values.push_back(assembler_data.solutionDirichlet.GetValue(local_dof.Global_Index));
    //                rhs_values.push_back(0.0);
    //                break;
    //            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
    //                solution_values.push_back(assembler_data.solution.GetValue(local_dof.Global_Index));
    //                rhs_values.push_back(assembler_data.rightHandSide.GetValue(local_dof.Global_Index));
    //                break;
    //            default:
    //                throw std::runtime_error("Unknown DOF Type");
    //            }
    //        }
    //    }

    //    for (unsigned int c = 0; c < mesh.Cell3DTotalNumber(); ++c)
    //    {
    //        const auto &boundary_info = mesh_dofs_info.CellsBoundaryInfo.at(2).at(c);

    //        const auto &local_dofs = dofs_data.CellsDOFs[2].at(c);

    //        const unsigned int num_loc_dofs = local_dofs.size();

    //        const auto local_polygon_coordinates = geometryUtilities.EquispaceCoordinates(num_loc_dofs + 1, 0.0, 1.0,
    //        true); const Eigen::Vector3d polygon_centroid = mesh_geometric_data.Cell2DsCentroids.at(c); const auto
    //        polygonCentroidEdgesDistance =
    //            geometryUtilities.PolygonCentroidEdgesDistance(mesh_geometric_data.Cell2DsVertices.at(c),
    //                                                           mesh_geometric_data.Cell2DsCentroids.at(c),
    //                                                           mesh_geometric_data.Cell2DsEdgeNormals.at(c));
    //        const double circle_diameter = 0.5 * geometryUtilities.PolygonInRadius(polygonCentroidEdgesDistance);

    //        for (unsigned int loc_i = 0; loc_i < num_loc_dofs; ++loc_i)
    //        {
    //            const auto &local_dof = local_dofs.at(loc_i);

    //            dof_cell_index_values.push_back(c);
    //            dof_dimension_values.push_back(2);
    //            dof_boundary_type_values.push_back(static_cast<double>(boundary_info.Type));
    //            dof_boundary_marker_values.push_back(boundary_info.Marker);
    //            dofs_coordinate.push_back(
    //                polygon_centroid + circle_diameter * Eigen::Vector3d(cos(2.0 * M_PI *
    //                local_polygon_coordinates.at(loc_i)),
    //                                                                     sin(2.0 * M_PI *
    //                                                                     local_polygon_coordinates.at(loc_i)), 0.0));
    //            dof_type_values.push_back(static_cast<double>(local_dof.Type));
    //            dof_global_index_values.push_back(local_dof.Global_Index);

    //            switch (local_dof.Type)
    //            {
    //            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
    //                solution_values.push_back(assembler_data.solutionDirichlet.GetValue(local_dof.Global_Index));
    //                rhs_values.push_back(0.0);
    //                break;
    //            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
    //                solution_values.push_back(assembler_data.solution.GetValue(local_dof.Global_Index));
    //                rhs_values.push_back(assembler_data.rightHandSide.GetValue(local_dof.Global_Index));
    //                break;
    //            default:
    //                throw std::runtime_error("Unknown DOF Type");
    //            }
    //        }
    //    }

    //    {
    //        Eigen::MatrixXd coordinates(3, dofs_coordinate.size());
    //        unsigned int c = 0;
    //        for (const auto &dof_coordinate : dofs_coordinate)
    //            coordinates.col(c++) << dof_coordinate;
    //        const auto rhs_values_data = std::vector<double>(rhs_values.begin(), rhs_values.end());
    //        const auto solution_values_data = std::vector<double>(solution_values.begin(), solution_values.end());
    //        const auto dof_global_index_values_data =
    //            std::vector<double>(dof_global_index_values.begin(), dof_global_index_values.end());
    //        const auto dof_type_values_data = std::vector<double>(dof_type_values.begin(), dof_type_values.end());
    //        const auto dof_cell_index_values_data = std::vector<double>(dof_cell_index_values.begin(),
    //        dof_cell_index_values.end()); const auto dof_dimension_values_data =
    //        std::vector<double>(dof_dimension_values.begin(), dof_dimension_values.end()); const auto
    //        dof_boundary_type_values_data =
    //            std::vector<double>(dof_boundary_type_values.begin(), dof_boundary_type_values.end());
    //        const auto dof_boundary_marker_values_data =
    //            std::vector<double>(dof_boundary_marker_values.begin(), dof_boundary_marker_values.end());

    //        Gedim::VTKUtilities exporter;
    //        exporter.AddPoints(coordinates,
    //                           {{"cell_dimension",
    //                             Gedim::VTPProperty::Formats::Points,
    //                             static_cast<unsigned int>(dof_dimension_values_data.size()),
    //                             dof_dimension_values_data.data()},
    //                            {"cell_index",
    //                             Gedim::VTPProperty::Formats::Points,
    //                             static_cast<unsigned int>(dof_cell_index_values_data.size()),
    //                             dof_cell_index_values_data.data()},
    //                            {"boundary_type",
    //                             Gedim::VTPProperty::Formats::Points,
    //                             static_cast<unsigned int>(dof_boundary_type_values_data.size()),
    //                             dof_boundary_type_values_data.data()},
    //                            {"boundary_marker",
    //                             Gedim::VTPProperty::Formats::Points,
    //                             static_cast<unsigned int>(dof_boundary_marker_values_data.size()),
    //                             dof_boundary_marker_values_data.data()},
    //                            {"dof_global_index",
    //                             Gedim::VTPProperty::Formats::Points,
    //                             static_cast<unsigned int>(dof_global_index_values_data.size()),
    //                             dof_global_index_values_data.data()},
    //                            {"dof_type",
    //                             Gedim::VTPProperty::Formats::Points,
    //                             static_cast<unsigned int>(dof_type_values_data.size()),
    //                             dof_type_values_data.data()},
    //                            {"rhs",
    //                             Gedim::VTPProperty::Formats::Points,
    //                             static_cast<unsigned int>(rhs_values_data.size()),
    //                             rhs_values_data.data()},
    //                            {"solution",
    //                             Gedim::VTPProperty::Formats::Points,
    //                             static_cast<unsigned int>(solution_values_data.size()),
    //                             solution_values_data.data()}});

    //        exporter.Export(exportVtuFolder + "/dofs.vtu");
    //    }
}
// ***************************************************************************
void export_discrepancy_errors(const Polydim::examples::Stokes_DF_PCC_3D::Program_configuration &config,
                               const Gedim::MeshMatricesDAO &mesh,
                               const Polydim::examples::Stokes_DF_PCC_3D::Assembler::DiscrepancyErrors_Data &discrepancy_errors_data,
                               const std::string &exportSolutionFolder,
                               const std::string &exportVtuFolder)
{
    const unsigned int VEM_ID = static_cast<unsigned int>(config.VemType());
    const unsigned int TEST_ID = static_cast<unsigned int>(config.TestType());

    {
        const char separator = ';';
        std::cout << "ProgramType" << separator;
        std::cout << "VemType" << separator;
        std::cout << "VemOrder" << separator;
        std::cout << "Cell3Ds" << separator;
        std::cout << "VelocityDofsRatio" << separator;
        std::cout << "PressureDofsRatio" << separator;
        std::cout << "errorH1Velocity" << separator;
        std::cout << "errorL2Pressure" << separator;
        std::cout << "normH1FULLVelocity" << separator;
        std::cout << "normL2FULLPressure" << separator;
        std::cout << "reducedResidual" << separator;
        std::cout << "fullResidual" << endl;

        std::cout.precision(2);
        std::cout << scientific << TEST_ID << separator;
        std::cout << scientific << VEM_ID << separator;
        std::cout << scientific << config.VemOrder() << separator;
        std::cout << scientific << mesh.Cell3DTotalNumber() << separator;
        std::cout << scientific << discrepancy_errors_data.velocity_dofs_ratio << separator;
        std::cout << scientific << discrepancy_errors_data.pressure_dofs_ratio << separator;
        std::cout << scientific << discrepancy_errors_data.discrepancy_error_H1_velocity << separator;
        std::cout << scientific << discrepancy_errors_data.discrepancy_error_L2_pressure << separator;
        std::cout << scientific << discrepancy_errors_data.full_norm_H1_velocity << separator;
        std::cout << scientific << discrepancy_errors_data.full_norm_L2_pressure << separator;
        std::cout << scientific << discrepancy_errors_data.reduced_residual_norm << separator;
        std::cout << scientific << discrepancy_errors_data.residual_norm << endl;
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
            errorFile << "VelocityDofsRatio" << separator;
            errorFile << "PressureDofsRatio" << separator;
            errorFile << "errorH1Velocity" << separator;
            errorFile << "errorL2Pressure" << separator;
            errorFile << "normH1FULLVelocity" << separator;
            errorFile << "normL2FULLPressure" << separator;
            errorFile << "reducedResidual" << separator;
            errorFile << "fullResidual" << endl;
        }

        errorFile.precision(16);
        errorFile << scientific << TEST_ID << separator;
        errorFile << scientific << VEM_ID << separator;
        errorFile << scientific << config.VemOrder() << separator;
        errorFile << scientific << mesh.Cell3DTotalNumber() << separator;
        errorFile << scientific << discrepancy_errors_data.velocity_dofs_ratio << separator;
        errorFile << scientific << discrepancy_errors_data.pressure_dofs_ratio << separator;
        errorFile << scientific << discrepancy_errors_data.discrepancy_error_H1_velocity << separator;
        errorFile << scientific << discrepancy_errors_data.discrepancy_error_L2_pressure << separator;
        errorFile << scientific << discrepancy_errors_data.full_norm_H1_velocity << separator;
        errorFile << scientific << discrepancy_errors_data.full_norm_L2_pressure << separator;
        errorFile << scientific << discrepancy_errors_data.reduced_residual_norm << separator;
        errorFile << scientific << discrepancy_errors_data.residual_norm << endl;

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

            exporter.Export(exportVtuFolder + "/DiscrepancyErrors_" + to_string(TEST_ID) + "_" + to_string(VEM_ID) +
                            +"_" + to_string(config.VemOrder()) + ".vtu");
        }
    }
}
// ***************************************************************************
} // namespace program_utilities
} // namespace Stokes_DF_PCC_3D
} // namespace examples
} // namespace Polydim
