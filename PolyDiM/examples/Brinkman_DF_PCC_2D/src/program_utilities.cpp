#include "program_utilities.hpp"

#include "VTKUtilities.hpp"

namespace Polydim
{
namespace examples
{
namespace Brinkman_DF_PCC_2D
{
namespace program_utilities
{
// ***************************************************************************
std::unique_ptr<Polydim::examples::Brinkman_DF_PCC_2D::test::I_Test> create_test(const Polydim::examples::Brinkman_DF_PCC_2D::Program_configuration &config)
{
    switch (config.TestType())
    {
    case Polydim::examples::Brinkman_DF_PCC_2D::test::Test_Types::Patch_Test:
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
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::CsvImporter: {
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
void export_solution(const Polydim::examples::Brinkman_DF_PCC_2D::Program_configuration &config,
                     const Gedim::MeshMatricesDAO &mesh,
                     const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
                     const Polydim::PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                     const Polydim::examples::Brinkman_DF_PCC_2D::Assembler::Stokes_DF_PCC_2D_Problem_Data &assembler_data,
                     const Polydim::examples::Brinkman_DF_PCC_2D::Assembler::PostProcess_Data &post_process_data,
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
        std::cout << std::scientific << VEM_ID << separator;
        std::cout << std::scientific << config.VemOrder() << separator;
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
                                          std::to_string(VEM_ID) + +"_" + std::to_string(config.VemOrder()) + ".csv";
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
            errorFile << "residual" << std::endl;
        }

        errorFile.precision(16);
        errorFile << std::scientific << TEST_ID << separator;
        errorFile << std::scientific << VEM_ID << separator;
        errorFile << std::scientific << config.VemOrder() << separator;
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
                                         std::to_string(VEM_ID) + +"_" + std::to_string(config.VemOrder()) + ".csv";

        std::ofstream fluxFile(fluxFileName, std::ios_base::out);

        fluxFile << "EdgeMarker" << separator;
        fluxFile << "Flux" << std::endl;

        fluxFile.precision(16);
        for(const auto &f : post_process_data.flux)
            fluxFile << std::scientific << f.first << separator << f.second << std::endl;;

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

            exporter.Export(exportVtuFolder + "/Velocity_" + std::to_string(TEST_ID) + "_" + std::to_string(VEM_ID) +
                            +"_" + std::to_string(config.VemOrder()) + ".vtu");
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

            exporter.Export(exportVtuFolder + "/Pressure_" + std::to_string(TEST_ID) + "_" + std::to_string(VEM_ID) +
                            +"_" + std::to_string(config.VemOrder()) + ".vtu");
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

            exporter.Export(exportVtuFolder + "/Velocity_" + std::to_string(TEST_ID) + "_" + std::to_string(VEM_ID) +
                            +"_" + std::to_string(config.VemOrder()) + ".vtu");
        }

        {
            Gedim::VTKUtilities exporter;
            exporter.AddPolygons(post_process_data.repeated_vertices_coordinates,
                                 post_process_data.repeated_connectivity,
                                 {{"Numeric Pressure",
                                   Gedim::VTPProperty::Formats::Points,
                                   static_cast<unsigned int>(post_process_data.cell0Ds_numeric_pressure.size()),
                                   post_process_data.cell0Ds_numeric_pressure.data()}});

            exporter.Export(exportVtuFolder + "/Pressure_" + std::to_string(TEST_ID) + "_" + std::to_string(VEM_ID) +
                            +"_" + std::to_string(config.VemOrder()) + ".vtu");
        }
    }
}
// ***************************************************************************
void export_velocity_dofs(const Polydim::examples::Brinkman_DF_PCC_2D::Program_configuration &config,
                          const Gedim::MeshMatricesDAO &mesh,
                          const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                          const std::vector<Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo> &mesh_dofs_info,
                          const VEM::DF_PCC::VEM_DF_PCC_2D_Velocity_ReferenceElement_Data &vem_velocity_reference_element_data,
                          const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
                          const Polydim::PDETools::Assembler_Utilities::count_dofs_data &count_dofs,
                          const Polydim::examples::Brinkman_DF_PCC_2D::Assembler::Stokes_DF_PCC_2D_Problem_Data &assembler_data,
                          const Polydim::examples::Brinkman_DF_PCC_2D::Assembler::PostProcess_Data &post_process_data,
                          const std::string &exportVtuFolder)
{
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = config.GeometricTolerance1D();
    geometryUtilitiesConfig.Tolerance2D = config.GeometricTolerance2D();
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    std::list<Eigen::Vector3d> dofs_coordinate;
    std::list<double> dof_dimension_values;
    std::list<double> dof_cell_index_values;
    std::list<std::array<double, 2>> solution_values;
    std::list<std::array<double, 2>> rhs_values;
    std::list<std::array<double, 2>> dof_global_index_values;
    std::list<std::array<double, 2>> dof_type_values;
    std::list<std::array<double, 2>> dof_boundary_type_values;
    std::list<std::array<double, 2>> dof_boundary_marker_values;

    for (unsigned int c = 0; c < mesh.Cell0DTotalNumber(); ++c)
    {
        for (unsigned int loc_i = 0; loc_i < dofs_data[0].CellsDOFs[0].at(c).size(); ++loc_i)
        {
            dofs_coordinate.push_back(mesh.Cell0DCoordinates(c));
            dof_dimension_values.push_back(0);
            dof_cell_index_values.push_back(c);

            std::array<double, 2> sol;
            std::array<double, 2> rhs;
            std::array<double, 2> dof_global;
            std::array<double, 2> dof_type;
            std::array<double, 2> dof_boundary_type;
            std::array<double, 2> dof_boundary_marker;

            for (unsigned int h = 0; h < 2; h++)
            {
                const auto &boundary_info = mesh_dofs_info[h].CellsBoundaryInfo.at(0).at(c);

                const auto &local_dofs = dofs_data[h].CellsDOFs[0].at(c);

                const auto &local_dof = local_dofs.at(loc_i);

                dof_boundary_type[h] = static_cast<double>(boundary_info.Type);
                dof_boundary_marker[h] = boundary_info.Marker;
                dof_type[h] = static_cast<double>(local_dof.Type);

                switch (local_dof.Type)
                {
                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
                    dof_global[h] = local_dof.Global_Index + count_dofs.offsets_Strongs[h];
                    sol[h] = assembler_data.solutionDirichlet.GetValue(local_dof.Global_Index + count_dofs.offsets_Strongs[h]);
                    rhs[h] = std::nan("");
                    break;
                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                    dof_global[h] = local_dof.Global_Index + count_dofs.offsets_DOFs[h];
                    sol[h] = assembler_data.solution.GetValue(local_dof.Global_Index + count_dofs.offsets_DOFs[h]);
                    rhs[h] = assembler_data.rightHandSide.GetValue(local_dof.Global_Index + count_dofs.offsets_DOFs[h]);
                    break;
                default:
                    throw std::runtime_error("Unknown DOF Type");
                }
            }

            solution_values.push_back(sol);
            rhs_values.push_back(rhs);
            dof_global_index_values.push_back(dof_global);
            dof_type_values.push_back(dof_type);
            dof_boundary_type_values.push_back(dof_boundary_type);
            dof_boundary_marker_values.push_back(dof_boundary_marker);
        }
    }

    for (unsigned int c = 0; c < mesh.Cell1DTotalNumber(); ++c)
    {
        const std::vector<double> local_edge_coordinates =
            geometryUtilities.EquispaceCoordinates(dofs_data[0].CellsDOFs[1].at(c).size(), 0.0, 1.0, false);
        const Eigen::Vector3d edge_origin = mesh.Cell1DOriginCoordinates(c);
        const Eigen::Vector3d edge_tangent = mesh.Cell1DEndCoordinates(c) - edge_origin;

        for (unsigned int loc_i = 0; loc_i < dofs_data[0].CellsDOFs[1].at(c).size(); ++loc_i)
        {
            std::array<double, 2> sol;
            std::array<double, 2> rhs;
            std::array<double, 2> dof_global;
            std::array<double, 2> dof_type;
            std::array<double, 2> dof_boundary_type;
            std::array<double, 2> dof_boundary_marker;

            dofs_coordinate.push_back(edge_origin + local_edge_coordinates[loc_i] * edge_tangent);
            dof_cell_index_values.push_back(c);
            dof_dimension_values.push_back(1);

            for (unsigned int h = 0; h < 2; h++)
            {
                const auto &boundary_info = mesh_dofs_info[h].CellsBoundaryInfo.at(1).at(c);

                const auto &local_dofs = dofs_data[h].CellsDOFs[1].at(c);

                const auto &local_dof = local_dofs.at(loc_i);

                dof_boundary_type[h] = static_cast<double>(boundary_info.Type);
                dof_boundary_marker[h] = boundary_info.Marker;
                dof_type[h] = static_cast<double>(local_dof.Type);

                switch (local_dof.Type)
                {
                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
                    dof_global[h] = local_dof.Global_Index + count_dofs.offsets_Strongs[h];
                    sol[h] = assembler_data.solutionDirichlet.GetValue(local_dof.Global_Index + count_dofs.offsets_Strongs[h]);
                    rhs[h] = std::nan("");
                    break;
                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                    dof_global[h] = local_dof.Global_Index + count_dofs.offsets_DOFs[h];
                    sol[h] = assembler_data.solution.GetValue(local_dof.Global_Index + count_dofs.offsets_DOFs[h]);
                    rhs[h] = assembler_data.rightHandSide.GetValue(local_dof.Global_Index + count_dofs.offsets_DOFs[h]);
                    break;
                default:
                    throw std::runtime_error("Unknown DOF Type");
                }
            }

            solution_values.push_back(sol);
            rhs_values.push_back(rhs);
            dof_global_index_values.push_back(dof_global);
            dof_type_values.push_back(dof_type);
            dof_boundary_type_values.push_back(dof_boundary_type);
            dof_boundary_marker_values.push_back(dof_boundary_marker);
        }
    }

    for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); ++c)
    {
        const unsigned int num_loc_dofs = dofs_data[2].CellsDOFs[2].at(c).size();

        if (num_loc_dofs == 0)
            break;

        const auto local_polygon_coordinates = geometryUtilities.EquispaceCoordinates(num_loc_dofs + 1, 0.0, 1.0, true);
        const Eigen::Vector3d polygon_centroid = mesh_geometric_data.Cell2DsCentroids.at(c);
        const auto polygonCentroidEdgesDistance =
            geometryUtilities.PolygonCentroidEdgesDistance(mesh_geometric_data.Cell2DsVertices.at(c),
                                                           mesh_geometric_data.Cell2DsCentroids.at(c),
                                                           mesh_geometric_data.Cell2DsEdgeNormals.at(c));
        const double circle_diameter = 0.5 * geometryUtilities.PolygonInRadius(polygonCentroidEdgesDistance);

        for (unsigned int loc_i = 0; loc_i < num_loc_dofs; ++loc_i)
        {
            dof_cell_index_values.push_back(c);
            dof_dimension_values.push_back(2);
            dofs_coordinate.push_back(
                polygon_centroid + circle_diameter * Eigen::Vector3d(cos(2.0 * M_PI * local_polygon_coordinates.at(loc_i)),
                                                                     sin(2.0 * M_PI * local_polygon_coordinates.at(loc_i)),
                                                                     0.0));

            std::array<double, 2> sol = {std::nan(""), std::nan("")};
            std::array<double, 2> rhs = {std::nan(""), std::nan("")};
            std::array<double, 2> dof_global = {std::nan(""), std::nan("")};
            std::array<double, 2> dof_type = {std::nan(""), std::nan("")};
            std::array<double, 2> dof_boundary_type = {std::nan(""), std::nan("")};
            std::array<double, 2> dof_boundary_marker = {std::nan(""), std::nan("")};

            for (unsigned int h = 0; h < 1; h++)
            {
                const auto &boundary_info = mesh_dofs_info[h + 2].CellsBoundaryInfo.at(2).at(c);
                const auto &local_dofs = dofs_data[h + 2].CellsDOFs[2].at(c);

                const auto &local_dof = local_dofs.at(loc_i);

                dof_boundary_type[h] = static_cast<double>(boundary_info.Type);
                dof_boundary_marker[h] = boundary_info.Marker;
                dof_type[h] = static_cast<double>(local_dof.Type);

                switch (local_dof.Type)
                {
                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
                    dof_global[h] = local_dof.Global_Index + count_dofs.offsets_Strongs[h + 2];
                    sol[h] = assembler_data.solutionDirichlet.GetValue(local_dof.Global_Index + count_dofs.offsets_Strongs[h + 2]);
                    rhs[h] = std::nan("");
                    break;
                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                    dof_global[h] = local_dof.Global_Index + count_dofs.offsets_DOFs[h + 2];
                    sol[h] = assembler_data.solution.GetValue(local_dof.Global_Index + count_dofs.offsets_DOFs[h + 2]);
                    rhs[h] = assembler_data.rightHandSide.GetValue(local_dof.Global_Index + count_dofs.offsets_DOFs[h + 2]);
                    break;
                default:
                    throw std::runtime_error("Unknown DOF Type");
                }
            }

            solution_values.push_back(sol);
            rhs_values.push_back(rhs);
            dof_global_index_values.push_back(dof_global);
            dof_type_values.push_back(dof_type);
            dof_boundary_type_values.push_back(dof_boundary_type);
            dof_boundary_marker_values.push_back(dof_boundary_marker);
        }
    }

    {
        const unsigned int n = dofs_coordinate.size();
        Eigen::MatrixXd coordinates(3, dofs_coordinate.size());
        unsigned int c = 0;
        for (const auto &dof_coordinate : dofs_coordinate)
            coordinates.col(c++) << dof_coordinate;
        const auto dof_dimension_values_data = std::vector<double>(dof_dimension_values.begin(), dof_dimension_values.end());
        const auto dof_cell_index_values_data = std::vector<double>(dof_cell_index_values.begin(), dof_cell_index_values.end());
        Eigen::VectorXd rhs_values_data(2 * n);
        c = 0;
        for (const auto &v : rhs_values)
        {
            for (unsigned int d = 0; d < 2; d++)
                rhs_values_data[2 * c + d] = v[d];
            c++;
        }
        Eigen::VectorXd solution_values_data(2 * n);
        c = 0;
        for (const auto &v : solution_values)
        {
            for (unsigned int d = 0; d < 2; d++)
                solution_values_data[2 * c + d] = v[d];
            c++;
        }
        Eigen::VectorXd dof_global_index_values_data(2 * n);
        c = 0;
        for (const auto &v : dof_global_index_values)
        {
            for (unsigned int d = 0; d < 2; d++)
                dof_global_index_values_data[2 * c + d] = v[d];
            c++;
        }
        Eigen::VectorXd dof_type_values_data(2 * n);
        c = 0;
        for (const auto &v : dof_type_values)
        {
            for (unsigned int d = 0; d < 2; d++)
                dof_type_values_data[2 * c + d] = v[d];
            c++;
        }
        Eigen::VectorXd dof_boundary_type_values_data(2 * n);
        c = 0;
        for (const auto &v : dof_boundary_type_values)
        {
            for (unsigned int d = 0; d < 2; d++)
                dof_boundary_type_values_data[2 * c + d] = v[d];
            c++;
        }
        Eigen::VectorXd dof_boundary_marker_values_data(2 * n);
        c = 0;
        for (const auto &v : dof_boundary_marker_values)
        {
            for (unsigned int d = 0; d < 2; d++)
                dof_boundary_marker_values_data[2 * c + d] = v[d];
            c++;
        }

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
                             Gedim::VTPProperty::Formats::PointsArray2,
                             static_cast<unsigned int>(dof_boundary_type_values_data.size()),
                             dof_boundary_type_values_data.data()},
                            {"boundary_marker",
                             Gedim::VTPProperty::Formats::PointsArray2,
                             static_cast<unsigned int>(dof_boundary_marker_values_data.size()),
                             dof_boundary_marker_values_data.data()},
                            {"dof_global_index",
                             Gedim::VTPProperty::Formats::PointsArray2,
                             static_cast<unsigned int>(dof_global_index_values_data.size()),
                             dof_global_index_values_data.data()},
                            {"dof_type",
                             Gedim::VTPProperty::Formats::PointsArray2,
                             static_cast<unsigned int>(dof_type_values_data.size()),
                             dof_type_values_data.data()},
                            {"rhs",
                             Gedim::VTPProperty::Formats::PointsArray2,
                             static_cast<unsigned int>(rhs_values_data.size()),
                             rhs_values_data.data()},
                            {"solution",
                             Gedim::VTPProperty::Formats::PointsArray2,
                             static_cast<unsigned int>(solution_values_data.size()),
                             solution_values_data.data()}});

        const unsigned int VEM_ID = static_cast<unsigned int>(config.VemType());
        const unsigned int TEST_ID = static_cast<unsigned int>(config.TestType());

        exporter.Export(exportVtuFolder + "/dofs_" + std::to_string(TEST_ID) + "_" + std::to_string(VEM_ID) + +"_" +
                        std::to_string(config.VemOrder()) + ".vtu");
    }
}
// ***************************************************************************
void export_discrepancy_errors(const Polydim::examples::Brinkman_DF_PCC_2D::Program_configuration &config,
                               const Gedim::MeshMatricesDAO &mesh,
                               const Polydim::examples::Brinkman_DF_PCC_2D::Assembler::DiscrepancyErrors_Data &discrepancy_errors_data,
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
        std::cout << std::scientific << VEM_ID << separator;
        std::cout << std::scientific << config.VemOrder() << separator;
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
                                          std::to_string(VEM_ID) + +"_" + std::to_string(config.VemOrder()) + ".csv";
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
            errorFile << "discrepancyErrorH1Velocity" << separator;
            errorFile << "discrepancyErrorL2Pressure" << separator;
            errorFile << "normH1FULLVelocity" << separator;
            errorFile << "normL2FULLProjectedPressure" << separator;
            errorFile << "reducedResidual" << separator;
            errorFile << "fullResidual" << std::endl;
        }

        errorFile.precision(16);
        errorFile << std::scientific << TEST_ID << separator;
        errorFile << std::scientific << VEM_ID << separator;
        errorFile << std::scientific << config.VemOrder() << separator;
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
                            std::to_string(VEM_ID) + +"_" + std::to_string(config.VemOrder()) + ".vtu");
        }
    }
}
// ***************************************************************************
} // namespace program_utilities
} // namespace Brinkman_DF_PCC_2D
} // namespace examples
} // namespace Polydim
