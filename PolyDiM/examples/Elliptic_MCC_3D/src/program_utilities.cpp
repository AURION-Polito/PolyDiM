#include "program_utilities.hpp"

#include "VTKUtilities.hpp"

namespace Polydim
{
namespace examples
{
namespace Elliptic_MCC_3D
{
namespace program_utilities
{
// ***************************************************************************
std::unique_ptr<Polydim::examples::Elliptic_MCC_3D::test::I_Test> create_test(const Polydim::examples::Elliptic_MCC_3D::Program_configuration &config)
{
    switch (config.TestType())
    {
    case Polydim::examples::Elliptic_MCC_3D::test::Test_Types::Patch_Test:
        return std::make_unique<Polydim::examples::Elliptic_MCC_3D::test::Patch_Test>();
    case Polydim::examples::Elliptic_MCC_3D::test::Test_Types::Poisson_Polynomial_Problem:
        return std::make_unique<Polydim::examples::Elliptic_MCC_3D::test::Poisson_Polynomial_Problem>();
    default:
        throw runtime_error("Test type " + to_string((unsigned int)config.TestType()) + " not supported");
    }
}
// ***************************************************************************
void create_domain_mesh(const Polydim::examples::Elliptic_MCC_3D::Program_configuration &config,
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
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_3D::Polyhedral:
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_3D::Cubic: {
        Polydim::PDETools::Mesh::PDE_Mesh_Utilities::create_mesh_3D(geometryUtilities,
                                                                    meshUtilities,
                                                                    config.MeshGenerator(),
                                                                    domain,
                                                                    config.MeshMaxVolume(),
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
Gedim::MeshUtilities::MeshGeometricData3D create_domain_mesh_geometric_properties(const Polydim::examples::Elliptic_MCC_3D::Program_configuration &config,
                                                                                  Gedim::MeshMatricesDAO &mesh)
{
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = config.GeometricTolerance1D();
    geometryUtilitiesConfig.Tolerance2D = config.GeometricTolerance2D();
    geometryUtilitiesConfig.Tolerance3D = config.GeometricTolerance3D();
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    return Polydim::PDETools::Mesh::PDE_Mesh_Utilities::compute_mesh_3D_geometry_data(geometryUtilities, mesh);
}
// ***************************************************************************
void export_solution(const Polydim::examples::Elliptic_MCC_3D::Program_configuration &config,
                     const Gedim::MeshMatricesDAO &mesh,
                     const std::vector<Polydim::PDETools::DOFs::DOFsManager::DOFsData> &dofs_data,
                     const Polydim::examples::Elliptic_MCC_3D::Assembler::Elliptic_MCC_3D_Problem_Data &assembler_data,
                     const Polydim::examples::Elliptic_MCC_3D::Assembler::PostProcess_Data &post_process_data,
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
        std::cout << "Dofs" << separator;
        std::cout << "Strongs" << separator;
        std::cout << "h" << separator;
        std::cout << "errorL2Velocity" << separator;
        std::cout << "errorL2Pressure" << separator;
        std::cout << "superErrorL2Pressure" << separator;
        std::cout << "normL2Velocity" << separator;
        std::cout << "normL2Pressure" << separator;
        std::cout << "nnzA" << separator;
        std::cout << "residual" << endl;

        std::cout.precision(2);
        std::cout << scientific << TEST_ID << separator;
        std::cout << scientific << VEM_ID << separator;
        std::cout << scientific << config.VemOrder() << separator;
        std::cout << scientific << mesh.Cell3DTotalNumber() << separator;
        std::cout << scientific << dofs_data[0].NumberDOFs + dofs_data[1].NumberDOFs << separator;
        std::cout << scientific << dofs_data[0].NumberStrongs + dofs_data[1].NumberStrongs << separator;
        std::cout << scientific << post_process_data.mesh_size << separator;
        std::cout << scientific << post_process_data.error_L2_velocity << separator;
        std::cout << scientific << post_process_data.error_L2_pressure << separator;
        std::cout << scientific << post_process_data.super_error_L2_pressure << separator;
        std::cout << scientific << post_process_data.norm_L2_velocity << separator;
        std::cout << scientific << post_process_data.norm_L2_pressure << separator;
        std::cout << scientific << assembler_data.globalMatrixA.NonZeros() << separator;
        std::cout << scientific << post_process_data.residual_norm << endl;
    }

    {
        const char separator = ';';
        const string errorFileName = exportSolutionFolder + "/Errors_" + to_string(TEST_ID) + "_" + to_string(VEM_ID) +
                                     +"_" + to_string(config.VemOrder()) + ".csv";
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
            errorFile << "errorL2Velocity" << separator;
            errorFile << "errorL2Pressure" << separator;
            errorFile << "superErrorL2Pressure" << separator;
            errorFile << "normL2Velocity" << separator;
            errorFile << "normL2Pressure" << separator;
            errorFile << "nnzA" << separator;
            errorFile << "residual" << endl;
        }

        errorFile.precision(16);
        errorFile << scientific << TEST_ID << separator;
        errorFile << scientific << VEM_ID << separator;
        errorFile << scientific << config.VemOrder() << separator;
        errorFile << scientific << mesh.Cell3DTotalNumber() << separator;
        errorFile << scientific << dofs_data[0].NumberDOFs + dofs_data[1].NumberDOFs << separator;
        errorFile << scientific << dofs_data[0].NumberStrongs + dofs_data[1].NumberStrongs << separator;
        errorFile << scientific << post_process_data.mesh_size << separator;
        errorFile << scientific << post_process_data.error_L2_velocity << separator;
        errorFile << scientific << post_process_data.error_L2_pressure << separator;
        errorFile << scientific << post_process_data.super_error_L2_pressure << separator;
        errorFile << scientific << post_process_data.norm_L2_velocity << separator;
        errorFile << scientific << post_process_data.norm_L2_pressure << separator;
        errorFile << scientific << assembler_data.globalMatrixA.NonZeros() << separator;
        errorFile << scientific << post_process_data.residual_norm << endl;

        errorFile.close();
    }

    {
        {
            Gedim::VTKUtilities exporter;
            exporter.AddPolyhedrons(mesh.Cell0DsCoordinates(),
                                    mesh.Cell3DsFacesVertices(),
                                    {{"Numeric Mean Value Pressure",
                                      Gedim::VTPProperty::Formats::Cells,
                                      static_cast<unsigned int>(post_process_data.cell3Ds_numeric_pressure.size()),
                                      post_process_data.cell3Ds_numeric_pressure.data()},
                                     {"Exact Mean Value Pressure",
                                      Gedim::VTPProperty::Formats::Cells,
                                      static_cast<unsigned int>(post_process_data.cell3Ds_exact_pressure.size()),
                                      post_process_data.cell3Ds_exact_pressure.data()},
                                     {"Pressure - ErrorL2",
                                      Gedim::VTPProperty::Formats::Cells,
                                      static_cast<unsigned int>(post_process_data.cell3Ds_error_L2_pressure.size()),
                                      post_process_data.cell3Ds_error_L2_pressure.data()},
                                     {"Velocity - ErrorL2",
                                      Gedim::VTPProperty::Formats::Cells,
                                      static_cast<unsigned int>(post_process_data.cell3Ds_error_L2_velocity.size()),
                                      post_process_data.cell3Ds_error_L2_velocity.data()}});

            exporter.Export(exportVtuFolder + "/Solution_" + to_string(TEST_ID) + "_" + to_string(VEM_ID) + +"_" +
                            to_string(config.VemOrder()) + ".vtu");
        }
    }
}
// ***************************************************************************
void export_velocity_dofs(const Polydim::examples::Elliptic_MCC_3D::Program_configuration &config,
                          const Gedim::MeshMatricesDAO &mesh,
                          const Gedim::MeshUtilities::MeshGeometricData3D &mesh_geometric_data,
                          const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                          const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                          const Polydim::examples::Elliptic_MCC_3D::Assembler::Elliptic_MCC_3D_Problem_Data &assembler_data,
                          const Polydim::examples::Elliptic_MCC_3D::Assembler::PostProcess_Data &post_process_data,
                          const std::string &exportVtuFolder)
{
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = config.GeometricTolerance1D();
    geometryUtilitiesConfig.Tolerance2D = config.GeometricTolerance2D();
    geometryUtilitiesConfig.Tolerance3D = config.GeometricTolerance3D();
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

    for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); ++c)
    {
        const auto &boundary_info = mesh_dofs_info.CellsBoundaryInfo.at(2).at(c);

        const auto &local_dofs = dofs_data.CellsDOFs[2].at(c);

        const unsigned int num_loc_dofs = local_dofs.size();

        if (num_loc_dofs == 0)
            continue;

        const std::vector<double> local_polygon_coordinates =
            geometryUtilities.EquispaceCoordinates(num_loc_dofs + 1, 0.0, 1.0, true);

        unsigned int neigh = mesh.Cell2DNeighbourCell3D(c, 0);
        if (neigh == std::numeric_limits<unsigned int>::max())
            neigh = mesh.Cell2DNeighbourCell3D(c, 1);

        const auto local_index_face = mesh.Cell3DFindFace(neigh, c);

        const Eigen::Vector3d polygon_centroid = mesh_geometric_data.Cell3DsFaces2DCentroids.at(neigh)[local_index_face];
        const auto polygonCentroidEdgesDistance = geometryUtilities.PolygonCentroidEdgesDistance(
            mesh_geometric_data.Cell3DsFaces2DVertices.at(neigh)[local_index_face],
            polygon_centroid,
            mesh_geometric_data.Cell3DsFacesEdge2DNormals.at(neigh)[local_index_face]);

        const double circle_diameter = 0.5 * geometryUtilities.PolygonInRadius(polygonCentroidEdgesDistance);

        for (unsigned int loc_i = 0; loc_i < num_loc_dofs; ++loc_i)
        {
            const auto &local_dof = local_dofs.at(loc_i);

            dof_cell_index_values.push_back(c);
            dof_dimension_values.push_back(2);
            dof_boundary_type_values.push_back(static_cast<double>(boundary_info.Type));
            dof_boundary_marker_values.push_back(boundary_info.Marker);
            if (num_loc_dofs > 1)
                dofs_coordinate.push_back(geometryUtilities.RotatePointsFrom2DTo3D(
                    polygon_centroid + circle_diameter * Eigen::Vector3d(cos(2.0 * M_PI * local_polygon_coordinates.at(loc_i)),
                                                                         sin(2.0 * M_PI * local_polygon_coordinates.at(loc_i)),
                                                                         0.0),
                    mesh_geometric_data.Cell3DsFacesRotationMatrices.at(neigh)[local_index_face],
                    mesh_geometric_data.Cell3DsFacesTranslations.at(neigh)[local_index_face]));
            else
                dofs_coordinate.push_back(geometryUtilities.RotatePointsFrom2DTo3D(
                    polygon_centroid,
                    mesh_geometric_data.Cell3DsFacesRotationMatrices.at(neigh)[local_index_face],
                    mesh_geometric_data.Cell3DsFacesTranslations.at(neigh)[local_index_face]));

            dof_type_values.push_back(static_cast<double>(local_dof.Type));
            dof_global_index_values.push_back(local_dof.Global_Index);

            switch (local_dof.Type)
            {
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
                solution_values.push_back(assembler_data.solutionNeumann.GetValue(local_dof.Global_Index));
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

    for (unsigned int c = 0; c < mesh.Cell3DTotalNumber(); ++c)
    {
        const auto &boundary_info = mesh_dofs_info.CellsBoundaryInfo.at(3).at(c);

        const auto &local_dofs = dofs_data.CellsDOFs[3].at(c);

        const unsigned int num_loc_dofs = local_dofs.size();

        const auto local_polyhedron_coordinates = geometryUtilities.fibonacci_sphere(num_loc_dofs);
        const Eigen::Vector3d polyhedron_centroid = mesh_geometric_data.Cell3DsCentroids.at(c);
        const auto polyhedron_centroid_faces_distance =
            geometryUtilities.PolyhedronCentroidFacesDistance(polyhedron_centroid,
                                                              mesh_geometric_data.Cell3DsFacesNormals.at(c),
                                                              mesh_geometric_data.Cell3DsFaces3DVertices.at(c));
        const double polyhedron_in_radius = geometryUtilities.PolyhedronInRadius(polyhedron_centroid_faces_distance);

        const double sphere_diameter = 0.5 * polyhedron_in_radius;

        for (unsigned int loc_i = 0; loc_i < num_loc_dofs; ++loc_i)
        {
            const auto &local_dof = local_dofs.at(loc_i);

            dof_cell_index_values.push_back(c);
            dof_dimension_values.push_back(3);
            dof_boundary_type_values.push_back(static_cast<double>(boundary_info.Type));
            dof_boundary_marker_values.push_back(boundary_info.Marker);
            dofs_coordinate.push_back(polyhedron_centroid + sphere_diameter * local_polyhedron_coordinates.col(loc_i));

            dof_type_values.push_back(static_cast<double>(local_dof.Type));
            dof_global_index_values.push_back(local_dof.Global_Index);

            switch (local_dof.Type)
            {
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
                solution_values.push_back(assembler_data.solutionNeumann.GetValue(local_dof.Global_Index));
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

        const unsigned int VEM_ID = static_cast<unsigned int>(config.VemType());
        const unsigned int TEST_ID = static_cast<unsigned int>(config.TestType());
        exporter.Export(exportVtuFolder + "/dofs_" + to_string(TEST_ID) + "_" + to_string(VEM_ID) + +"_" +
                        to_string(config.VemOrder()) + ".vtu");
    }
}
// ***************************************************************************
} // namespace program_utilities
} // namespace Elliptic_MCC_3D
} // namespace examples
} // namespace Polydim
