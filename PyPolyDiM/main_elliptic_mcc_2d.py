import argparse
import os.path

from pypolydim import polydim, gedim
from Elliptic_MCC_2D.program_utilities import create_test, create_mesh, export_errors
from Elliptic_MCC_2D.assembler import Assembler
from pypolydim.export_vtk_utilities import ExportVTKUtilities
from pypolydim.assembler_utilities import assembler_utilities
import cProfile
from typing import List


def main():

    parser =argparse.ArgumentParser()
    parser.add_argument('-order','--method-order',dest='method_order', default=0, type=int, help="Method order: 0, 1, 2, ...")
    parser.add_argument('-method','--method-type',dest='method_type', default=1, type=int, help="Method type: 1 - VEM, 3 - VEM Ortho")
    parser.add_argument('-test', '--test-id', dest='test_id', default=1, type=int, help="Test type: 1 - Elliptic")
    parser.add_argument('-mesh', '--mesh-type', dest='mesh_type', default=0, type=int,
                        help="Mesh type: 0 - Triangular; 1 - Minimal; 2 - Polygonal; 3 - OFF Importer; 4 - CsvImporter (; separator); 5 - Squared ")
    parser.add_argument('-tol1', '--tolerance-1-d', dest='tolerance1_d', default=1.0e-12, type=float, help="Geometric Tolerance 1D")
    parser.add_argument('-tol2', '--tolerance-2-d', dest='tolerance2_d', default=1.0e-14, type=float, help="Geometric Tolerance 2D")
    parser.add_argument('-area', '--mesh-max-relative-area', dest='max_relative_area', default=0.1, type=float, help="Mesh max relative area")
    parser.add_argument('-export', '--export-path', dest='export_path', default='./Export/Elliptic_MCC_2D', type=str, help="Export Path")
    parser.add_argument('-import', '--import-path', dest='import_path', default='./', type=str, help="Mesh Import Path")
    args = parser.parse_args()

    export_file_path = args.export_path
    if not os.path.exists(export_file_path):
        os.makedirs(export_file_path)

    # Mesh file path
    export_mesh_path = args.export_path + "/Mesh"
    if not os.path.exists(export_mesh_path):
        os.makedirs(export_mesh_path)

    mesh_type = polydim.pde_tools.mesh.pde_mesh_utilities.MeshGenerator_Types_2D(args.mesh_type)
    method_type = polydim.pde_tools.local_space_mcc_2_d.MethodTypes(args.method_type)
    method_order = args.method_order

    print("Set problem...")
    test = create_test(args.test_id)
    pde_domain = test.domain()
    boundary_info = test.boundary_info()

    print("Create mesh...")
    geometry_utilities_config = gedim.GeometryUtilitiesConfig()
    geometry_utilities_config.tolerance1_d = args.tolerance1_d
    geometry_utilities_config.tolerance2_d = args.tolerance2_d
    geometry_utilities = gedim.GeometryUtilities(geometry_utilities_config)
    mesh_utilities = gedim.MeshUtilities()
    vtk_utilities = ExportVTKUtilities()

    mesh_data = gedim.MeshMatrices()
    mesh = gedim.MeshMatricesDAO(mesh_data)

    create_mesh(geometry_utilities, mesh_utilities, mesh_type, args.max_relative_area, args.import_path, pde_domain, mesh)

    print("Export Mesh...")
    vtk_utilities.export_mesh(export_mesh_path, mesh)

    print("Compute Geometric Properties...")
    mesh_geometric_data = polydim.pde_tools.mesh.pde_mesh_utilities.compute_mesh_2_d_geometry_data(geometry_utilities, mesh_utilities, mesh)

    print("Create Discrete Local Space...")
    reference_element_data = polydim.pde_tools.local_space_mcc_2_d.create_reference_element(method_type, method_order)
    reference_element_num_do_fs = polydim.pde_tools.local_space_mcc_2_d.reference_element_num_do_fs(reference_element_data)
    mesh_connectivity_data = polydim.pde_tools.mesh.MeshMatricesDAO_mesh_connectivity_data(mesh)

    dof_manager = polydim.pde_tools.do_fs.DOFsManager()

    mesh_do_fs_info: List[polydim.pde_tools.do_fs.DOFsManager.MeshDOFsInfo] = []
    do_fs_data: List[polydim.pde_tools.do_fs.DOFsManager.DOFsData] = []
    do_fs_data_indices: List[polydim.pde_tools.do_fs.DOFsManager.CellsDOFsIndicesData] = []

    for h in range(2):
        constant_do_fs_info = polydim.pde_tools.do_fs.DOFsManager.ConstantDOFsInfo()
        constant_do_fs_info.num_do_fs = reference_element_num_do_fs[h]
        constant_do_fs_info.boundary_info = boundary_info

        mesh_do_fs_info.append(dof_manager.create_constant_do_fs_info_2_d(mesh_connectivity_data, constant_do_fs_info))
        do_fs_data.append(dof_manager.create_do_fs_2_d(mesh_do_fs_info[h], mesh_connectivity_data))
        do_fs_data_indices.append(dof_manager.compute_cells_do_fs_indices(do_fs_data[h], 2))

    assembler_utilities_obj = assembler_utilities()
    count_do_fs_data = assembler_utilities_obj.count_do_fs(do_fs_data)


    print('\x1b[6;30;42m' + "Created discrete space with ", count_do_fs_data.num_total_do_fs, " DOFs and ", count_do_fs_data.num_total_strongs, " STRONGs" + '\x1b[0m')
    print("Assemble...")
    assembler = Assembler()
    assembler_data = assembler.assemble(geometry_utilities_config,
                                        mesh,
                                        mesh_geometric_data,
                                        mesh_do_fs_info,
                                        do_fs_data,
                                        do_fs_data_indices,
                                        count_do_fs_data,
                                        reference_element_data,
                                        test)

    print("Solve...")
    assembler.solve(count_do_fs_data, assembler_data)

    print("Compute Errors...")
    post_process_data = assembler.post_process_solution(geometry_utilities_config,
                                                        mesh,
                                                        mesh_geometric_data,
                                                        do_fs_data,
                                                        do_fs_data_indices,
                                                        count_do_fs_data,
                                                        reference_element_data,
                                                        assembler_data,
                                                        test)

    print("Export Solution...")

    export_errors(export_file_path, args.test_id, args.mesh_type, args.method_type, args.method_order, mesh, count_do_fs_data, post_process_data)

    point_data = {}
    cell_data = {"Numeric Mean Pressure": post_process_data.cell2_ds_numeric_pressure,
                 "Exact Mean Pressure": post_process_data.cell2_ds_exact_pressure,
                 "Error L2 Pressure": post_process_data.cell2_ds_error_l2_pressure,
                 "Error L2 Velocity": post_process_data.cell2_ds_error_l2_velocity}

    vtk_utilities.export_cells2_d(export_file_path + '/Solution_' + str(args.test_id) + '_' + str(args.method_type)
                                    + '_' + str(method_order), mesh, point_data=point_data, cell_data=cell_data)

    print('\x1b[6;30;42m' + "Finish" + '\x1b[0m')

if __name__=='__main__':

    export_file_path = "./Export"
    if not os.path.exists(export_file_path):
        os.makedirs(export_file_path)

    pr = cProfile.Profile()
    pr.enable()
    main()
    pr.disable()
    pr.dump_stats("./Export/program.prof")  # call "snakeviz program.proof" from command line to visualize the result
