import argparse
import os.path

from pypolydim import polydim, gedim
from Elliptic_PCC_3D.program_utilities import create_test, create_mesh
from Elliptic_PCC_3D.assembler import Assembler
from pypolydim.vtk_utilities import VTKUtilities
from pypolydim.assembler_utilities import assembler_utilities
import cProfile


def main():

    parser =argparse.ArgumentParser()
    parser.add_argument('-order','--method-order',dest='method_order', default=1, type=int, help="Method order")
    parser.add_argument('-method','--method-type',dest='method_type', default=1, type=int, help="Method type")
    parser.add_argument('-test', '--test-id', dest='test_id', default=1, type=int, help="Test type")
    parser.add_argument('-mesh', '--mesh-type', dest='mesh_type', default=0, type=int, help="Mesh type")
    parser.add_argument('-volume', '--mesh-max-relative-volume', dest='max_relative_volume', default=0.1, type=float, help="Mesh max relative volume")
    parser.add_argument('-export', '--export-path', dest='export_path', default='./Export/Elliptic_PCC_3D', type=str, help="Export Path")
    parser.add_argument('-import', '--import-path', dest='import_path', default='./', type=str, help="Mesh Import Path")
    args = parser.parse_args()

    export_file_path = args.export_path
    if not os.path.exists(export_file_path):
        os.makedirs(export_file_path)

    # Mesh file path
    export_mesh_path = args.export_path + "/Mesh"
    if not os.path.exists(export_mesh_path):
        os.makedirs(export_mesh_path)

    mesh_type = polydim.pde_tools.mesh.pde_mesh_utilities.MeshGenerator_Types_3D(args.mesh_type)
    method_type = polydim.pde_tools.local_space_pcc_3_d.MethodTypes(args.method_type)
    method_order = args.method_order

    print("Set problem...")
    test = create_test(args.test_id)
    pde_domain = test.domain()
    boundary_info = test.boundary_info()

    print("Create mesh...")
    geometry_utilities_config = gedim.GeometryUtilitiesConfig()
    geometry_utilities = gedim.GeometryUtilities(geometry_utilities_config)
    mesh_utilities = gedim.MeshUtilities()
    vtk_utilities = VTKUtilities()

    mesh_data = gedim.MeshMatrices()
    mesh = gedim.MeshMatricesDAO(mesh_data)

    create_mesh(geometry_utilities, mesh_utilities, mesh_type, args.max_relative_area, args.import_path, pde_domain, mesh)

    print("Export Mesh...")
    vtk_utilities.export_mesh(export_mesh_path, mesh)

    print("Compute Geometric Properties...")
    mesh_geometric_data = polydim.pde_tools.mesh.pde_mesh_utilities.compute_mesh_2_d_geometry_data(geometry_utilities, mesh_utilities, mesh)

    print("Create Discrete Local Space...")
    reference_element_data = polydim.pde_tools.local_space_pcc_2_d.create_reference_element(method_type, method_order)

    mesh_connectivity_data = polydim.pde_tools.mesh.MeshMatricesDAO_mesh_connectivity_data(mesh)

    dof_manager = polydim.pde_tools.do_fs.DOFsManager()
    mesh_do_fs_info = polydim.pde_tools.local_space_pcc_2_d.set_mesh_do_fs_info(reference_element_data, mesh, boundary_info)
    do_fs_data = dof_manager.create_do_fs_2_d(mesh_do_fs_info, mesh_connectivity_data)
    assembler_utilities_obj = assembler_utilities()
    count_do_fs_data = assembler_utilities_obj.count_do_fs([do_fs_data])


    print('\x1b[6;30;42m' + "Created discrete space with ", do_fs_data.number_do_fs, " DOFs and ", do_fs_data.number_strongs, " STRONGs" + '\x1b[0m')
    print("Assemble...")
    assembler = Assembler()
    assembler_data = assembler.assemble(mesh,
                                        mesh_geometric_data,
                                        mesh_do_fs_info,
                                        do_fs_data,
                                        reference_element_data,
                                        test)

    print("Solve...")
    assembler.solve(do_fs_data, assembler_data)

    print("Compute Errors...")
    post_process_data = assembler.post_process_solution(mesh,
                                                        mesh_geometric_data,
                                                        do_fs_data,
                                                        count_do_fs_data,
                                                        reference_element_data,
                                                        assembler_data,
                                                        test)

    print()
    print("L2 error: ", post_process_data.error_l2, " H1 error: ", post_process_data.error_h1)

    print("Export Solution...")
    vtk_utilities.export_solution(export_file_path + '/solution', mesh, post_process_data.cell0_ds_numeric, cell0_d_exact_solution=post_process_data.cell0_ds_exact)

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
