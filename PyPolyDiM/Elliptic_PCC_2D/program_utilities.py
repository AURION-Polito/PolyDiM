from Elliptic_PCC_2D.test_definition import ITest, EllipticPolynomialProblem
from pypolydim import gedim, polydim
from pypolydim.assembler_utilities import assembler_utilities
from Elliptic_PCC_2D.assembler import Assembler
import os

def create_test(test_id: int) -> ITest:

    match test_id:
        case 1:
            return EllipticPolynomialProblem()
        case _:
            raise ValueError("not valid test id")


def create_mesh(geometry_utilities: gedim.GeometryUtilities, mesh_utilities: gedim.MeshUtilities,
                mesh_type: polydim.pde_tools.mesh.pde_mesh_utilities.MeshGenerator_Types_2D,
                mesh_max_relative_area: float, import_path: str,
                pde_domain: polydim.pde_tools.mesh.pde_mesh_utilities.PDE_Domain_2D,
                mesh: gedim.MeshMatricesDAO):

    if (polydim.pde_tools.mesh.pde_mesh_utilities.MeshGenerator_Types_2D.triangular == mesh_type or
            polydim.pde_tools.mesh.pde_mesh_utilities.MeshGenerator_Types_2D.minimal == mesh_type or
            polydim.pde_tools.mesh.pde_mesh_utilities.MeshGenerator_Types_2D.polygonal == mesh_type or
            polydim.pde_tools.mesh.pde_mesh_utilities.MeshGenerator_Types_2D.squared == mesh_type or
        polydim.pde_tools.mesh.pde_mesh_utilities.MeshGenerator_Types_2D.random_distorted == mesh_type):

        polydim.pde_tools.mesh.pde_mesh_utilities.create_mesh_2_d(geometry_utilities,
                                                                  mesh_utilities,
                                                                  mesh_type,
                                                                  pde_domain,
                                                                  mesh_max_relative_area,
                                                                  mesh)

    elif (polydim.pde_tools.mesh.pde_mesh_utilities.MeshGenerator_Types_2D.csv_importer == mesh_type or
            polydim.pde_tools.mesh.pde_mesh_utilities.MeshGenerator_Types_2D.off_importer == mesh_type):

        polydim.pde_tools.mesh.pde_mesh_utilities.import_mesh_2_d(mesh_utilities,
                                                                  mesh_type,
                                                                  import_path,
                                                                  mesh)
    else:
        raise ValueError("MeshGenerator " + str(mesh_type) + " not supported")


def export_errors(file_path: str, test_id: int, mesh_type: int, method_id: int, method_order: int,
                  mesh:  gedim.MeshMatricesDAO,
                  count_do_fs_data: assembler_utilities.CountDOFsData,
                  post_process_data: Assembler.PostProcessData,
                  file_separator = ';') -> None:

    print("{:<5} {:<5} {:<7} {:<5} {:<10} {:<10} {:<10} {:<10} {:<10} {:<10} {:<10} {:<10} {:<10}"
          .format('Test', 'Mesh', 'Method', 'Order', 'Cell2Ds', 'DOFs',
                  'Strongs', 'h', 'errorL2',  'errorH1', 'normL2', 'normH1', 'residual'))

    print(
        "{:<5d} {:<5d} {:<7d} {:<5d} {:<10d} {:<10d} {:<10d} {:<10.2e} {:<10.2e} {:<10.2e} {:<10.2e} {:<10.2e} {:<10.2e}"
        .format(test_id, mesh_type, method_id, method_order, mesh.cell2_d_total_number(),
                count_do_fs_data.num_total_do_fs, count_do_fs_data.num_total_strongs,
                post_process_data.mesh_size, post_process_data.error_l2, post_process_data.error_h1,
                post_process_data.norm_l2, post_process_data.norm_h1, post_process_data.residual_norm))

    file_name = file_path + "/Errors_" + str(test_id) + "_" + str(method_id) + "_" + str(method_order) + ".csv"
    if not os.path.exists(file_name):
        with open(file_name, 'w') as fd:
            fd.write('Test' + file_separator + 'Mesh' + file_separator + 'Method' + file_separator + 'Order' + file_separator + 'Cell2Ds' +
                     file_separator + 'DOFs' + file_separator + 'Strongs' + file_separator + 'h' + file_separator
             + 'errorL2' + file_separator + 'errorH1' + file_separator + 'normL2' + file_separator + 'normH1' + file_separator + 'residual\n')

    with open(file_name, 'a') as fd:
        fd.write("{:<d}{:<s}{:<d}{:<s}{:<d}{:<s}{:<d}{:<s}{:<d}{:<s}{:<d}{:<s}{:<d}{:<s}{:<.16e}{:<s}{:<.16e}{:<s}{:<.16e}{:<s}{:<.16e}{:<s}{:<.16e}{:<s}{:<.16e}\n"
        .format(test_id, file_separator, mesh_type, file_separator, method_id, file_separator, method_order, file_separator,
                mesh.cell2_d_total_number(), file_separator,
                count_do_fs_data.num_total_do_fs, file_separator, count_do_fs_data.num_total_strongs, file_separator,
                post_process_data.mesh_size, file_separator, post_process_data.error_l2, file_separator, post_process_data.error_h1, file_separator,
                post_process_data.norm_l2, file_separator, post_process_data.norm_h1, file_separator, post_process_data.residual_norm))

