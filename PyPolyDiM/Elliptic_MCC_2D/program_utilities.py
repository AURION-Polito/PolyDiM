from Elliptic_MCC_2D.test_definition import ITest, EllipticPolynomialProblem
from pypolydim import gedim, polydim
from pypolydim.assembler_utilities import assembler_utilities
from Elliptic_MCC_2D.assembler import Assembler
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

    print("{:<5} {:<5} {:<7} {:<5} {:<10} {:<10} {:<10} {:<10} {:<20} {:<20} {:<20} {:<20} {:<20} {:<20}"
          .format('Test', 'Mesh', 'Method', 'Order', 'Cell2Ds',
                  'DOFs', 'Strongs', 'h', 'errorL2Pressure',  'errorL2Velocity',
                  'superErrorL2Pressure', 'normL2Pressure', 'normL2Velocity', 'residual'))

    print(
        "{:<5d} {:<5d} {:<7d} {:<5d} {:<10d} {:<10d} {:<10d} {:<10.2e} {:<20.2e} {:<20.2e} {:<20.2e} {:<20.2e} {:<20.2e} {:<20.2e}"
        .format(test_id, mesh_type, method_id, method_order, mesh.cell2_d_total_number(),
                count_do_fs_data.num_total_do_fs, count_do_fs_data.num_total_strongs,
                post_process_data.mesh_size, post_process_data.error_l2_pressure, post_process_data.error_l2_velocity, post_process_data.super_error_l2_pressure,
                post_process_data.norm_l2_pressure, post_process_data.norm_l2_velocity, post_process_data.residual_norm))

    file_name = file_path + "/Errors_" + str(test_id) + "_" + str(method_id) + "_" + str(method_order) + ".csv"
    if not os.path.exists(file_name):
        with open(file_name, 'w') as fd:
            fd.write('Test' + file_separator + 'Mesh' + file_separator + 'Method' + file_separator + 'Order' + file_separator + 'Cell2Ds' +
                     file_separator + 'DOFs' + file_separator + 'Strongs' + file_separator + 'h' + file_separator
             + 'errorL2Pressure' + file_separator + 'errorL2Velocity' + file_separator + 'superErrorL2Pressure' + file_separator + 'normL2Pressure' + file_separator + 'normL2Velocity' + file_separator + 'residual\n')

    with open(file_name, 'a') as fd:
        fd.write("{:<d}{:<s}{:<d}{:<s}{:<d}{:<s}{:<d}{:<s}{:<d}{:<s}{:<d}{:<s}{:<d}{:<s}{:<.16e}{:<s}{:<.16e}{:<s}{:<.16e}{:<s}{:<.16e}{:<s}{:<.16e}{:<s}{:<.16e}{:<s}{:<.16e}\n"
        .format(test_id, file_separator, mesh_type, file_separator, method_id, file_separator, method_order, file_separator,
                mesh.cell2_d_total_number(), file_separator,
                count_do_fs_data.num_total_do_fs, file_separator, count_do_fs_data.num_total_strongs, file_separator,
                post_process_data.mesh_size, file_separator, post_process_data.error_l2_pressure, file_separator, post_process_data.error_l2_velocity, file_separator,
                post_process_data.super_error_l2_pressure, file_separator,
                post_process_data.norm_l2_pressure, file_separator, post_process_data.norm_l2_velocity, file_separator, post_process_data.residual_norm))

