import os
import csv
import math
import numpy as np


def run_program(program_folder,
                program_path,
                run_folder,
                vem_type,
                vem_order,
                test_type,
                mesh_generator,
                mesh_max_volume,
                compute_discrepancy,
                solver_type,
                num_it_solver,
                tol_solver):
    export_path = os.path.join(program_folder,
                               export_folder,
                               "{0}_TT{1}".format(
                                   run_folder,
                                   test_type),
                               "{0}_TT{1}_VT{2}".format(
                                   run_folder,
                                   test_type,
                                   vem_type),
                               "{0}_TT{1}_VT{2}_VO{3}".format(
                                   run_folder,
                                   test_type,
                                   vem_type,
                                   vem_order))

    program_parameters = "VemType:uint={0}".format(vem_type)
    program_parameters += " VemOrder:uint={0}".format(vem_order)
    program_parameters += " ExportFolder:string={0}".format(export_path)
    program_parameters += " TestType:uint={0}".format(test_type)
    program_parameters += " MeshGenerator:uint={0}".format(mesh_generator)
    program_parameters += " MeshMaxVolume:double={0}".format(mesh_max_volume)
    program_parameters += " ComputeVEMPerformance:bool={0}".format(0)
    program_parameters += " ComputeDiscrepancyError:bool={0}".format(compute_discrepancy)
    program_parameters += " SolverType:uint={0}".format(solver_type)
    program_parameters += " MaxNumberIterations:uint={0}".format(num_it_solver)
    program_parameters += " RelResidualTolerance:double={0}".format(tol_solver)

    output_file = os.path.join(program_folder,
                               "terminal.log")

    run_label = "VemType {0}".format(vem_type)
    run_label += " VemOrder {0}".format(vem_order)
    run_label += " TestType {0}".format(test_type)
    run_label += " MeshGenerator {0}".format(mesh_generator)
    run_label += " MeshMaxVolume {0}".format(mesh_max_volume)
    print("Run " + run_label + "...")
    os.system(program_path + " " + program_parameters + "> " + output_file)
    os.system("mv " + output_file + " " + export_path)
    print("Run SUCCESS")

    return export_path


def import_errors(export_path, vem_type, vem_order, test_type):
    errors_file = os.path.join(export_path,
                               "Solution",
                               "Errors_" + str(test_type) + "_" + str(vem_type) + "_" + str(vem_order) + ".csv")
    errors = []
    with open(errors_file, newline='') as csvfile:
        file_reader = csv.reader(csvfile, delimiter=';')
        data = list(file_reader)

        counter = 0
        for col in data:
            errors_row = []
            if counter == 0:
                errors_row.append(col[6])
                errors_row.append(col[7])
                errors_row.append(col[8])
                errors_row.append(col[9])
                errors_row.append(col[10])
            else:
                errors_row.append(float(col[6]))
                errors_row.append(float(col[7]))
                errors_row.append(float(col[8]))
                errors_row.append(float(col[9]))
                errors_row.append(float(col[10]))
            errors.append(errors_row)
            counter += 1

    return errors


def import_discrepancy_errors(export_path, vem_type, vem_order, test_type):
    errors_file = os.path.join(export_path,
                               "Solution",
                               "DiscrepancyErrors_" + str(test_type) + "_" + str(vem_type) + "_" + str(
                                   vem_order) + ".csv")
    errors = []
    with open(errors_file, newline='') as csvfile:
        file_reader = csv.reader(csvfile, delimiter=';')
        data = list(file_reader)

        counter = 0
        for col in data:
            errors_row = []
            if counter == 0:
                errors_row.append(col[6])
                errors_row.append(col[7])
                errors_row.append(col[8])
                errors_row.append(col[9])
            else:
                errors_row.append(float(col[6]))
                errors_row.append(float(col[7]))
                errors_row.append(float(col[8]))
                errors_row.append(float(col[9]))
            errors.append(errors_row)
            counter += 1

    return errors


def test_errors(errors,
                vem_order,
                tol,
                test_type):
    num_rows = len(errors)

    if (num_rows == 2):
        print("Patch: ", abs(errors[1][1]) / abs(errors[1][3]), abs(errors[1][2]) / abs(errors[1][4]))
        assert abs(errors[1][1]) < tol * abs(errors[1][3])
        assert abs(errors[1][2]) < tol * abs(errors[1][4])
    else:
        if test_type == 3:
            errors = np.array(errors[1:])
            slope_L2_pres = np.polyfit(np.log(errors[:, 0]), np.log(errors[:, 2]), 1)[0]
            print("Num. Ref. ", str(num_rows-1), " : ", errors[:, 1], " , ", slope_L2_pres)
            assert round(slope_L2_pres) == round(float(vem_order))
            for nr in range(num_rows-1):
                assert abs(errors[nr, 1]) < tol
        elif test_type == 4:
            errors = np.array(errors[1:])
            slope_L2_pres = np.polyfit(np.log(errors[:, 0]), np.log(errors[:, 2]), 1)[0]
            slope_H1_vel = np.polyfit(np.log(errors[:, 0]), np.log(errors[:, 1]), 1)[0]
            print("Num. Ref. ", str(num_rows-1), " : ", slope_H1_vel, slope_L2_pres)
            assert round(slope_L2_pres) == round(float(vem_order))
            assert round(slope_H1_vel) >= round(float(vem_order + 2))
        else:
            errors = np.array(errors[1:])
            slope_L2_pres = np.polyfit(np.log(errors[:, 0]), np.log(errors[:, 2]), 1)[0]
            slope_H1_vel = np.polyfit(np.log(errors[:, 0]), np.log(errors[:, 1]), 1)[0]
            print("Num. Ref. ", str(num_rows - 1), " : ", slope_H1_vel, slope_L2_pres)
            assert round(slope_L2_pres) >= round(float(vem_order))
            assert round(slope_H1_vel) >= round(float(vem_order))


def test_discrepancy_errors(errors, tol):
    num_rows = len(errors)

    for r in range(num_rows - 1):
        print("Discrepancy - Id. Ref. ", str(r), " : ", abs(errors[r + 1][0]),
              abs(errors[r + 1][1]))
        assert abs(errors[r + 1][0]) < tol
        assert abs(errors[r + 1][1]) < tol


if __name__ == "__main__":
    program_folder = os.path.dirname(os.path.realpath(__file__))
    program_path = os.path.join(".", program_folder, "Stokes_DF_PCC_3D")

    remove_folder = False

    export_folder = "integration_tests"
    os.system("rm -rf " + os.path.join(program_folder, export_folder))
    tol = 5.0e-3

    print("RUN TESTS...")

    test_type = 1
    mesh_generator = 1
    mesh_max_volume = 0.0
    vem_types = [1]
    vem_orders = [2, 3, 4]
    solver_type = 0
    tol_solver = 1.0e-12
    num_it_solver = 50
    for vem_type in vem_types:
        for vem_order in vem_orders:
            export_path = run_program(program_folder,
                                      program_path,
                                      "Run_MG{0}".format(mesh_generator),
                                      vem_type,
                                      vem_order,
                                      test_type,
                                      mesh_generator,
                                      mesh_max_volume,
                                      0,
                                      solver_type,
                                      num_it_solver,
                                      tol_solver)
            errors = import_errors(export_path, vem_type, vem_order, test_type)
            test_errors(errors,
                        vem_order,
                        tol,
                        test_type)
            if remove_folder:
                os.system("rm -rf " + os.path.join(program_folder, export_path))

    vem_types = [2]
    for vem_type in vem_types:
        for vem_order in vem_orders:
            export_path = run_program(program_folder,
                                      program_path,
                                      "Run_MG{0}".format(mesh_generator),
                                      vem_type,
                                      vem_order,
                                      test_type,
                                      mesh_generator,
                                      mesh_max_volume,
                                      1,
                                      solver_type,
                                      num_it_solver,
                                      tol_solver)
            errors = import_discrepancy_errors(export_path, vem_type, vem_order, test_type)
            test_discrepancy_errors(errors, tol)
            if remove_folder:
                os.system("rm -rf " + os.path.join(program_folder, export_path))

    # test_type = 2
    # mesh_generator = 6
    # mesh_max_volumes = [0.125 * 0.125 * 0.125, 0.0005787037]
    # vem_types = [1]
    # solver_type = 1
    # tol_solver = 1.0e-12
    # num_it_solver = 50
    # for vem_type in vem_types:
    #     for vem_order in vem_orders:
    #         for mesh_max_volume in mesh_max_volumes:
    #             export_path = run_program(program_folder,
    #                                       program_path,
    #                                       "Run_MG{0}".format(mesh_generator),
    #                                       vem_type,
    #                                       vem_order,
    #                                       test_type,
    #                                       mesh_generator,
    #                                       mesh_max_volume,
    #                                       0,
    #                                       solver_type,
    #                                       num_it_solver,
    #                                       tol_solver)
    #         errors = import_errors(export_path, vem_type, vem_order, test_type)
    #         test_errors(errors,
    #                     vem_order,
    #                     tol,
    #                     test_type)
    #         if remove_folder:
    #             os.system("rm -rf " + os.path.join(program_folder, export_path))
    #
    # vem_types = [2]
    # for vem_type in vem_types:
    #     for vem_order in vem_orders:
    #         for mesh_max_volume in mesh_max_volumes:
    #             export_path = run_program(program_folder,
    #                                       program_path,
    #                                       "Run_MG{0}".format(mesh_generator),
    #                                       vem_type,
    #                                       vem_order,
    #                                       test_type,
    #                                       mesh_generator,
    #                                       mesh_max_volume,
    #                                       1,
    #                                       solver_type,
    #                                       num_it_solver,
    #                                       tol_solver)
    #         errors = import_discrepancy_errors(export_path, vem_type, vem_order, test_type)
    #         test_discrepancy_errors(errors, 1.0e-03)
    #         if remove_folder:
    #             os.system("rm -rf " + os.path.join(program_folder, export_path))

    test_type = 3
    mesh_generator = 6
    mesh_max_volumes = [0.125 * 0.125 * 0.125, 0.0005787037]
    vem_types = [1]
    vem_orders = [2, 3]
    solver_type = 1
    tol_solver = 1.0e-12
    num_it_solver = 30
    for vem_type in vem_types:
        for vem_order in vem_orders:
            for mesh_max_volume in mesh_max_volumes:
                export_path = run_program(program_folder,
                                          program_path,
                                          "Run_MG{0}".format(mesh_generator),
                                          vem_type,
                                          vem_order,
                                          test_type,
                                          mesh_generator,
                                          mesh_max_volume,
                                          0,
                                          solver_type,
                                          num_it_solver,
                                          tol_solver)
            errors = import_errors(export_path, vem_type, vem_order, test_type)
            test_errors(errors,
                        vem_order,
                        tol,
                        test_type)
            if remove_folder:
                os.system("rm -rf " + os.path.join(program_folder, export_path))

    vem_types = [2]
    for vem_type in vem_types:
        for vem_order in vem_orders:
            for mesh_max_volume in mesh_max_volumes:
                export_path = run_program(program_folder,
                                          program_path,
                                          "Run_MG{0}".format(mesh_generator),
                                          vem_type,
                                          vem_order,
                                          test_type,
                                          mesh_generator,
                                          mesh_max_volume,
                                          1,
                                          solver_type,
                                          num_it_solver,
                                          tol_solver)
            errors = import_discrepancy_errors(export_path, vem_type, vem_order, test_type)
            test_discrepancy_errors(errors, 1.0e-03)
            if remove_folder:
                os.system("rm -rf " + os.path.join(program_folder, export_path))

    test_type = 4
    mesh_generator = 6
    mesh_max_volumes = [0.125 * 0.125 * 0.125, 0.0005787037]
    vem_types = [1]
    vem_orders = [2, 3]
    solver_type = 1
    tol_solver = 1.0e-12
    num_it_solver = 30
    for vem_type in vem_types:
        for vem_order in vem_orders:
            for mesh_max_volume in mesh_max_volumes:
                export_path = run_program(program_folder,
                                          program_path,
                                          "Run_MG{0}".format(mesh_generator),
                                          vem_type,
                                          vem_order,
                                          test_type,
                                          mesh_generator,
                                          mesh_max_volume,
                                          0,
                                          solver_type,
                                          num_it_solver,
                                          tol_solver)
            errors = import_errors(export_path, vem_type, vem_order, test_type)
            test_errors(errors,
                        vem_order,
                        tol,
                        test_type)
            if remove_folder:
                os.system("rm -rf " + os.path.join(program_folder, export_path))

    vem_types = [2]
    for vem_type in vem_types:
        for vem_order in vem_orders:
            for mesh_max_volume in mesh_max_volumes:
                export_path = run_program(program_folder,
                                          program_path,
                                          "Run_MG{0}".format(mesh_generator),
                                          vem_type,
                                          vem_order,
                                          test_type,
                                          mesh_generator,
                                          mesh_max_volume,
                                          1,
                                          solver_type,
                                          num_it_solver,
                                          tol_solver)
            errors = import_discrepancy_errors(export_path, vem_type, vem_order, test_type)
            test_discrepancy_errors(errors, 1.0e-03)
            if remove_folder:
                os.system("rm -rf " + os.path.join(program_folder, export_path))

    if remove_folder:
        os.system("rm -rf " + os.path.join(program_folder, export_folder))

    print("TESTS SUCCESS")

