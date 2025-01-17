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
                mesh_max_area):
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

    program_parameters = "MethodType:uint={0}".format(vem_type)
    program_parameters += " MethodOrder:uint={0}".format(vem_order)
    program_parameters += " ExportFolder:string={0}".format(export_path)
    program_parameters += " TestType:uint={0}".format(test_type)
    program_parameters += " MeshGenerator:uint={0}".format(mesh_generator)
    program_parameters += " MeshMaxArea:double={0}".format(mesh_max_area)
    program_parameters += " ComputeMethodPerformance:bool={0}".format(0)

    output_file = os.path.join(program_folder,
                               "terminal.log")

    run_label = "MethodType {0}".format(vem_type)
    run_label += " MethodOrder {0}".format(vem_order)
    run_label += " TestType {0}".format(test_type)
    run_label += " MeshGenerator {0}".format(mesh_generator)
    run_label += " MeshMaxArea {0}".format(mesh_max_area)
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
        for row in data:
            errors_row = []
            if counter == 0:
                errors_row.append(row[6])
                errors_row.append(row[7])
                errors_row.append(row[8])
                errors_row.append(row[9])
                errors_row.append(row[10])
                errors_row.append(row[11])
            else:
                errors_row.append(float(row[6]))
                errors_row.append(float(row[7]))
                errors_row.append(float(row[8]))
                errors_row.append(float(row[9]))
                errors_row.append(float(row[10]))
                errors_row.append(float(row[11]))
            errors.append(errors_row)
            counter += 1

    return errors


def test_errors(errors,
                vem_order,
                vem_type,
                tol):
    num_rows = len(errors)

    if num_rows == 2:
        print("Num. Ref. 1: ", abs(errors[1][1]) / abs(errors[1][4]), abs(errors[1][2]) / abs(errors[1][5]),
              abs(errors[1][3]) / abs(errors[1][5]))
        assert abs(errors[1][1]) < tol * abs(errors[1][4])
        assert abs(errors[1][2]) < tol * abs(errors[1][5])
        assert abs(errors[1][3]) < tol * abs(errors[1][5])
    else:
        errors = np.array(errors[1:])
        slope_L2_vel = np.polyfit(np.log(errors[:, 0]), np.log(errors[:, 1]), 1)[0]
        slope_L2_pres = np.polyfit(np.log(errors[:, 0]), np.log(errors[:, 2]), 1)[0]
        slope_super_L2_pres = np.polyfit(np.log(errors[:, 0]), np.log(errors[:, 3]), 1)[0]
        print("Num. Ref. ", str(num_rows - 1), ": ", slope_L2_vel, slope_L2_pres, slope_super_L2_pres)
        if vem_order != 3:
            assert round(slope_L2_vel) == round(float(vem_order + 1.0))
        assert round(slope_L2_pres) == round(float(vem_order + 1.0))
        if vem_type != 5 and vem_type != 4:
            assert round(slope_super_L2_pres) >= round(float(vem_order + 2.0))


if __name__ == "__main__":
    program_folder = os.path.dirname(os.path.realpath(__file__))
    program_path = os.path.join(".", program_folder, "Elliptic_MCC_2D")

    remove_folder = True

    vem_types = [1, 2, 3, 4, 5]
    vem_orders = [0, 1, 2, 3]
    export_folder = "integration_tests"
    os.system("rm -rf " + os.path.join(program_folder, export_folder))
    tol = 1.0e-12

    print("RUN TESTS...")

    test_type = 1
    mesh_generator = 1
    mesh_max_area = 0.0
    for vem_type in vem_types:
        for vem_order in vem_orders:
            export_path = run_program(program_folder,
                                      program_path,
                                      "Run_MG{0}".format(mesh_generator),
                                      vem_type,
                                      vem_order,
                                      test_type,
                                      mesh_generator,
                                      mesh_max_area)
            errors = import_errors(export_path, vem_type, vem_order, test_type)
            test_errors(errors,
                        vem_order,
                        vem_type,
                        tol)
            if remove_folder:
                os.system("rm -rf " + os.path.join(program_folder, export_path))

    test_type = 2
    mesh_generator = 0
    vem_orders = [0, 1, 2]
    mesh_max_areas = [0.01, 0.001]
    for vem_type in vem_types:
        for vem_order in vem_orders:
            for mesh_max_area in mesh_max_areas:
                export_path = run_program(program_folder,
                                          program_path,
                                          "Run_MG{0}".format(mesh_generator),
                                          vem_type,
                                          vem_order,
                                          test_type,
                                          mesh_generator,
                                          mesh_max_area)
            errors = import_errors(export_path, vem_type, vem_order, test_type)
            test_errors(errors,
                        vem_order,
                        vem_type,
                        tol)
            if remove_folder:
                os.system("rm -rf " + os.path.join(program_folder, export_path))

    test_type = 2
    vem_types = [1, 2, 3, 4, 5]
    mesh_generator = 5
    mesh_max_areas = [0.01, 0.001]
    for vem_type in vem_types:
        for vem_order in vem_orders:
            for mesh_max_area in mesh_max_areas:
                export_path = run_program(program_folder,
                                          program_path,
                                          "Run_MG{0}".format(mesh_generator),
                                          vem_type,
                                          vem_order,
                                          test_type,
                                          mesh_generator,
                                          mesh_max_area)
            errors = import_errors(export_path, vem_type, vem_order, test_type)
            test_errors(errors,
                        vem_order,
                        vem_type,
                        tol)
            if remove_folder:
                os.system("rm -rf " + os.path.join(program_folder, export_path))

    if remove_folder:
        os.system("rm -rf " + os.path.join(program_folder, export_folder))

    print("TESTS SUCCESS")


