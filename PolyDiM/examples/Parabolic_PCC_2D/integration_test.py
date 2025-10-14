import os
import csv
import matplotlib.pyplot as plt
import numpy as np


def time_order(theta):
    if theta == 0.5:
        return 2
    else:
        return 1


def run_program(program_folder,
                program_path,
                run_folder,
                method_type,
                method_order,
                test_type,
                mesh_generator,
                num_ref,
                mesh_max_area=0.1,
                time_step=0.5,
                theta=0.0,
                mesh_import_path="./",
                ):
    export_path = os.path.join(program_folder,
                               export_folder,
                               "{0}_TT{1}".format(
                                   run_folder,
                                   test_type),
                               "{0}_TT{1}_VT{2}".format(
                                   run_folder,
                                   test_type,
                                   method_type),
                               "{0}_TT{1}_VT{2}_VO{3}".format(
                                   run_folder,
                                   test_type,
                                   method_type,
                                   method_order),
                               "{0}_TT{1}_VT{2}_VO{3}_TO{4}".format(
                                   run_folder,
                                   test_type,
                                   method_type,
                                   method_order,
                                   time_order(theta)))

    program_parameters = "MethodType:uint={0}".format(method_type)
    program_parameters += " MethodOrder:uint={0}".format(method_order)
    program_parameters += " ExportFolder:string={0}".format(export_path)
    program_parameters += " TestType:uint={0}".format(test_type)
    program_parameters += " MeshGenerator:uint={0}".format(mesh_generator)
    program_parameters += " MeshMaxArea:double={0}".format(mesh_max_area)
    program_parameters += " TimeStep:double={0}".format(time_step)
    program_parameters += " Theta:double={0}".format(theta)
    program_parameters += " MeshImportFilePath:string={0}".format(mesh_import_path)

    output_file = os.path.join(program_folder,
                               "terminal.log")

    run_label = "MethodType {0}".format(method_type)
    run_label += " SpaceMethodOrder {0}".format(method_order)
    run_label += " TimeMethodOrder {0}".format(time_order(theta))
    run_label += " TestType {0}".format(test_type)
    run_label += " MeshGenerator {0}".format(mesh_generator)
    run_label += " NumRefinement {0}".format(num_ref)
    print("Run " + run_label + "...")
    os.system(program_path + " " + program_parameters + " > " + output_file)
    os.system("mv " + output_file + " " + export_path)
    print("Run SUCCESS")

    return export_path


def import_errors(export_path, method_type, method_order, time_order, test_type):
    errors_file = os.path.join(export_path,
                               "Solution",
                               "Errors_" + str(test_type) + "_" + str(method_type) + "_" + str(
                                   method_order) + "_" + str(time_order) + ".csv")
    errors = []
    with open(errors_file, newline='') as csvfile:
        file_reader = csv.reader(csvfile, delimiter=';')
        data = list(file_reader)

        counter = 0
        for row in data:
            errors_row = []
            if counter == 0:
                errors_row.append(row[7])
                errors_row.append(row[8])
                errors_row.append(row[9])
                errors_row.append(row[10])
                errors_row.append(row[11])
                errors_row.append(row[12])
                errors_row.append(row[13])
            else:
                errors_row.append(float(row[7]))
                errors_row.append(float(row[8]))
                errors_row.append(float(row[9]))
                errors_row.append(float(row[10]))
                errors_row.append(float(row[11]))
                errors_row.append(float(row[12]))
                errors_row.append(float(row[13]))
            errors.append(errors_row)
            counter += 1

    return errors


def test_space_errors(errors,
                      method_order,
                      max_time,
                      tol):
    errors = np.array(errors[1:])
    T_rows = np.where((errors[:, 2] == max_time))[0]
    num_rows = T_rows.size

    if num_rows == 1:
        row_number = T_rows[0]
        print("Num. Ref. 1: ", abs(errors[row_number, 3]) / abs(errors[row_number, 5]),
              abs(errors[row_number, 4]) / abs(errors[row_number, 6]))
        assert abs(errors[row_number, 3]) < tol * abs(errors[row_number, 5])
        assert abs(errors[row_number, 4]) < tol * abs(errors[row_number, 6])
    else:
        slope_L2 = np.polyfit(np.log(errors[T_rows, 0]), np.log(errors[T_rows, 3]), 1)[0]
        slope_H1 = np.polyfit(np.log(errors[T_rows, 0]), np.log(errors[T_rows, 4]), 1)[0]
        print("Num. Ref. ", str(num_rows - 1), ": ", slope_L2, slope_H1)
        assert round(slope_L2) == round(float(method_order + 1.0))
        assert round(slope_H1) == round(float(method_order))


def test_time_errors(errors,
                     method_order,
                     max_time,
                     tol):
    errors = np.array(errors[1:])
    T_rows = np.where((errors[:, 2] == max_time))[0]
    num_rows = T_rows.size

    if num_rows == 1:
        row_number = T_rows[0]
        print("Num. Ref. 1: ", abs(errors[row_number, 3]) / abs(errors[row_number, 5]),
              abs(errors[row_number, 4]) / abs(errors[row_number, 6]))
        assert abs(errors[row_number, 3]) < tol * abs(errors[row_number, 5])
        assert abs(errors[row_number, 4]) < tol * abs(errors[row_number, 6])
    else:
        slope_L2 = np.polyfit(np.log(1.0 / errors[T_rows, 1]), np.log(errors[T_rows, 3]), 1)[0]
        slope_H1 = np.polyfit(np.log(1.0 / errors[T_rows, 1]), np.log(errors[T_rows, 4]), 1)[0]
        print("Num. Ref. ", str(num_rows - 1), ": ", slope_L2, slope_H1)
        assert round(slope_L2) == round(float(method_order))
        assert round(slope_H1) == round(float(method_order))


if __name__ == "__main__":
    program_folder = os.path.dirname(os.path.realpath(__file__))
    program_path = os.path.join(".", program_folder, "Parabolic_PCC_2D")

    remove_folder = True

    export_folder = "integration_tests"
    os.system("rm -rf " + os.path.join(program_folder, export_folder))
    tol = 1.0e-10

    print("RUN TESTS...")

    test_type = 1
    mesh_generator = 1
    mesh_max_area = 0.0
    method_types = [0, 1, 2, 3]
    method_orders = [1, 2, 3]
    thetas = [1, 0.5]
    for theta in thetas:
        for method_type in method_types:
            for method_order in method_orders:
                export_path = run_program(program_folder,
                                          program_path,
                                          "Run_MG{0}".format(mesh_generator),
                                          method_type,
                                          method_order,
                                          test_type,
                                          mesh_generator,
                                          0,
                                          mesh_max_area=mesh_max_area,
                                          theta=theta)
                errors = import_errors(export_path, method_type, method_order, time_order(theta), test_type)
                test_space_errors(errors,
                                  method_order,
                                  1.0,
                                  tol)

                if remove_folder:
                    os.system("rm -rf " + os.path.join(program_folder, export_path))

    test_type = 1
    mesh_generator = 0
    mesh_max_area = 0.1
    method_types = [0]
    method_orders = [1, 2, 3]
    thetas = [1, 0.5]
    for theta in thetas:
        for method_type in method_types:
            for method_order in method_orders:
                export_path = run_program(program_folder,
                                          program_path,
                                          "Run_MG{0}".format(mesh_generator),
                                          method_type,
                                          method_order,
                                          test_type,
                                          mesh_generator,
                                          0,
                                          mesh_max_area=mesh_max_area,
                                          theta=theta)
                errors = import_errors(export_path, method_type, method_order, time_order(theta), test_type)
                test_space_errors(errors,
                                  method_order,
                                  1.0,
                                  tol)

                if remove_folder:
                    os.system("rm -rf " + os.path.join(program_folder, export_path))

    test_type = 3
    mesh_generator = 0
    mesh_max_area = 0.1
    method_types = [0]
    method_orders = [1, 2, 3]
    time_steps = [1. / 5., 1. / 10., 1. / 20]
    thetas = [1, 0.5]
    for theta in thetas:
        for method_type in method_types:
            for method_order in method_orders:
                for time_step in time_steps:
                    export_path = run_program(program_folder,
                                              program_path,
                                              "Run_MG{0}".format(mesh_generator),
                                              method_type,
                                              method_order,
                                              test_type,
                                              mesh_generator,
                                              0,
                                              mesh_max_area=mesh_max_area,
                                              time_step=time_step,
                                              theta=theta)

                errors = import_errors(export_path, method_type, method_order, time_order(theta), test_type)
                test_time_errors(errors,
                                 time_order(theta),
                                 1.0,
                                 tol)

                if remove_folder:
                    os.system("rm -rf " + os.path.join(program_folder, export_path))

    test_type = 4
    mesh_generator = 5
    method_types = [0, 1, 2, 3]
    mesh_max_areas = [1. / 16., 1. / 64., 1. / (16. * 16), 1. / (32. * 32)]
    method_orders = [1, 2]
    time_steps = [1. / 40.]
    thetas = [0.5]
    for theta in thetas:
        for time_step in time_steps:
            for method_type in method_types:
                for method_order in method_orders:
                    num_ref = 0
                    for mesh_max_area in mesh_max_areas:
                        export_path = run_program(program_folder,
                                                  program_path,
                                                  "Run_MG{0}".format(mesh_generator),
                                                  method_type,
                                                  method_order,
                                                  test_type,
                                                  mesh_generator,
                                                  num_ref,
                                                  mesh_max_area=mesh_max_area,
                                                  time_step=time_step,
                                                  theta=theta)
                        num_ref += 1

                    errors = import_errors(export_path, method_type, method_order, time_order(theta), test_type)
                    test_space_errors(errors,
                                      method_order,
                                      1.0,
                                      tol)

                    if remove_folder:
                        os.system("rm -rf " + os.path.join(program_folder, export_path))

    print("TESTS SUCCESS")
