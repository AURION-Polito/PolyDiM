import os
import csv
import matplotlib.pyplot as plt
import numpy as np


def run_program(program_folder,
                program_path,
                run_folder,
                method_type,
                method_order,
                test_type,
                mesh_generator,
                num_ref,
                time_step = 1.0,
                mesh_max_area = 0.1,
                mesh_import_path = "./",
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
                                   method_order))

    program_parameters = "MethodType:uint={0}".format(method_type)
    program_parameters += " ComputeMethodPerformance:bool={0}".format(0)
    program_parameters += " MethodOrder:uint={0}".format(method_order)
    program_parameters += " ExportFolder:string={0}".format(export_path)
    program_parameters += " TestType:uint={0}".format(test_type)
    program_parameters += " MeshGenerator:uint={0}".format(mesh_generator)
    program_parameters += " MeshMaxArea:double={0}".format(mesh_max_area)
    program_parameters += " TimeStep:double={0}".format(time_step)
    program_parameters += " MeshImportFilePath:string={0}".format(mesh_import_path)

    output_file = os.path.join(program_folder,
                               "terminal.log")

    run_label = "MethodType {0}".format(method_type)
    run_label += " MethodOrder {0}".format(method_order)
    run_label += " TestType {0}".format(test_type)
    run_label += " TimeStep {:<.4e}".format(time_step)
    run_label += " NumRefinement {0}".format(num_ref)
    print("Run " + run_label + "...")
    os.system(program_path + " " + program_parameters + " > " + output_file)
    os.system("mv " + output_file + " " + export_path)
    print("Run SUCCESS")

    return export_path


def import_errors(export_path, method_type, method_order, test_type):
    errors_file = os.path.join(export_path,
                               "Solution",
                               "Errors_" + str(test_type) + "_" + str(method_type) + "_" + str(method_order) + ".csv")
    errors = []
    with open(errors_file, newline='') as csvfile:
        file_reader = csv.reader(csvfile, delimiter=';')
        data = list(file_reader)

        counter = 0
        for row in data:
            errors_row = []
            if counter == 0:
                errors_row.append(row[5]) # h
                errors_row.append(row[6]) # delta_time
                errors_row.append(row[7]) # errorL2
                errors_row.append(row[8]) # errorH1
                errors_row.append(row[9]) # normL2
                errors_row.append(row[10]) # normH1
            else:
                errors_row.append(float(row[5]))
                errors_row.append(float(row[6]))
                errors_row.append(float(row[7]))
                errors_row.append(float(row[8]))
                errors_row.append(float(row[9]))
                errors_row.append(float(row[10]))
            errors.append(errors_row)
            counter += 1

    return errors


def test_errors(errors,
                method_order,
                tol):
    num_rows = len(errors)

    print("{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}".format("h", "deltaT", "errorL2", "EOC-L2", "errorH1", "EOC-H1"))
    print("{:<10.2e}{:<10.2e}{:<10.2e}{:<10}{:<10.2e}{:<10}".format(errors[0, 0], errors[0, 1],  errors[0, 2], "-", errors[0, 3], "-"))
    for i in range(num_rows-1):
        slope_L2 = np.polyfit(np.log(errors[0:i+2, 0]), np.log(errors[0:i+2, 2]), 1)[0]
        slope_H1 = np.polyfit(np.log(errors[0:i+2, 0]), np.log(errors[0:i+2, 3]), 1)[0]
        print("{:<10.2e}{:<10.2e}{:<10.2e}{:<10.2e}{:<10.2e}{:<10.2e}".format(errors[i+1, 0], errors[i+1, 1],  errors[i+1, 2], slope_L2, errors[i+1, 3], slope_H1))
        assert round(slope_L2) >= round(float(method_order + 1.0))
        assert round(slope_H1) >= round(float(method_order))


if __name__ == "__main__":
    program_folder = os.path.dirname(os.path.realpath(__file__))
    program_path = os.path.join(".", program_folder, "Parabolic_PCC_BulkFace_2D")

    remove_folder = False

    export_folder = "integration_tests"
    os.system("rm -rf " + os.path.join(program_folder, export_folder))
    tol = 1.0e-12

    print("RUN TESTS...")

    test_type = 1
    mesh_generator = 4
    mesh_import_paths = [program_folder + "/../../../../Mesh/2D/CircleTriangularMesh/CircleTriangularMesh_R1",
                         program_folder + "/../../../../Mesh/2D/CircleTriangularMesh/CircleTriangularMesh_R2",
                         program_folder + "/../../../../Mesh/2D/CircleTriangularMesh/CircleTriangularMesh_R3",
                         program_folder + "/../../../../Mesh/2D/CircleTriangularMesh/CircleTriangularMesh_R4"]
    method_types = [0, 1, 2, 3]
    method_order = 1
    list_errors = []
    time_step = 1.0
    for method_type in method_types:
        num_ref = 0
        for mesh_import_path in mesh_import_paths:
            export_path = run_program(program_folder,
                                      program_path,
                                      "Run_MG{0}".format(mesh_generator),
                                      method_type,
                                      method_order,
                                      test_type,
                                      mesh_generator,
                                      num_ref,
                                      time_step = time_step,
                                      mesh_import_path=mesh_import_path)
            num_ref += 1
        errors = import_errors(export_path, method_type, method_order, test_type)
        tab_errors = np.array(errors[1:])

        test_errors(tab_errors,  method_order,  tol)

        if remove_folder:
            os.system("rm -rf " + os.path.join(program_folder, export_path))

        list_errors.append(tab_errors)

    test_type = 2
    mesh_generator = 4
    mesh_import_paths = [program_folder + "/../../../../Mesh/2D/CircleTriangularMesh/CircleTriangularMesh_R1",
                         program_folder + "/../../../../Mesh/2D/CircleTriangularMesh/CircleTriangularMesh_R2",
                         program_folder + "/../../../../Mesh/2D/CircleTriangularMesh/CircleTriangularMesh_R3",
                         program_folder + "/../../../../Mesh/2D/CircleTriangularMesh/CircleTriangularMesh_R4"]
    method_types = [0, 1, 2, 3]
    method_order = 1
    list_errors = []
    time_steps = [2.500e-3, 6.250e-4, 1.563e-4, 1.563e-4]
    for method_type in method_types:
        num_ref = 0
        for mesh_import_path in mesh_import_paths:
            export_path = run_program(program_folder,
                                      program_path,
                                      "Run_MG{0}".format(mesh_generator),
                                      method_type,
                                      method_order,
                                      test_type,
                                      mesh_generator,
                                      num_ref,
                                      time_step = time_steps[num_ref],
                                      mesh_import_path=mesh_import_path)

            num_ref += 1
        errors = import_errors(export_path, method_type, method_order, test_type)
        tab_errors = np.array(errors[1:])

        test_errors(tab_errors,  method_order,  tol)

        if remove_folder:
            os.system("rm -rf " + os.path.join(program_folder, export_path))

        list_errors.append(tab_errors)



    if remove_folder:
        os.system("rm -rf " + os.path.join(program_folder, export_folder))

    print("TESTS SUCCESS")