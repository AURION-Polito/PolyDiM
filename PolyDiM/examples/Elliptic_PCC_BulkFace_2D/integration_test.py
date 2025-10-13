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
                mesh_max_area = 0.1,
                supg = False,
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
    program_parameters += " MethodOrder:uint={0}".format(method_order)
    program_parameters += " ExportFolder:string={0}".format(export_path)
    program_parameters += " TestType:uint={0}".format(test_type)
    program_parameters += " MeshGenerator:uint={0}".format(mesh_generator)
    program_parameters += " MeshMaxArea:double={0}".format(mesh_max_area)
    program_parameters += " ComputeMethodPerformance:bool={0}".format(0)
    program_parameters += " SUPG:bool={0}".format(supg)
    program_parameters += " PecletConstant:double={0}".format(1.0/3.0)
    program_parameters += " MeshImportFilePath:string={0}".format(mesh_import_path)

    output_file = os.path.join(program_folder,
                               "terminal.log")

    run_label = "MethodType {0}".format(method_type)
    run_label += " MethodOrder {0}".format(method_order)
    run_label += " TestType {0}".format(test_type)
    run_label += " MeshGenerator {0}".format(mesh_generator)
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
                errors_row.append(row[6])
                errors_row.append(row[7])
                errors_row.append(row[8])
                errors_row.append(row[9])
                errors_row.append(row[10])
            else:
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

    if num_rows == 2:
        print("Num. Ref. 1: ", abs(errors[1][1]) / abs(errors[1][3]), abs(errors[1][2]) / abs(errors[1][4]))
        assert abs(errors[1][1]) < tol * abs(errors[1][3])
        assert abs(errors[1][2]) < tol * abs(errors[1][4])
    else:
        errors = np.array(errors[1:])
        slope_L2 = np.polyfit(np.log(errors[:, 0]), np.log(errors[:, 1]), 1)[0]
        slope_H1 = np.polyfit(np.log(errors[:, 0]), np.log(errors[:, 2]), 1)[0]
        print("Num. Ref. ", str(num_rows-1), ": ", slope_L2, slope_H1)
        assert round(slope_L2) == round(float(method_order + 1.0))
        assert round(slope_H1) == round(float(method_order))


if __name__ == "__main__":
    program_folder = os.path.dirname(os.path.realpath(__file__))
    program_path = os.path.join(".", program_folder, "Elliptic_PCC_2D")

    remove_folder = False

    export_folder = "integration_tests"
    os.system("rm -rf " + os.path.join(program_folder, export_folder))
    tol = 1.0e-12

    print("RUN TESTS...")

    test_type = 1
    mesh_generator = 1
    mesh_max_area = 0.0
    method_types = [0, 1, 2, 3]
    method_orders = [1, 2, 3]
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
                                      mesh_max_area=mesh_max_area)
            errors = import_errors(export_path, method_type, method_order, test_type)
            test_errors(errors,
                        method_order,
                        tol)

            if remove_folder:
                os.system("rm -rf " + os.path.join(program_folder, export_path))

    test_type = 1
    mesh_generator = 0
    mesh_max_area = 0.1
    method_types = [0]
    method_orders = [1, 2, 3]
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
                                      mesh_max_area=mesh_max_area)
            errors = import_errors(export_path, method_type, method_order, test_type)
            test_errors(errors,
                        method_order,
                        tol)

            if remove_folder:
                os.system("rm -rf " + os.path.join(program_folder, export_path))

    test_type = 2
    mesh_generator = 0
    method_types = [0, 1, 2, 3]
    mesh_max_areas = [0.01, 0.001]
    method_orders = [1, 2, 3]
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
                                          mesh_max_area=mesh_max_area)
                num_ref += 1
            errors = import_errors(export_path, method_type, method_order, test_type)
            test_errors(errors,
                        method_order,
                        tol)
            if remove_folder:
                os.system("rm -rf " + os.path.join(program_folder, export_path))

    test_type = 2
    mesh_generator = 2
    method_types = [1, 2, 3]
    mesh_max_areas = [0.01, 0.001]
    method_orders = [1, 2, 3]
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
                                          mesh_max_area=mesh_max_area)
                num_ref += 1

            errors = import_errors(export_path, method_type, method_order, test_type)
            test_errors(errors,
                        method_order,
                        tol)
            if remove_folder:
                os.system("rm -rf " + os.path.join(program_folder, export_path))

    test_type = 4
    mesh_generator = 5
    method_types = [0, 1, 2, 3]
    mesh_max_areas = [0.01, 0.001]
    method_orders = [1, 2, 3]
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
                                          mesh_max_area=mesh_max_area)
                num_ref += 1

            errors = import_errors(export_path, method_type, method_order, test_type)
            test_errors(errors,
                        method_order,
                        tol)
            if remove_folder:
                os.system("rm -rf " + os.path.join(program_folder, export_path))

    test_type = 2
    mesh_generator = 4
    mesh_import_paths = ["../../../../Mesh/2D/GenericPolyMesh"]
    method_types = [1, 2, 3]
    method_orders = np.arange(1, 13)
    list_errors = []
    for method_type in method_types:
        tab_errors = np.zeros([len(method_orders), 5])
        for method_order in method_orders:
            export_path = run_program(program_folder,
                                      program_path,
                                      "Run_MG{0}".format(mesh_generator),
                                      method_type,
                                      method_order,
                                      test_type,
                                      mesh_generator,
                                      0,
                                      mesh_import_path=mesh_import_paths[0])

            errors = import_errors(export_path, method_type, method_order, test_type)
            tab_errors[method_order-1, :] = np.array(errors[1:])
            if remove_folder:
                os.system("rm -rf " + os.path.join(program_folder, export_path))

        list_errors.append(tab_errors)

    fig, ax = plt.subplots(figsize=(12, 12))
    for h in range(len(method_types)):
        errors = list_errors[h]
        num_rows = len(errors)

        if method_types[h] == 1:
            ax.plot(method_orders, errors[:, 2], '-k^', linewidth=2, markersize=12,
                    label="Mon")
        elif method_types[h] == 2:
            ax.plot(method_orders, errors[:, 2], '-ro', linewidth=2, markersize=12,
                    label="Inrt")
        elif method_types[h] == 3:
            ax.plot(method_orders, errors[:, 2], '-bs', linewidth=2, markersize=12,
                    label="Ortho")
        else:
            raise ValueError("Not valid method type")

    plt.legend(bbox_to_anchor=(0., 1.02, 1.0, 0.2), loc="lower left",
               mode="expand", borderaxespad=0, ncol=3, fontsize=30)

    plt.xlabel('$k$', fontsize=30)
    plt.ylabel('$e_1$', fontsize=30)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.yscale('log')
    plt.grid(True, which="both", ls="--")
    plt.ylim(None, 10)
    plt.savefig(export_folder + "/{}_decay_plot.png".format(test_type), bbox_inches='tight', dpi=300)
    plt.show()

    test_type = 3
    mesh_generator = 0
    method_types = [0, 1, 2, 3]
    mesh_max_areas = [0.01, 0.001]
    method_orders = [1, 2]
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
                                          supg = True)
                num_ref += 1

            errors = import_errors(export_path, method_type, method_order, test_type)
            test_errors(errors,
                        method_order,
                        tol)
            if remove_folder:
                os.system("rm -rf " + os.path.join(program_folder, export_path))

    test_type = 3
    mesh_generator = 2
    method_orders = [1, 2, 3]
    method_types = [1, 2, 3]
    mesh_max_areas = [0.01, 0.001]
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
                                          supg = True)
                num_ref += 1
            errors = import_errors(export_path, method_type, method_order, test_type)
            test_errors(errors,
                        method_order,
                        tol)
            if remove_folder:
                os.system("rm -rf " + os.path.join(program_folder, export_path))

    if remove_folder:
        os.system("rm -rf " + os.path.join(program_folder, export_folder))

    print("TESTS SUCCESS")