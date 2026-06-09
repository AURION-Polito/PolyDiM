import os
import csv
import numpy as np




def run_program(program_folder,
                program_path,
                run_folder,
                method_type,
                method_order,
                test_type,
                mesh_generator,
                num_ref,
                mesh_max_area=0.0,
                mesh_import_path="./",
                compute_discrepancy=False):
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
    program_parameters += " ComputeDiscrepancyError:bool={0}".format(compute_discrepancy)
    program_parameters += " MeshImportFilePath:string={0}".format(mesh_import_path)

    output_file = os.path.join(program_folder,
                               "terminal.log")

    run_label = "MethodType {0}".format(method_type)
    run_label += " MethodOrder {0}".format(method_order)
    run_label += " TestType {0}".format(test_type)
    run_label += " MeshGenerator {0}".format(mesh_generator)
    run_label += " NumRefinement {0}".format(num_ref)
    print("Run " + run_label + "...")
    os.system(program_path + " " + program_parameters + "> " + output_file)
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


def import_discrepancy_errors(export_path, method_type, method_order, test_type):
    errors_file = os.path.join(export_path,
                               "Solution",
                               "DiscrepancyErrors_" + str(test_type) + "_" + str(method_type) + "_" + str(
                                   method_order) + ".csv")
    errors = []
    with open(errors_file, newline='') as csvfile:
        file_reader = csv.reader(csvfile, delimiter=';')
        data = list(file_reader)

        counter = 0
        for col in data:
            errors_row = []
            if counter == 0:
                errors_row.append(col[6])  # discrepancy h1 velocity
                errors_row.append(col[7])  # discrepancy l2 pressure
                errors_row.append(col[8])  # norm h1 full velocity
                errors_row.append(col[9])  # norm l2 full pressure
            else:
                errors_row.append(float(col[6]))
                errors_row.append(float(col[7]))
                errors_row.append(float(col[8]))
                errors_row.append(float(col[9]))
            errors.append(errors_row)
            counter += 1

    return errors


def read_flux(export_path, method_type, method_order, test_type):
    errors_file = os.path.join(export_path, "Solution",
                               "Flux_" + str(test_type) + "_" + str(method_type) + "_" + str(method_order) + ".csv")

    flux = []
    with open(errors_file, newline='') as csvfile:
        file_reader = csv.reader(csvfile, delimiter=';')
        data = list(file_reader)

        data.pop(0)  # remove header

        for row in data:
            flux.append(float(row[1]))

    return np.array(flux)


def errors_test(errors,
                method_order,
                method_type,
                tol,
                test_type):
    num_rows = len(errors)

    if num_rows == 2:
        print("Patch: ", abs(errors[1][1]) / abs(errors[1][3]), abs(errors[1][2]) / abs(errors[1][4]))
        assert abs(errors[1][1]) < tol * abs(errors[1][3])
        if method_type != 2:
            assert abs(errors[1][2]) < tol * abs(errors[1][4])
    else:
        if test_type == 3:
            errors = np.array(errors[1:])
            if method_order >= 4:
                print("Num. Ref. ", str(num_rows - 1), " : [", abs(errors[0, 1]), " , ", abs(errors[1, 1]), " ], ",
                      "[", abs(errors[0, 2]), " , ", abs(errors[1, 2]), " ] ")
                assert abs(errors[0, 2]) < 1.0e-11
                assert abs(errors[1, 2]) < 1.0e-11
                assert abs(errors[0, 1]) < 1.0e-11
                assert abs(errors[1, 1]) < 1.0e-11
            else:
                slope_l2_pres = float(np.polyfit(np.log(errors[:, 0]), np.log(errors[:, 2]), 1)[0])
                print("Num. Ref. ", str(num_rows - 1), " : [", abs(errors[0, 1]), " , ", abs(errors[1, 1]), " ], ",
                      slope_l2_pres)
                assert round(slope_l2_pres) == round(float(method_order))
                assert abs(errors[0, 1]) < tol
                assert abs(errors[1, 1]) < tol
        elif test_type == 4:
            errors = np.array(errors[1:])
            slope_l2_pres = float(np.polyfit(np.log(errors[:, 0]), np.log(errors[:, 2]), 1)[0])
            slope_h1_vel = float(np.polyfit(np.log(errors[:, 0]), np.log(errors[:, 1]), 1)[0])
            print("Num. Ref. ", str(num_rows - 1), " : ", slope_h1_vel, slope_l2_pres)
            assert round(slope_l2_pres) == round(float(method_order))
            assert round(slope_h1_vel) >= round(float(method_order + 2))
        else:
            errors = np.array(errors[1:])
            slope_l2_pres = float(np.polyfit(np.log(errors[:, 0]), np.log(errors[:, 2]), 1)[0])
            slope_h1_vel = float(np.polyfit(np.log(errors[:, 0]), np.log(errors[:, 1]), 1)[0])
            print("Num. Ref. ", str(num_rows - 1), " : ", slope_h1_vel, slope_l2_pres)
            assert round(slope_l2_pres) >= round(float(method_order))
            assert round(slope_h1_vel) >= round(float(method_order))


def discrepancy_errors_test(errors, tol):
    num_rows = len(errors)

    for r in range(num_rows - 1):
        print("Discrepancy - Id. Ref. ", str(r), " : ", abs(errors[r + 1][0]) / abs(errors[r + 1][2]),
              abs(errors[r + 1][1]) / abs(errors[r + 1][3]))
        assert abs(errors[r + 1][0]) < tol * abs(errors[r + 1][2])
        assert abs(errors[r + 1][1]) < tol * abs(errors[r + 1][3])


def main():
    program_folder = os.path.dirname(os.path.realpath(__file__))
    program_path = os.path.join(".", program_folder, "Brinkman_DF_PCC_2D")

    remove_folder = False

    method_orders = [2, 3, 4]
    os.system("rm -rf " + os.path.join(program_folder, export_folder))
    tol = 5.0e-8

    print("RUN TESTS...")

    test_type = 1
    mesh_generator = 1
    mesh_max_area = 0.0
    method_types = [1]
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
                                      compute_discrepancy=False)
            errors = import_errors(export_path, method_type, method_order, test_type)
            errors_test(errors,
                        method_order,
                        method_type,
                        tol,
                        test_type)
            if remove_folder:
                os.system("rm -rf " + os.path.join(program_folder, export_path))

    test_type = 2
    mesh_generator = 5
    mesh_max_areas = [0.125 * 0.125, 0.0625 * 0.0625]
    method_types = [1]
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
                                          compute_discrepancy=False)
                num_ref += 1
            errors = import_errors(export_path, method_type, method_order, test_type)
            errors_test(errors,
                        method_order,
                        method_type,
                        tol,
                        test_type)
            if remove_folder:
                os.system("rm -rf " + os.path.join(program_folder, export_path))

    method_types = [2]
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
                                          compute_discrepancy=True)
                num_ref += 1
            errors = import_discrepancy_errors(export_path, method_type, method_order, test_type)
            discrepancy_errors_test(errors, tol)
            if remove_folder:
                os.system("rm -rf " + os.path.join(program_folder, export_path))

    test_type = 2
    mesh_generator = 2
    mesh_max_areas = [0.01, 0.005]
    method_types = [1]
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
                                          compute_discrepancy=False)
                num_ref += 1
            errors = import_errors(export_path, method_type, method_order, test_type)
            errors_test(errors,
                        method_order,
                        method_type,
                        tol,
                        test_type)
            if remove_folder:
                os.system("rm -rf " + os.path.join(program_folder, export_path))

    method_types = [2]
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
                                          compute_discrepancy=True)
                num_ref += 1
            errors = import_discrepancy_errors(export_path, method_type, method_order, test_type)
            discrepancy_errors_test(errors, tol)
            if remove_folder:
                os.system("rm -rf " + os.path.join(program_folder, export_path))

    test_type = 3
    mesh_generator = 5
    mesh_max_areas = [0.125 * 0.125, 0.0625 * 0.0625]
    method_types = [1]
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
                                          compute_discrepancy=False)
                num_ref += 1
            errors = import_errors(export_path, method_type, method_order, test_type)
            errors_test(errors,
                        method_order,
                        method_type,
                        tol,
                        test_type)
            if remove_folder:
                os.system("rm -rf " + os.path.join(program_folder, export_path))

    test_type = 3
    mesh_generator = 2
    mesh_max_areas = [0.01, 0.005]
    method_types = [1]
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
                                          compute_discrepancy=False)

                num_ref += 1
            errors = import_errors(export_path, method_type, method_order, test_type)
            errors_test(errors,
                        method_order,
                        method_type,
                        tol,
                        test_type)
            if remove_folder:
                os.system("rm -rf " + os.path.join(program_folder, export_path))

    test_type = 4
    mesh_generator = 5
    mesh_max_areas = [0.125 * 0.125, 0.0625 * 0.0625]
    method_types = [1]
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
                                          compute_discrepancy=False)
                num_ref += 1
            errors = import_errors(export_path, method_type, method_order, test_type)
            errors_test(errors,
                        method_order,
                        method_type,
                        tol,
                        test_type)
            if remove_folder:
                os.system("rm -rf " + os.path.join(program_folder, export_path))

    test_type = 4
    mesh_generator = 2
    mesh_max_areas = [0.01, 0.005]
    method_types = [1]
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
                                          compute_discrepancy=False)
                num_ref += 1
            errors = import_errors(export_path, method_type, method_order, test_type)
            errors_test(errors,
                        method_order,
                        method_type,
                        tol,
                        test_type)
            if remove_folder:
                os.system("rm -rf " + os.path.join(program_folder, export_path))

    test_type = 5
    mesh_generator = 5
    mesh_max_areas = [0.125 * 0.125, 0.0625 * 0.0625]
    method_types = [1]
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
                                          compute_discrepancy=False)
                num_ref += 1
            errors = import_errors(export_path, method_type, method_order, test_type)
            errors_test(errors,
                        method_order,
                        method_type,
                        tol,
                        test_type)
            if remove_folder:
                os.system("rm -rf " + os.path.join(program_folder, export_path))

    test_type = 5
    mesh_generator = 2
    mesh_max_areas = [0.01, 0.005]
    method_types = [1]
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
                                          compute_discrepancy=False)
                num_ref += 1
            errors = import_errors(export_path, method_type, method_order, test_type)
            errors_test(errors,
                        method_order,
                        method_type,
                        tol,
                        test_type)
            if remove_folder:
                os.system("rm -rf " + os.path.join(program_folder, export_path))

    test_type = 6
    mesh_generator = 5
    mesh_max_areas = [0.125 * 0.125, 0.0625 * 0.0625]
    method_types = [1]
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
                                          compute_discrepancy=False)
                num_ref += 1
            errors = import_errors(export_path, method_type, method_order, test_type)
            errors_test(errors,
                        method_order,
                        method_type,
                        tol,
                        test_type)
            if remove_folder:
                os.system("rm -rf " + os.path.join(program_folder, export_path))

    method_types = [2]
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
                                          compute_discrepancy=True)
                num_ref += 1
            errors = import_discrepancy_errors(export_path, method_type, method_order, test_type)
            discrepancy_errors_test(errors, tol)
            if remove_folder:
                os.system("rm -rf " + os.path.join(program_folder, export_path))

    test_type = 6
    mesh_generator = 2
    mesh_max_areas = [0.005, 0.001]
    method_types = [1]
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
                                          compute_discrepancy=False)
                num_ref += 1
            errors = import_errors(export_path, method_type, method_order, test_type)
            errors_test(errors,
                        method_order,
                        method_type,
                        tol,
                        test_type)
            if remove_folder:
                os.system("rm -rf " + os.path.join(program_folder, export_path))

    method_types = [2]
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
                                          compute_discrepancy=True)
                num_ref += 1
            errors = import_discrepancy_errors(export_path, method_type, method_order, test_type)
            discrepancy_errors_test(errors, tol)
            if remove_folder:
                os.system("rm -rf " + os.path.join(program_folder, export_path))

    test_type = 7
    mesh_generator = 4
    mesh_names = [program_folder + "/../../../../Mesh/2D/DarcyStokesMesh/TriMesh_0_01",
                  program_folder + "/../../../../Mesh/2D/DarcyStokesMesh/Square_1_8"]
    method_types = [1, 2]
    for method_type in method_types:
        for method_order in method_orders:
            for idx_mesh in range(len(mesh_names)):
                export_path = run_program(program_folder,
                                          program_path,
                                          "Run_MG{0}_{1}".format(mesh_generator, idx_mesh),
                                          method_type,
                                          method_order,
                                          test_type,
                                          mesh_generator,
                                          0,
                                          mesh_import_path=mesh_names[idx_mesh],
                                          compute_discrepancy=False)

                errors = import_errors(export_path, method_type, method_order, test_type)
                errors_test(errors,
                            method_order,
                            method_type,
                            tol,
                            test_type)

                flux = read_flux(export_path, method_type, method_order, test_type)
                exact_flux = np.array([0., 0., 2., -2., 0., -1., 0., -1., 2.])

                assert np.linalg.norm(flux - exact_flux) < tol

                if remove_folder:
                    os.system("rm -rf " + os.path.join(program_folder, export_path))

    if remove_folder:
        os.system("rm -rf " + os.path.join(program_folder, export_folder))

    print("TESTS SUCCESS")


if __name__ == "__main__":
    export_folder = "integration_tests"
    main()