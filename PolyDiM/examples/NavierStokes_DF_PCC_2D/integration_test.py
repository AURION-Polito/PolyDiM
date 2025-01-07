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
                mesh_max_area,
                num_nl_iterations):
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
    program_parameters += " MeshMaxArea:double={0}".format(mesh_max_area)
    program_parameters += " ComputeVEMPerformance:bool={0}".format(0)
    program_parameters += " NLMaxNumberIterations:uint={0}".format(num_nl_iterations)
    program_parameters += " NLAbsChangeInSolutionTolerance:double={0}".format(1.0e-10)
    program_parameters += " NLAbsResidualTolerance:double={0}".format(1.0e-12)
    program_parameters += " NLRelChangeInSolutionTolerance:double={0}".format(1.0e-10)
    program_parameters += " NLRelResidualTolerance:double={0}".format(1.0e-12)

    output_file = os.path.join(program_folder,
                               "terminal.log")

    run_label = "VemType {0}".format(vem_type)
    run_label += " VemOrder {0}".format(vem_order)
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
        for col in data:
            errors_row = []
            if counter == 0:
                errors_row.append(col[6])
                errors_row.append(col[7])
                errors_row.append(col[8])
                errors_row.append(col[9])
                errors_row.append(col[10])
                errors_row.append(col[12])
                errors_row.append(col[13])
            else:
                errors_row.append(float(col[6]))
                errors_row.append(float(col[7]))
                errors_row.append(float(col[8]))
                errors_row.append(float(col[9]))
                errors_row.append(float(col[10]))
                errors_row.append(float(col[12]))
                errors_row.append(float(col[13]))
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
            if vem_order >= 4:
                print("Num. Ref. ", str(num_rows-1), " : [", abs(errors[0, 1]), " , ", abs(errors[1, 1]), " ], ",
                      "[", abs(errors[0, 2]), " , ", abs(errors[1, 2]), " ] ")
                assert abs(errors[0, 2]) < 1.0e-11
                assert abs(errors[1, 2]) < 1.0e-11
                assert abs(errors[0, 1]) < 1.0e-11
                assert abs(errors[1, 1]) < 1.0e-11
            else:
                slope_L2_pres = np.polyfit(np.log(errors[:, 0]), np.log(errors[:, 2]), 1)[0]
                print("Num. Ref. ", str(num_rows-1), " : [", abs(errors[0, 1]), " , " , abs(errors[1, 1]), " ], ", slope_L2_pres)
                assert round(slope_L2_pres) == round(float(vem_order))
                assert abs(errors[0, 1]) < tol
                assert abs(errors[1, 1]) < tol
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
            print("Num. Ref. ", str(num_rows-1), " : ", slope_H1_vel, slope_L2_pres)
            assert round(slope_L2_pres) >= round(float(vem_order))
            assert round(slope_H1_vel) >= round(float(vem_order))


if __name__ == "__main__":
    program_folder = os.path.dirname(os.path.realpath(__file__))
    program_path = os.path.join(".", program_folder, "NavierStokes_DF_PCC_2D")

    remove_folder = False

    vem_orders = [2, 3, 4]
    export_folder = "integration_tests"
    os.system("rm -rf " + os.path.join(program_folder, export_folder))
    tol = 5.0e-8

    print("RUN TESTS...")

    test_type = 1
    mesh_generator = 1
    mesh_max_area = 0.0
    vem_types = [1]
    for vem_type in vem_types:
        for vem_order in vem_orders:
            export_path = run_program(program_folder,
                                      program_path,
                                      "Run_MG{0}".format(mesh_generator),
                                      vem_type,
                                      vem_order,
                                      test_type,
                                      mesh_generator,
                                      mesh_max_area,
                                      1)
            errors = import_errors(export_path, vem_type, vem_order, test_type)
            test_errors(errors,
                        vem_order,
                        tol,
                        test_type)
            if remove_folder:
                os.system("rm -rf " + os.path.join(program_folder, export_path))

    if remove_folder:
        os.system("rm -rf " + os.path.join(program_folder, export_folder))

    print("TESTS SUCCESS")

