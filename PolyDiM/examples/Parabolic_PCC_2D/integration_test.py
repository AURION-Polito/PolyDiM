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
                theta = 0.0,
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
    program_parameters += " Theta:double={0}".format(theta)
    program_parameters += " MeshImportFilePath:string={0}".format(mesh_import_path)

    output_file = os.path.join(program_folder,
                               "terminal.log")
    
    time_order = 1                       
    if theta == 0.5:
    	time_order = 2
    	

    run_label = "MethodType {0}".format(method_type)
    run_label += " SpaceMethodOrder {0}".format(method_order)
    run_label += " TimeMethodOrder {0}".format(time_order)
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
                               "Errors_" + str(test_type) + "_" + str(method_type) + "_" + str(method_order) + "_" + str(time_order) + ".csv")
    errors = []
    with open(errors_file, newline='') as csvfile:
        file_reader = csv.reader(csvfile, delimiter=';')
        data = list(file_reader)

        counter = 0
        for row in data:
            errors_row = []
            if counter == 0:
                errors_row.append(row[7])
                errors_row.append(row[10])
                errors_row.append(row[11])
                errors_row.append(row[12])
                errors_row.append(row[13])
            else:
                errors_row.append(float(row[7]))
                errors_row.append(float(row[10]))
                errors_row.append(float(row[11]))
                errors_row.append(float(row[12]))
                errors_row.append(float(row[13]))
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
    program_path = os.path.join(".", program_folder, "Parabolic_PCC_2D")

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
    thetas = [0, 1, 0.5]
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
    thetas = [0, 1, 0.5]
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
		        errors = import_errors(export_path, method_type, method_order, test_type)
		        test_errors(errors,
		                    method_order,
		                    tol)

		        if remove_folder:
		            os.system("rm -rf " + os.path.join(program_folder, export_path))

    test_type = 4
    mesh_generator = 5
    method_types = [0, 1, 2, 3]
    mesh_max_areas = [1./16., 1./64., 1./(16.*16), 1./(32.*32)]
    method_orders = [1, 2]
    thetas = [0, 1, 0.5]
    for theta in thetas:
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
		                                      theta=theta)
		            num_ref += 1

		        errors = import_errors(export_path, method_type, method_order, test_type)
		        test_errors(errors,
		                    method_order,
		                    tol)
		        if remove_folder:
		            os.system("rm -rf " + os.path.join(program_folder, export_path))

    print("TESTS SUCCESS")
