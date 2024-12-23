import os
import csv
import math

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
    
    program_parameters = "VemType:uint={0}".format(vem_type)
    program_parameters += " VemOrder:uint={0}".format(vem_order)
    program_parameters += " ExportFolder:string={0}".format(export_path)
    program_parameters += " TestType:uint={0}".format(test_type)
    program_parameters += " MeshGenerator:uint={0}".format(mesh_generator)
    program_parameters += " MeshMaxArea:double={0}".format(mesh_max_area)
    program_parameters += " ComputeVEMPerformance:bool={0}".format(0)

    output_file = os.path.join(program_folder,
                               "terminal.log")

    print("Run " + export_path + "...")
    os.system(program_path + " " + program_parameters + "> " + output_file)
    os.system("mv " + output_file + " " + export_path)
    print("Run SUCCESS")

    return export_path

def import_errors(export_path):
    errors_file = os.path.join(export_path,
                               "Solution",
                                "Errors.csv")
    errors = []
    with open(errors_file, newline='') as csvfile:
        file_reader = csv.reader(csvfile, delimiter=';')
        data = list(file_reader)        

        counter = 0
        for row in data:
            errors_row = []
            if counter == 0:
                errors_row.append(row[4])
                errors_row.append(row[7])
                errors_row.append(row[8])
                errors_row.append(row[9])
                errors_row.append(row[10])
            else:
                errors_row.append(float(row[4]))
                errors_row.append(float(row[7]))
                errors_row.append(float(row[8]))
                errors_row.append(float(row[9]))
                errors_row.append(float(row[10]))
            errors.append(errors_row)
            counter += 1

    return errors  

def test_errors(errors,
                vem_order,
                tol):
    print(errors)
    num_rows = len(errors)

    if (num_rows == 2):
        assert abs(errors[1][1]) < tol * abs(errors[1][3])
        assert abs(errors[1][2]) < tol * abs(errors[1][4])
    elif (num_rows == 3):
        slope_L2 = abs(math.log(errors[2][1] / errors[2][3])  - math.log(errors[1][1] / errors[2][3])) / \
        abs(math.log(errors[2][0]) - math.log(errors[1][0])) * (-2.0)
        slope_H1 = abs(math.log(errors[2][2] / errors[2][4])  - math.log(errors[1][2] / errors[2][4])) / \
        abs(math.log(errors[2][0]) - math.log(errors[1][0])) * (-2.0)
        print(slope_L2, slope_H1)
        assert slope_L2 < tol * float(vem_order + 1.0)
        assert slope_H1 < tol * float(vem_order)
    else:
        raise Exception("Case {0} not managed".format(num_rows))

if __name__ == "__main__":
    program_folder = os.path.dirname(os.path.realpath(__file__))
    program_path = os.path.join(".",program_folder, "Elliptic_PCC_2D")

    vem_types = [1, 2, 3]
    vem_orders = [1, 2, 3]
    export_folder = "integration_tests"
    tol = 1.0e-8

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
            errors = import_errors(export_path)
            test_errors(errors,
                        vem_order,
                        tol)
            
    test_type = 2
    mesh_generator = 0
    mesh_max_areas = [0.01, 0.001]
    for vem_type in vem_types:
        for vem_order in vem_orders:
            for mesh_max_area in mesh_max_areas:
                run_program(program_folder,
                            program_path,
                            "Run_MG{0}".format(mesh_generator),
                            vem_type,
                            vem_order, 
                            test_type,
                            mesh_generator,
                            mesh_max_area)
                errors = import_errors(export_path)
                test_errors(errors,
                            vem_order,
                            tol)
    
    test_type = 2
    mesh_generator = 2
    mesh_max_areas = [0.01, 0.001]
    for vem_type in vem_types:
        for vem_order in vem_orders:
            for mesh_max_area in mesh_max_areas:
                run_program(program_folder,
                            program_path,
                            "Run_MG{0}".format(mesh_generator),
                            vem_type,
                            vem_order, 
                            test_type,
                            mesh_generator,
                            mesh_max_area)
                errors = import_errors(export_path)
                test_errors(errors,
                            vem_order,
                            tol)
    
    os.system("rm -rf " + os.path.join(program_folder, export_folder))
    print("TESTS SUCCESS")

