import os
import csv

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

        for row in data:
            errors_row = []
            errors_row.append(row[4])
            errors_row.append(row[7])
            errors_row.append(row[8])
            errors_row.append(row[9])
            errors_row.append(row[10])
            errors.append(errors_row)

    return errors  

def test_errors(errors):
    exit()  

if __name__ == "__main__":
    program_folder = os.path.dirname(os.path.realpath(__file__))
    program_path = os.path.join(".",program_folder, "Elliptic_PCC_2D")

    vem_types = [1, 2, 3]
    vem_orders = [1, 2, 3]
    export_folder = "integration_tests"

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
            test_errors(errors)
            
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
    
    os.system("rm -rf " + os.path.join(program_folder, export_folder))

