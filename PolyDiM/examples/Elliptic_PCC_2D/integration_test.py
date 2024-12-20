import os

program_folder = os.path.dirname(os.path.realpath(__file__))
program_path = os.path.join(".",program_folder, "Elliptic_PCC_2D")

test_types = [1, 2]
vem_types = [1, 2, 3]
vem_orders = [1, 2, 3]
export_folder = "Run"

for test_type in test_types:
    for vem_type in vem_types:
        for vem_order in vem_orders:
            export_path = os.path.join(".",
                                    program_folder, 
                                    "{0}_TT{1}_VT{2}".format(
                                        export_folder, 
                                        test_type,
                                        vem_type))
            program_parameters = "VemType:uint={0} VemOrder:uint={1} ExportFolder:string={2} TestType:uint={3}".format(
                vem_type,
                vem_order, 
                export_path,
                test_type)
            os.system(program_path + " " + program_parameters)
