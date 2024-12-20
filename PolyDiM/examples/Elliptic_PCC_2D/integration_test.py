import os

program_folder = os.path.dirname(os.path.realpath(__file__))
program_path = os.path.join(".",program_folder, "Elliptic_PCC_2D")

vem_types = [1, 2, 3]
vem_orders = [1, 2, 3]
export_folder = "Run"

for vem_type in vem_types:
    for vem_order in vem_orders:
        export_path = os.path.join(".",
                                program_folder, 
                                "{0}_VT{1}".format(export_folder, vem_type))
        program_parameters = "VemOrder:uint={0} ExportFolder:string={1}".format(vem_order, 
                                                                                export_path)
        os.system(program_path + " " + program_parameters)
