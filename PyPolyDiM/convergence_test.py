import os
import csv
import matplotlib.pyplot as plt
import numpy as np
import shutil


def import_errors(export_path, method_type, method_order, test_type):
    errors_file = os.path.join(export_path,
                               "Errors_" + str(test_type) + "_" + str(method_type) + "_" + str(method_order) + ".csv")
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
            else:
                errors_row.append(float(row[7]))
                errors_row.append(float(row[8]))
                errors_row.append(float(row[9]))
                errors_row.append(float(row[10]))
                errors_row.append(float(row[11]))
            errors.append(errors_row)
            counter += 1

    return errors

def test_errors(errors,
                method_order,
                tol):
    num_rows = len(errors)

    if num_rows == 2:
        print('\x1b[0;31;40m' + "Num. Ref. 1: ", abs(errors[1][1]) / abs(errors[1][3]), abs(errors[1][2]) / abs(errors[1][4]),  '\x1b[0m')
        assert abs(errors[1][1]) < tol * abs(errors[1][3])
        assert abs(errors[1][2]) < tol * abs(errors[1][4])
    else:
        errors = np.array(errors[1:])
        slope_l2 = np.polyfit(np.log(errors[:, 0]), np.log(errors[:, 1]), 1)[0]
        slope_h1 = np.polyfit(np.log(errors[:, 0]), np.log(errors[:, 2]), 1)[0]
        print('\x1b[0;31;40m' + "Num. Ref. ", str(num_rows-1), ": ", slope_l2, slope_h1,  '\x1b[0m')
        assert round(slope_l2) >= round(float(method_order + 1.0))
        assert round(slope_h1) >= round(float(method_order))

dirpath = "./Export/Export2D"
if os.path.exists(dirpath) and os.path.isdir(dirpath):
    shutil.rmtree(dirpath)

dirpath = "./Export/Export3D"
if os.path.exists(dirpath) and os.path.isdir(dirpath):
    shutil.rmtree(dirpath)

tol = 1.0e-12
test_type = 1
method_orders = [1, 2, 3]
method_types = [0, 1, 2, 3]
mesh_types = [0]
for mesh_type in mesh_types:
    for method_type in method_types:
        for order in method_orders:
            export_path = "./Export/Export2D/Export_" +  str(method_type) + "_" +  str(order) + "_" + str(mesh_type)
            os.system("python ./main_elliptic_pcc_2d.py --method-order={0} --method-type={1} --test-id=1 --mesh-type={3} --mesh-max-relative-area=0.1 --export-path={2}".format(order, method_type, export_path, mesh_type))
            os.system("python ./main_elliptic_pcc_2d.py --method-order={0} --method-type={1} --test-id=1 --mesh-type={3} --mesh-max-relative-area=0.05 --export-path={2}".format(order, method_type, export_path, mesh_type))
            os.system("python ./main_elliptic_pcc_2d.py --method-order={0} --method-type={1} --test-id=1 --mesh-type={3} --mesh-max-relative-area=0.01 --export-path={2}".format(order, method_type, export_path, mesh_type))
            os.system("python ./main_elliptic_pcc_2d.py --method-order={0} --method-type={1} --test-id=1 --mesh-type={3} --mesh-max-relative-area=0.005 --export-path={2}".format(order, method_type, export_path, mesh_type))

            errors = import_errors(export_path, method_type, order, test_type)
            test_errors(errors, order, tol)

tol = 1.0e-12
test_type = 1
method_orders = [1]
method_types = [0, 1, 2, 3]
mesh_types = [5]
for mesh_type in mesh_types:
    for method_type in method_types:
        for order in method_orders:
            export_path = "./Export/Export2D/Export_" +  str(method_type) + "_" +  str(order) + "_" + str(mesh_type)
            os.system("python ./main_elliptic_pcc_2d.py --method-order={0} --method-type={1} --test-id=1 --mesh-type={3} --mesh-max-relative-area=0.1 --export-path={2}".format(order, method_type, export_path, mesh_type))
            os.system("python ./main_elliptic_pcc_2d.py --method-order={0} --method-type={1} --test-id=1 --mesh-type={3} --mesh-max-relative-area=0.05 --export-path={2}".format(order, method_type, export_path, mesh_type))
            os.system("python ./main_elliptic_pcc_2d.py --method-order={0} --method-type={1} --test-id=1 --mesh-type={3} --mesh-max-relative-area=0.01 --export-path={2}".format(order, method_type, export_path, mesh_type))
            os.system("python ./main_elliptic_pcc_2d.py --method-order={0} --method-type={1} --test-id=1 --mesh-type={3} --mesh-max-relative-area=0.005 --export-path={2}".format(order, method_type, export_path, mesh_type))

            errors = import_errors(export_path, method_type, order, test_type)
            test_errors(errors, order, tol)

tol = 1.0e-12
test_type = 1
method_orders = [1, 2, 3]
method_types = [1, 2, 3]
mesh_types = [2]
for mesh_type in mesh_types:
    for method_type in method_types:
        for order in method_orders:
            export_path = "./Export/Export2D/Export_" + str(method_type) + "_" + str(order) + "_" + str(mesh_type)
            os.system(
                "python ./main_elliptic_pcc_2d.py --method-order={0} --method-type={1} --test-id=1 --mesh-type={3} --mesh-max-relative-area=0.1 --export-path={2}".format(
                    order, method_type, export_path, mesh_type))
            os.system(
                "python ./main_elliptic_pcc_2d.py --method-order={0} --method-type={1} --test-id=1 --mesh-type={3} --mesh-max-relative-area=0.05 --export-path={2}".format(
                    order, method_type, export_path, mesh_type))
            os.system(
                "python ./main_elliptic_pcc_2d.py --method-order={0} --method-type={1} --test-id=1 --mesh-type={3} --mesh-max-relative-area=0.01 --export-path={2}".format(
                    order, method_type, export_path, mesh_type))
            os.system(
                "python ./main_elliptic_pcc_2d.py --method-order={0} --method-type={1} --test-id=1 --mesh-type={3} --mesh-max-relative-area=0.005 --export-path={2}".format(
                    order, method_type, export_path, mesh_type))

            errors = import_errors(export_path, method_type, order, test_type)
            test_errors(errors, order, tol)


tol = 1.0e-12
test_type = 1
method_orders = [1, 2, 3]
method_types = [0, 1, 2, 3]
mesh_types = [0]
for mesh_type in mesh_types:
    for method_type in method_types:
        for order in method_orders:
            export_path = "./Export/Export3D/Export_" +  str(method_type) + "_" +  str(order) + "_" + str(mesh_type)
            os.system("python ./main_elliptic_pcc_3d.py --method-order={0} --method-type={1} --test-id=1 --mesh-type={3} --mesh-max-relative-volume=0.005 --export-path={2}".format(order, method_type, export_path, mesh_type))
            os.system("python ./main_elliptic_pcc_3d.py --method-order={0} --method-type={1} --test-id=1 --mesh-type={3} --mesh-max-relative-volume=0.001 --export-path={2}".format(order, method_type, export_path, mesh_type))

            errors = import_errors(export_path, method_type, order, test_type)
            test_errors(errors, order, tol)
            # errors = np.array(errors[1:])
            # fig, ax = plt.subplots(figsize=(12, 12))
            # ax.plot(errors[:, 0], errors[:, 1], '-k^', linewidth=2, markersize=12)
            # plt.xlabel('$h$', fontsize=30)
            # plt.ylabel('$e_0$', fontsize=30)
            # plt.xticks(fontsize=20)
            # plt.yticks(fontsize=20)
            # plt.yscale('log')
            # plt.xscale('log')
            # plt.grid(True, which="both", ls="--")
            # plt.ylim(None, 10)
            # plt.show()

tol = 1.0e-12
test_type = 1
method_orders = [1, 2, 3]
method_types = [0, 1, 2, 3]
mesh_types = [6]
for mesh_type in mesh_types:
    for method_type in method_types:
        for order in method_orders:
            export_path = "./Export/Export3D/Export_" +  str(method_type) + "_" +  str(order) + "_" + str(mesh_type)
            os.system("python ./main_elliptic_pcc_3d.py --method-order={0} --method-type={1} --test-id=1 --mesh-type={3} --mesh-max-relative-volume=0.0025 --export-path={2}".format(order, method_type, export_path, mesh_type))
            os.system("python ./main_elliptic_pcc_3d.py --method-order={0} --method-type={1} --test-id=1 --mesh-type={3} --mesh-max-relative-volume=0.0005 --export-path={2}".format(order, method_type, export_path, mesh_type))

            errors = import_errors(export_path, method_type, order, test_type)
            test_errors(errors, order, tol)
