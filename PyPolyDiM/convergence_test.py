import os
import csv
import matplotlib.pyplot as plt
import numpy as np
import shutil


def import_errors_mcc(export_path, method_type, method_order, test_type):
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
                errors_row.append(row[7]) # h
                errors_row.append(row[8]) # error l2 pressure
                errors_row.append(row[9]) # error l2 velocity
                errors_row.append(row[10]) # error l2 super pressure
                errors_row.append(row[11]) # norm l2 pressure
                errors_row.append(row[12]) # norm l2 velocity
            else:
                errors_row.append(float(row[7]))
                errors_row.append(float(row[8]))
                errors_row.append(float(row[9]))
                errors_row.append(float(row[10]))
                errors_row.append(float(row[11]))
                errors_row.append(float(row[12]))
            errors.append(errors_row)
            counter += 1

    return errors

def import_errors_pcc(export_path, method_type, method_order, test_type):
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

def test_errors_pcc(errors, method_order, tol):
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

def test_errors_mcc(errors, method_order, tol):
    num_rows = len(errors)

    if num_rows == 2:
        print('\x1b[0;31;40m' + "Num. Ref. 1: ", abs(errors[1][1]) / abs(errors[1][3]), abs(errors[1][2]) / abs(errors[1][4]),  '\x1b[0m')
        assert abs(errors[1][1]) < tol * abs(errors[1][3])
        assert abs(errors[1][2]) < tol * abs(errors[1][4])
    else:
        errors = np.array(errors[1:])
        slope_l2_press = np.polyfit(np.log(errors[:, 0]), np.log(errors[:, 1]), 1)[0]
        slope_l2_vel = np.polyfit(np.log(errors[:, 0]), np.log(errors[:, 2]), 1)[0]
        slope_l2_press_super = np.polyfit(np.log(errors[:, 0]), np.log(errors[:, 3]), 1)[0]
        print('\x1b[0;31;40m' + "Num. Ref. ", str(num_rows-1), ": ", slope_l2_press, slope_l2_vel, slope_l2_press_super, '\x1b[0m')
        assert round(slope_l2_press) >= round(float(method_order))
        assert round(slope_l2_press_super) >= round(float(method_order + 1))
        assert round(slope_l2_vel) >= round(float(method_order))

dirpath_pcc_2d = "./Export/Export_PCC_2D"
if os.path.exists(dirpath_pcc_2d) and os.path.isdir(dirpath_pcc_2d):
    shutil.rmtree(dirpath_pcc_2d)

dirpath_pcc_3d = "./Export/Export_PCC_3D"
if os.path.exists(dirpath_pcc_3d) and os.path.isdir(dirpath_pcc_3d):
    shutil.rmtree(dirpath_pcc_3d)

dirpath_mcc_2d = "./Export/Export_MCC_2D"
if os.path.exists(dirpath_mcc_2d) and os.path.isdir(dirpath_mcc_2d):
    shutil.rmtree(dirpath_mcc_2d)

execute_test_pcc = False
execute_test_pcc_3d = False
execute_test_mcc = True

# Test PCC 2D
if execute_test_pcc:
    tol = 1.0e-12
    test_type = 1
    method_orders = [1, 2, 3]
    method_types = [0, 1, 2, 3]
    mesh_types = [0]
    for mesh_type in mesh_types:
        for method_type in method_types:
            for order in method_orders:
                export_path = dirpath_pcc_2d + "/Export_" +  str(method_type) + "_" +  str(order) + "_" + str(mesh_type)
                os.system("python ./main_elliptic_pcc_2d.py --method-order={0} --method-type={1} --test-id=1 --mesh-type={3} --mesh-max-relative-area=0.1 --export-path={2}".format(order, method_type, export_path, mesh_type))
                os.system("python ./main_elliptic_pcc_2d.py --method-order={0} --method-type={1} --test-id=1 --mesh-type={3} --mesh-max-relative-area=0.05 --export-path={2}".format(order, method_type, export_path, mesh_type))
                os.system("python ./main_elliptic_pcc_2d.py --method-order={0} --method-type={1} --test-id=1 --mesh-type={3} --mesh-max-relative-area=0.01 --export-path={2}".format(order, method_type, export_path, mesh_type))
                os.system("python ./main_elliptic_pcc_2d.py --method-order={0} --method-type={1} --test-id=1 --mesh-type={3} --mesh-max-relative-area=0.005 --export-path={2}".format(order, method_type, export_path, mesh_type))

                errors = import_errors_pcc(export_path, method_type, order, test_type)
                test_errors_pcc(errors, order, tol)

    tol = 1.0e-12
    test_type = 1
    method_orders = [1]
    method_types = [0, 1, 2, 3]
    mesh_types = [5]
    for mesh_type in mesh_types:
        for method_type in method_types:
            for order in method_orders:
                export_path = dirpath_pcc_2d + "/Export_" +  str(method_type) + "_" +  str(order) + "_" + str(mesh_type)
                os.system("python ./main_elliptic_pcc_2d.py --method-order={0} --method-type={1} --test-id=1 --mesh-type={3} --mesh-max-relative-area=0.1 --export-path={2}".format(order, method_type, export_path, mesh_type))
                os.system("python ./main_elliptic_pcc_2d.py --method-order={0} --method-type={1} --test-id=1 --mesh-type={3} --mesh-max-relative-area=0.05 --export-path={2}".format(order, method_type, export_path, mesh_type))
                os.system("python ./main_elliptic_pcc_2d.py --method-order={0} --method-type={1} --test-id=1 --mesh-type={3} --mesh-max-relative-area=0.01 --export-path={2}".format(order, method_type, export_path, mesh_type))
                os.system("python ./main_elliptic_pcc_2d.py --method-order={0} --method-type={1} --test-id=1 --mesh-type={3} --mesh-max-relative-area=0.005 --export-path={2}".format(order, method_type, export_path, mesh_type))

                errors = import_errors_pcc(export_path, method_type, order, test_type)
                test_errors_pcc(errors, order, tol)

    tol = 1.0e-12
    test_type = 1
    method_orders = [1, 2, 3]
    method_types = [1, 2, 3]
    mesh_types = [2]
    for mesh_type in mesh_types:
        for method_type in method_types:
            for order in method_orders:
                export_path = dirpath_pcc_2d + "/Export_" + str(method_type) + "_" + str(order) + "_" + str(mesh_type)
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

                errors = import_errors_pcc(export_path, method_type, order, test_type)
                test_errors_pcc(errors, order, tol)

# Test MCC 2D
if execute_test_mcc:
    tol = 1.0e-12
    test_type = 1
    method_orders = [0, 1, 2]
    method_types = [1, 3]
    mesh_types = [0]
    for mesh_type in mesh_types:
        for method_type in method_types:
            for order in method_orders:
                export_path = dirpath_mcc_2d + "/Export_" +  str(method_type) + "_" +  str(order) + "_" + str(mesh_type)
                os.system("python ./main_elliptic_mcc_2d.py --method-order={0} --method-type={1} --test-id=1 --mesh-type={3} --mesh-max-relative-area=0.1 --export-path={2}".format(order, method_type, export_path, mesh_type))
                os.system("python ./main_elliptic_mcc_2d.py --method-order={0} --method-type={1} --test-id=1 --mesh-type={3} --mesh-max-relative-area=0.05 --export-path={2}".format(order, method_type, export_path, mesh_type))
                os.system("python ./main_elliptic_mcc_2d.py --method-order={0} --method-type={1} --test-id=1 --mesh-type={3} --mesh-max-relative-area=0.01 --export-path={2}".format(order, method_type, export_path, mesh_type))
                os.system("python ./main_elliptic_mcc_2d.py --method-order={0} --method-type={1} --test-id=1 --mesh-type={3} --mesh-max-relative-area=0.005 --export-path={2}".format(order, method_type, export_path, mesh_type))

                errors = import_errors_mcc(export_path, method_type, order, test_type)
                test_errors_mcc(errors, order, tol)

    tol = 1.0e-12
    test_type = 1
    method_orders = [0, 1, 2]
    method_types = [1, 3]
    mesh_types = [5]
    for mesh_type in mesh_types:
        for method_type in method_types:
            for order in method_orders:
                export_path = dirpath_mcc_2d + "/Export_" +  str(method_type) + "_" +  str(order) + "_" + str(mesh_type)
                os.system("python ./main_elliptic_mcc_2d.py --method-order={0} --method-type={1} --test-id=1 --mesh-type={3} --mesh-max-relative-area=0.1 --export-path={2}".format(order, method_type, export_path, mesh_type))
                os.system("python ./main_elliptic_mcc_2d.py --method-order={0} --method-type={1} --test-id=1 --mesh-type={3} --mesh-max-relative-area=0.05 --export-path={2}".format(order, method_type, export_path, mesh_type))
                os.system("python ./main_elliptic_mcc_2d.py --method-order={0} --method-type={1} --test-id=1 --mesh-type={3} --mesh-max-relative-area=0.01 --export-path={2}".format(order, method_type, export_path, mesh_type))
                os.system("python ./main_elliptic_mcc_2d.py --method-order={0} --method-type={1} --test-id=1 --mesh-type={3} --mesh-max-relative-area=0.005 --export-path={2}".format(order, method_type, export_path, mesh_type))

                errors = import_errors_mcc(export_path, method_type, order, test_type)
                test_errors_mcc(errors, order, tol)

    tol = 1.0e-12
    test_type = 1
    method_orders = [0, 1, 2]
    method_types = [1, 3]
    mesh_types = [2]
    for mesh_type in mesh_types:
        for method_type in method_types:
            for order in method_orders:
                export_path = dirpath_mcc_2d + "/Export_" + str(method_type) + "_" + str(order) + "_" + str(mesh_type)
                os.system(
                    "python ./main_elliptic_mcc_2d.py --method-order={0} --method-type={1} --test-id=1 --mesh-type={3} --mesh-max-relative-area=0.1 --export-path={2}".format(
                        order, method_type, export_path, mesh_type))
                os.system(
                    "python ./main_elliptic_mcc_2d.py --method-order={0} --method-type={1} --test-id=1 --mesh-type={3} --mesh-max-relative-area=0.05 --export-path={2}".format(
                        order, method_type, export_path, mesh_type))
                os.system(
                    "python ./main_elliptic_mcc_2d.py --method-order={0} --method-type={1} --test-id=1 --mesh-type={3} --mesh-max-relative-area=0.01 --export-path={2}".format(
                        order, method_type, export_path, mesh_type))
                os.system(
                    "python ./main_elliptic_mcc_2d.py --method-order={0} --method-type={1} --test-id=1 --mesh-type={3} --mesh-max-relative-area=0.005 --export-path={2}".format(
                        order, method_type, export_path, mesh_type))

                errors = import_errors_mcc(export_path, method_type, order, test_type)
                test_errors_mcc(errors, order, tol)


# Test PCC 3D

if execute_test_pcc_3d:
    tol = 1.0e-12
    test_type = 1
    method_orders = [1, 2, 3]
    method_types = [0, 1, 2, 3]
    mesh_types = [0]
    for mesh_type in mesh_types:
        for method_type in method_types:
            for order in method_orders:
                export_path = dirpath_pcc_3d + "/Export_" +  str(method_type) + "_" +  str(order) + "_" + str(mesh_type)
                os.system("python ./main_elliptic_pcc_3d.py --method-order={0} --method-type={1} --test-id=1 --mesh-type={3} --mesh-max-relative-volume=0.005 --export-path={2}".format(order, method_type, export_path, mesh_type))
                os.system("python ./main_elliptic_pcc_3d.py --method-order={0} --method-type={1} --test-id=1 --mesh-type={3} --mesh-max-relative-volume=0.001 --export-path={2}".format(order, method_type, export_path, mesh_type))

                errors = import_errors_pcc(export_path, method_type, order, test_type)
                test_errors_pcc(errors, order, tol)
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
                export_path = dirpath_pcc_3d + "/Export_" +  str(method_type) + "_" +  str(order) + "_" + str(mesh_type)
                os.system("python ./main_elliptic_pcc_3d.py --method-order={0} --method-type={1} --test-id=1 --mesh-type={3} --mesh-max-relative-volume=0.0025 --export-path={2}".format(order, method_type, export_path, mesh_type))
                os.system("python ./main_elliptic_pcc_3d.py --method-order={0} --method-type={1} --test-id=1 --mesh-type={3} --mesh-max-relative-volume=0.0005 --export-path={2}".format(order, method_type, export_path, mesh_type))

                errors = import_errors_pcc(export_path, method_type, order, test_type)
                test_errors_pcc(errors, order, tol)
