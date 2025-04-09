import os
import csv
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pickle


def draw_loglog_slope(fig, ax, origin, width_inches, slope, inverted=False, color=None, polygon_kwargs=None,
                      label=True, labelcolor=None, label_kwargs=None, zorder=None):
    """
    This function draws slopes or "convergence triangles" into loglog plots.
    @param fig: The figure
    @param ax: The axes object to draw to
    @param origin: The 2D origin (usually lower-left corner) coordinate of the triangle
    @param width_inches: The width in inches of the triangle
    @param slope: The slope of the triangle, i.e. order of convergence
    @param inverted: Whether to mirror the triangle around the origin, i.e. whether
        it indicates the slope towards the lower left instead of upper right (defaults to false)
    @param color: The color of the of the triangle edges (defaults to default color)
    @param polygon_kwargs: Additional kwargs to the Polygon draw call that creates the slope
    @param label: Whether to enable labeling the slope (defaults to true)
    @param labelcolor: The color of the slope labels (defaults to the edge color)
    @param label_kwargs: Additional kwargs to the Annotation draw call that creates the labels
    @param zorder: The z-order value of the triangle and labels, defaults to a high value
    """

    if polygon_kwargs is None:
        polygon_kwargs = {}
    if label_kwargs is None:
        label_kwargs = {}

    if color is not None:
        polygon_kwargs["color"] = color
    if "linewidth" not in polygon_kwargs:
        polygon_kwargs["linewidth"] = 0.75 * mpl.rcParams["lines.linewidth"]
    if labelcolor is not None:
        label_kwargs["color"] = labelcolor
    if "color" not in label_kwargs:
        label_kwargs["color"] = polygon_kwargs["color"]
    if "fontsize" not in label_kwargs:
        label_kwargs["fontsize"] = 3 * mpl.rcParams["font.size"]

    if inverted:
        width_inches = -width_inches
    if zorder is None:
        zorder = 10

    # For more information on coordinate transformations in Matplotlib see
    # https://matplotlib.org/3.1.1/tutorials/advanced/transforms_tutorial.html

    # Convert the origin into figure coordinates in inches
    origin_disp = ax.transData.transform(origin)
    origin_dpi = fig.dpi_scale_trans.inverted().transform(origin_disp)

    # Obtain the bottom-right corner in data coordinates
    corner_dpi = origin_dpi + width_inches * np.array([1.0, 0.0])
    corner_disp = fig.dpi_scale_trans.transform(corner_dpi)
    corner = ax.transData.inverted().transform(corner_disp)

    (x1, y1) = (origin[0], origin[1])
    x2 = corner[0]

    # The width of the triangle in data coordinates
    width = x2 - x1
    # Compute offset of the slope
    log_offset = y1 / (x1 ** slope)

    y2 = log_offset * ((x1 + width) ** slope)
    height = y2 - y1

    # The vertices of the slope
    a = origin
    b = corner
    c = [x2, y2]

    # Draw the slope triangle
    X = np.array([a, b, c])
    triangle = plt.Polygon(X[:3, :], fill=False, zorder=zorder, **polygon_kwargs)
    ax.add_patch(triangle)

    # Convert vertices into display space
    a_disp = ax.transData.transform(a)
    b_disp = ax.transData.transform(b)
    c_disp = ax.transData.transform(c)

    # Figure out the center of the triangle sides in display space
    bottom_center_disp = a_disp + 0.5 * (b_disp - a_disp)
    bottom_center = ax.transData.inverted().transform(bottom_center_disp)

    right_center_disp = b_disp + 0.5 * (c_disp - b_disp)
    right_center = ax.transData.inverted().transform(right_center_disp)

    # Label alignment depending on inversion parameter
    va_xlabel = "top" if not inverted else "bottom"
    ha_ylabel = "left" if not inverted else "right"

    # Label offset depending on inversion parameter
    offset_xlabel = [0.0, -0.33 * label_kwargs["fontsize"]] if not inverted else [0.0, 0.33 * label_kwargs["fontsize"]]
    offset_ylabel = [0.33 * label_kwargs["fontsize"], 0.0] if not inverted else [-0.33 * label_kwargs["fontsize"], 0.0]

    # Draw the slope labels
    # ax.annotate("$1$", bottom_center, xytext=offset_xlabel, textcoords='offset points', ha="center", va=va_xlabel, zorder=zorder, **label_kwargs)
    ax.annotate(f"${slope}$", right_center, xytext=offset_ylabel, textcoords='offset points', ha=ha_ylabel, va="center", zorder=zorder, **label_kwargs)

def run_program(program_folder,
                program_path,
                run_folder,
                vem_type,
                vem_order,
                test_type,
                mesh_generator,
                mesh_max_area,
                num_ref,
                conv_form,
                num_nl_iterations,
                mesh_import_path = "./"):
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
    program_parameters += " NLAbsChangeInSolutionTolerance:double={0}".format(1.0e-9)
    program_parameters += " NLAbsResidualTolerance:double={0}".format(1.0e-10)
    program_parameters += " NLRelChangeInSolutionTolerance:double={0}".format(1.0e-9)
    program_parameters += " NLRelResidualTolerance:double={0}".format(1.0e-10)
    program_parameters += " ConvectiveForm:uint={0}".format(conv_form)
    program_parameters += " MeshImportFilePath:string={0}".format(mesh_import_path)

    output_file = os.path.join(program_folder,
                               "terminal.log")

    run_label = "VemType {0}".format(vem_type)
    run_label += " VemOrder {0}".format(vem_order)
    run_label += " TestType {0}".format(test_type)
    run_label += " ConvectiveForm {0}".format(conv_form)
    run_label += " MeshGenerator {0}".format(mesh_generator)
    run_label += " NumRefinement {0}".format(num_ref)
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
                errors_row.append(col[6]) # h
                errors_row.append(col[7]) # error h1
                errors_row.append(col[8]) # error l2
                errors_row.append(col[9]) # norm h1
                errors_row.append(col[10]) # norm l2
                errors_row.append(col[12]) # nl it
                errors_row.append(col[13]) # residual
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
                test_type,
                conv_form):
    num_rows = len(errors)

    if (num_rows == 2):
        print("Patch: ", abs(errors[1][1]) / abs(errors[1][3]), abs(errors[1][2]) / abs(errors[1][4]), int(errors[1][5]), errors[1][6])
        assert abs(errors[1][1]) < tol * abs(errors[1][3])
        assert abs(errors[1][2]) < tol * abs(errors[1][4])
        if test_type == 1 and conv_form != 1:
            assert abs(errors[1][5]) == 1
        assert abs(errors[1][6]) < tol
    else:
        errors = np.array(errors[1:])
        slope_L2_pres = np.polyfit(np.log(errors[:, 0]), np.log(errors[:, 2]), 1)[0]
        slope_H1_vel = np.polyfit(np.log(errors[:, 0]), np.log(errors[:, 1]), 1)[0]
        print("Num. Ref. ", str(num_rows-1), " : ", slope_H1_vel, slope_L2_pres, errors[:, 5], errors[:, 6])
        assert round(slope_L2_pres) >= round(float(vem_order))
        if test_type == 4 and conv_form == 0:
            assert round(slope_H1_vel) >= round(float(vem_order + 2))
        else:
            assert round(slope_H1_vel) >= round(float(vem_order))


if __name__ == "__main__":
    program_folder = os.path.dirname(os.path.realpath(__file__))
    program_path = os.path.join(".", program_folder, "NavierStokes_DF_PCC_2D")

    remove_folder = False

    vem_type = 1
    vem_orders = [2, 3, 4]
    export_folder = "integration_tests"
    tol = 5.0e-8

    print("RUN TESTS...")

    test_type = 1
    mesh_generator = 1
    mesh_max_area = 0.0
    conv_forms = [0, 2]
    for conv_form in conv_forms:
        for vem_order in vem_orders:
            export_path = run_program(program_folder,
                                      program_path,
                                      "Run_MG{0}_CF{1}".format(mesh_generator, conv_form),
                                      vem_type,
                                      vem_order,
                                      test_type,
                                      mesh_generator,
                                      mesh_max_area,
                                      0,
                                      conv_form,
                                      10)
            errors = import_errors(export_path, vem_type, vem_order, test_type)
            test_errors(errors,
                        vem_order,
                        tol,
                        test_type,
                        conv_form)
            if remove_folder:
                os.system("rm -rf " + os.path.join(program_folder, export_path))

    test_type = 2
    mesh_generator = 5
    mesh_max_areas = [0.125*0.125, 0.0625*0.0625]
    conv_forms = [2]
    for conv_form in conv_forms:
        for vem_order in vem_orders:
            num_ref = 0
            for mesh_max_area in mesh_max_areas:
                export_path = run_program(program_folder,
                                          program_path,
                                          "Run_MG{0}_CF{1}".format(mesh_generator, conv_form),
                                          vem_type,
                                          vem_order,
                                          test_type,
                                          mesh_generator,
                                          mesh_max_area,
                                          num_ref,
                                          conv_form,
                                          10)
                num_ref += 1
            errors = import_errors(export_path, vem_type, vem_order, test_type)
            test_errors(errors,
                        vem_order,
                        tol,
                        test_type,
                        conv_form)
            if remove_folder:
                os.system("rm -rf " + os.path.join(program_folder, export_path))

    test_type = 3
    mesh_generator = 5
    mesh_max_areas = [0.125*0.125, 0.0625*0.0625]
    conv_forms = [0, 1]
    for conv_form in conv_forms:
        for vem_order in vem_orders:
            num_ref = 0
            for mesh_max_area in mesh_max_areas:
                export_path = run_program(program_folder,
                                          program_path,
                                          "Run_MG{0}_CF{1}".format(mesh_generator, conv_form),
                                          vem_type,
                                          vem_order,
                                          test_type,
                                          mesh_generator,
                                          mesh_max_area,
                                          num_ref,
                                          conv_form,
                                          10)
                num_ref += 1
            errors = import_errors(export_path, vem_type, vem_order, test_type)
            test_errors(errors,
                        vem_order,
                        tol,
                        test_type,
                        conv_form)
            if remove_folder:
                os.system("rm -rf " + os.path.join(program_folder, export_path))

    test_type = 3
    mesh_generator = 2
    mesh_max_areas = [0.01, 0.005]
    conv_forms = [0, 1]
    for conv_form in conv_forms:
        for vem_order in vem_orders:
            num_ref = 0
            for mesh_max_area in mesh_max_areas:
                export_path = run_program(program_folder,
                                          program_path,
                                          "Run_MG{0}_CF{1}".format(mesh_generator, conv_form),
                                          vem_type,
                                          vem_order,
                                          test_type,
                                          mesh_generator,
                                          mesh_max_area,
                                          num_ref,
                                          conv_form,
                                          10)
                num_ref += 1

            errors = import_errors(export_path, vem_type, vem_order, test_type)
            test_errors(errors,
                        vem_order,
                        tol,
                        test_type,
                        conv_form)
            if remove_folder:
                os.system("rm -rf " + os.path.join(program_folder, export_path))

    test_type = 4
    mesh_generator = 4
    mesh_import_paths = ["../../../../Mesh/2D/CircleTriangularMesh/CircleTriangularMesh_R1",
                         "../../../../Mesh/2D/CircleTriangularMesh/CircleTriangularMesh_R2",
                         "../../../../Mesh/2D/CircleTriangularMesh/CircleTriangularMesh_R3"]
    conv_forms = [0, 1]
    vem_orders = [2]
    list_errors = []
    for conv_form in conv_forms:
        for vem_order in vem_orders:
            num_ref = 0
            for mesh_imp_path in mesh_import_paths:
                export_path = run_program(program_folder,
                                          program_path,
                                          "Run_MG{0}_CF{1}".format(mesh_generator, conv_form),
                                          vem_type,
                                          vem_order,
                                          test_type,
                                          mesh_generator,
                                          0.0,
                                          num_ref,
                                          conv_form,
                                          10,
                                          mesh_imp_path)
                num_ref += 1

            errors = import_errors(export_path, vem_type, vem_order, test_type)
            test_errors(errors,
                        vem_order,
                        tol,
                        test_type,
                        conv_form)

            list_errors.append(np.array(errors[1:]))

            if remove_folder:
                os.system("rm -rf " + os.path.join(program_folder, export_path))

    fig, ax = plt.subplots(figsize=(12, 12))
    for h in range(len(conv_forms)):
        errors = list_errors[h]
        num_rows = len(errors)

        if h == 0:
            ax.plot(errors[:, 0], errors[:, 1], '-k^', linewidth=2, markersize=20,
                    label="non-skew")
        elif h == 1:
            ax.plot(errors[:, 0], errors[:, 1], '-ro', linewidth=2, markersize=20,
                    label="skew")
        else:
            raise ValueError("Not valid method type")

    plt.legend(bbox_to_anchor=(0., 1.02, 1.0, 0.2), loc="lower left",
               mode="expand", borderaxespad=0, ncol=3, fontsize=30)

    plt.xlabel('$h$', fontsize=30)
    plt.ylabel('$e_{\\nabla \\mathbf{u}}$', fontsize=30)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.yscale('log')
    plt.xscale('log')
    plt.grid(True, which="both", ls="--")
    t_pos = [0.1, 1.0e-03, 0.1, 5.0e-08]
    draw_loglog_slope(fig, ax, (t_pos[0], t_pos[1]), 3, slope=2, color='k')
    draw_loglog_slope(fig, ax, (t_pos[2], t_pos[3]), 3, slope=4, color='k')
    plt.savefig(export_folder + "/{}_h1_velocity_decay_plot.png".format(test_type), bbox_inches='tight', dpi=300)
    plt.show()

    fig, ax = plt.subplots(figsize=(12, 12))
    for h in range(len(conv_forms)):
        errors = list_errors[h]
        num_rows = len(errors)

        if h == 0:
            ax.plot(errors[:, 0], errors[:, 2], '-k^', linewidth=2, markersize=20,
                    label="non-skew")
        elif h == 1:
            ax.plot(errors[:, 0], errors[:, 2], '-ro', linewidth=2, markersize=20,
                    label="skew")
        else:
            raise ValueError("Not valid method type")

    plt.legend(bbox_to_anchor=(0., 1.02, 1.0, 0.2), loc="lower left",
               mode="expand", borderaxespad=0, ncol=3, fontsize=30)

    plt.xlabel('$h$', fontsize=30)
    plt.ylabel('$e_{p}$', fontsize=30)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.yscale('log')
    plt.xscale('log')
    plt.grid(True, which="both", ls="--")
    t_pos = [0.1, 0.008]
    draw_loglog_slope(fig, ax, (t_pos[0], t_pos[1]), 3, slope=2, color='k')
    plt.savefig(export_folder + "/{}_l2_pressure_decay_plot.png".format(test_type), bbox_inches='tight', dpi=300)
    plt.show()

    if remove_folder:
        os.system("rm -rf " + os.path.join(program_folder, export_folder))

    print("TESTS SUCCESS")

