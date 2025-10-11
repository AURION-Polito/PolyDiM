import numpy as np
from pypolydim import polydim, gedim
from Elliptic_PCC_2D.test_definition import ITest
import scipy.sparse.linalg as sla
from pypolydim.assembler_utilities import assembler_utilities
from scipy.sparse import coo_array


class Assembler:

    class EllipticPCC2DProblemData:

        def __init__(self):
            self.global_matrix_a_data = assembler_utilities.SparseMatrix()
            self.dirichlet_matrix_a_data = assembler_utilities.SparseMatrix()
            self.global_matrix_a: coo_array
            self.dirichlet_matrix_a: coo_array
            self.right_hand_side: np.ndarray = np.ndarray(0)
            self.solution: np.ndarray = np.ndarray(0)
            self.solution_dirichlet: np.ndarray = np.ndarray(0)


    class PostProcessData:

        def __init__(self):
            self.cell0_ds_numeric: np.ndarray
            self.cell0_ds_exact: np.ndarray

            self.cell2_ds_error_l2: np.ndarray
            self.cell2_ds_norm_l2: np.ndarray
            self.error_l2: float = 0.0
            self.norm_l2: float = 0.0
            self.cell2_ds_error_h1: np.ndarray
            self.cell2_ds_norm_h1: np.ndarray
            self.error_h1: float = 0.0
            self.norm_h1: float = 0.0

            self.mesh_size: float = np.inf
            self.residual_norm: float = np.inf

    @staticmethod
    def compute_strong_term(cell2_d_index: int,
                            mesh: gedim.MeshMatricesDAO,
                            mesh_do_fs_info: polydim.pde_tools.do_fs.DOFsManager.MeshDOFsInfo,
                            do_fs_data: polydim.pde_tools.do_fs.DOFsManager.DOFsData,
                            reference_element_data: polydim.pde_tools.local_space_pcc_2_d.ReferenceElement_Data,
                            local_space_data: polydim.pde_tools.local_space_pcc_2_d.LocalSpace_Data,
                            test: ITest,
                            assembler_data: EllipticPCC2DProblemData) -> None:

        # Assemble strong boundary condition on Cell0Ds
        for v in range(mesh.cell2_d_number_vertices(cell2_d_index)):

            cell0_d_index = mesh.cell2_d_vertex(cell2_d_index, v)
            boundary_info = mesh_do_fs_info.cells_boundary_info[0][cell0_d_index]

            if boundary_info.type != polydim.pde_tools.do_fs.DOFsManager.MeshDOFsInfo.BoundaryInfo.BoundaryTypes.strong:
                continue

            coordinates = np.expand_dims(mesh.cell0_d_coordinates(cell0_d_index), axis=1)

            strong_boundary_values = test.strong_boundary_condition(boundary_info.marker, coordinates)

            local_dofs = do_fs_data.cells_do_fs[0][cell0_d_index]

            assert len(local_dofs) == len(strong_boundary_values)

            for loc_i in range(len(local_dofs)):

                local_dof_i = local_dofs[loc_i]

                match local_dof_i.type:
                    case polydim.pde_tools.do_fs.DOFsManager.DOFsData.DOF.Types.strong:
                        assembler_data.solution_dirichlet[local_dof_i.global_index] = strong_boundary_values[loc_i]
                        pass
                    case polydim.pde_tools.do_fs.DOFsManager.DOFsData.DOF.Types.dof:
                        pass
                    case _:
                        raise ValueError("Unknown DOF Type")

        # Assemble strong boundary condition on Cell1Ds
        for ed in range(mesh.cell2_d_number_edges(cell2_d_index)):

            cell1_d_index = mesh.cell2_d_edge(cell2_d_index, ed)

            boundary_info = mesh_do_fs_info.cells_boundary_info[1][cell1_d_index]
            local_dofs = do_fs_data.cells_do_fs[1][cell1_d_index]

            if (boundary_info.type != polydim.pde_tools.do_fs.DOFsManager.MeshDOFsInfo.BoundaryInfo.BoundaryTypes.strong
                    or len(local_dofs) == 0):
                continue

            edge_do_fs_coordinates = polydim.pde_tools.local_space_pcc_2_d.edge_dofs_coordinates(reference_element_data, local_space_data, ed)

            strong_boundary_values = test.strong_boundary_condition(boundary_info.marker, edge_do_fs_coordinates)

            assert len(local_dofs) == len(strong_boundary_values)

            for loc_i in range(len(local_dofs)):

                local_dof_i = local_dofs[loc_i]

                match local_dof_i.type:
                    case polydim.pde_tools.do_fs.DOFsManager.DOFsData.DOF.Types.strong:
                        assembler_data.solution_dirichlet[local_dof_i.global_index] = strong_boundary_values[loc_i]
                        pass
                    case polydim.pde_tools.do_fs.DOFsManager.DOFsData.DOF.Types.dof:
                        pass
                    case _:
                        raise ValueError("Unknown DOF Type")

    def compute_weak_term(self,
                          cell2_d_index: int,
                          mesh: gedim.MeshMatricesDAO,
                          mesh_geometric_data: gedim.MeshUtilities.MeshGeometricData2D,
                          mesh_do_fs_info: polydim.pde_tools.do_fs.DOFsManager.MeshDOFsInfo,
                          do_fs_data: polydim.pde_tools.do_fs.DOFsManager.DOFsData,
                          reference_element_data: polydim.pde_tools.local_space_pcc_2_d.ReferenceElement_Data,
                          local_space_data: polydim.pde_tools.local_space_pcc_2_d.LocalSpace_Data,
                          test: ITest,
                          assembler_data: EllipticPCC2DProblemData) -> None:

        num_vertices = mesh.cell2_d_number_vertices(cell2_d_index)

        for ed in range(num_vertices):
            cell1_d_index = mesh.cell2_d_edge(cell2_d_index, ed)

            boundary_info = mesh_do_fs_info.cells_boundary_info[1][cell1_d_index]

            if boundary_info.type != polydim.pde_tools.do_fs.DOFsManager.MeshDOFsInfo.BoundaryInfo.BoundaryTypes.weak:
                continue

            # compute vem values
            weak_reference_segment = gedim.quadrature.Quadrature_Gauss1D.fill_points_and_weights(2 * reference_element_data.order)
            points_curvilinear_coordinates = weak_reference_segment.points[0, :]

            # map edge internal quadrature points
            edge_start = mesh_geometric_data.cell2_ds_vertices[cell2_d_index][:, ed] if mesh_geometric_data.cell2_ds_edge_directions[cell2_d_index][ed] else mesh_geometric_data.cell2_ds_vertices[cell2_d_index][:, (ed + 1) % num_vertices]
            edge_tangent = mesh_geometric_data.cell2_ds_edge_tangents[cell2_d_index][:, ed]
            direction = 1.0 if mesh_geometric_data.cell2_ds_edge_directions[cell2_d_index][ed] else -1.0
            num_edge_weak_quadrature_points = weak_reference_segment.points.shape[1]

            weak_quadrature_points = np.zeros([3, num_edge_weak_quadrature_points])
            for q in range(num_edge_weak_quadrature_points):
                weak_quadrature_points[:, q] = edge_start + direction * weak_reference_segment.points[0, q] * edge_tangent

            edge_length = mesh_geometric_data.cell2_ds_edge_lengths[cell2_d_index][ed]
            weak_quadrature_weights = weak_reference_segment.weights * edge_length

            neumann_values = test.weak_boundary_condition(boundary_info.marker, weak_quadrature_points)
            weak_basis_function_values = (polydim.pde_tools.
                                          local_space_pcc_2_d.
                                          basis_functions_values_on_edge(ed, reference_element_data,
                                                                         local_space_data,
                                                                         points_curvilinear_coordinates))

            # compute values of Neumann
            neumann_contributions = (weak_basis_function_values.T @
                                     np.diag(weak_quadrature_weights) @
                                     neumann_values)

            for p in range(2):
                cell0_d_index = mesh.cell1_d_vertex(cell1_d_index, p)
                local_do_fs = do_fs_data.cells_do_fs[0][cell0_d_index]

                for loc_i in range(len(local_do_fs)):
                    local_dof_i = local_do_fs[loc_i]

                    match local_dof_i.type:
                        case polydim.pde_tools.do_fs.DOFsManager.DOFsData.DOF.Types.strong:
                            continue
                        case polydim.pde_tools.do_fs.DOFsManager.DOFsData.DOF.Types.dof:
                            assembler_data.right_hand_side[local_dof_i.global_index] += neumann_contributions[p]
                            pass
                        case _:
                            raise ValueError("Unknown DOF Type")

            local_do_fs = do_fs_data.cells_do_fs[1][cell1_d_index]
            for loc_i in range(len(local_do_fs)):
                local_dof_i = local_do_fs[loc_i]

                match local_dof_i.type:
                    case polydim.pde_tools.do_fs.DOFsManager.DOFsData.DOF.Types.strong:
                        continue
                    case polydim.pde_tools.do_fs.DOFsManager.DOFsData.DOF.Types.dof:
                        assembler_data.right_hand_side[local_dof_i.global_index] += neumann_contributions[loc_i + 2]
                        pass
                    case _:
                        raise ValueError("Unknown DOF Type")

    @staticmethod
    def solve(do_fs_data: polydim.pde_tools.do_fs.DOFsManager.DOFsData,
              assembler_data: EllipticPCC2DProblemData) -> None:

        if do_fs_data.number_strongs > 0:
            assembler_data.right_hand_side -= assembler_data.dirichlet_matrix_a @ assembler_data.solution_dirichlet

        assembler_data.solution = sla.spsolve(assembler_data.global_matrix_a.tocsc(), assembler_data.right_hand_side)


    def assemble(self,
                 geometry_utilities_confi: gedim.GeometryUtilitiesConfig,
                 mesh: gedim.MeshMatricesDAO,
                 mesh_geometric_data: gedim.MeshUtilities.MeshGeometricData2D,
                 mesh_do_fs_info: polydim.pde_tools.do_fs.DOFsManager.MeshDOFsInfo,
                 do_fs_data: polydim.pde_tools.do_fs.DOFsManager.DOFsData,
                 do_fs_data_indices: polydim.pde_tools.do_fs.DOFsManager.CellsDOFsIndicesData,
                 reference_element_data: polydim.pde_tools.local_space_pcc_2_d.ReferenceElement_Data,
                 test: ITest) -> EllipticPCC2DProblemData:

        result = self.EllipticPCC2DProblemData()

        result.right_hand_side = np.zeros([do_fs_data.number_do_fs])
        result.solution_dirichlet = np.zeros([do_fs_data.number_strongs])

        for c in range(mesh.cell2_d_total_number()):

            local_space_data = polydim.pde_tools.local_space_pcc_2_d.create_local_space(geometry_utilities_confi.tolerance1_d,
                                                                                        geometry_utilities_confi.tolerance2_d,
                                                                                        mesh_geometric_data,
                                                                                        c,
                                                                                        reference_element_data)

            basis_functions_values = polydim.pde_tools.local_space_pcc_2_d.basis_functions_values(reference_element_data, local_space_data)
            basis_functions_derivative_values = polydim.pde_tools.local_space_pcc_2_d.basis_functions_derivative_values(reference_element_data, local_space_data)

            cell2_d_internal_quadrature = polydim.pde_tools.local_space_pcc_2_d.internal_quadrature(reference_element_data, local_space_data)

            diffusion_term_values = test.diffusion_term(cell2_d_internal_quadrature.points)
            source_term_values = test.source_term(cell2_d_internal_quadrature.points)
            weights = cell2_d_internal_quadrature.weights


            equation = polydim.pde_tools.equations.EllipticEquation()
            local_a = equation.compute_cell_diffusion_matrix(diffusion_term_values,
                                                             basis_functions_derivative_values,
                                                             weights)

            local_rhs = equation.compute_cell_forcing_term(source_term_values,
                                                           basis_functions_values,
                                                           weights)

            k_max = np.max(abs(diffusion_term_values))
            local_a_stab = k_max * polydim.pde_tools.local_space_pcc_2_d.stabilization_matrix(reference_element_data, local_space_data)

            global_do_fs = do_fs_data.cells_global_do_fs[2][c]

            assert polydim.pde_tools.local_space_pcc_2_d.size(reference_element_data, local_space_data) == len(global_do_fs)

            local_to_global_data = assembler_utilities.LocalMatrixToGlobalMatrixDOFsData()
            local_to_global_data.do_fs_data_indices = [do_fs_data_indices]
            local_to_global_data.local_offsets = [0]
            local_to_global_data.global_offsets_do_fs = [0]
            local_to_global_data.global_offsets_strongs = [0]

            assembler_utilities.assemble_local_matrix_to_global_matrix(2,
                                                                       c,
                                                                       local_to_global_data,
                                                                       local_to_global_data,
                                                                       local_a + local_a_stab,
                                                                       result.global_matrix_a_data,
                                                                       result.dirichlet_matrix_a_data,
                                                                       local_rhs,
                                                                       result.right_hand_side)


            self.compute_strong_term(c, mesh, mesh_do_fs_info, do_fs_data, reference_element_data, local_space_data,
                                     test, result)
            self.compute_weak_term(c, mesh, mesh_geometric_data, mesh_do_fs_info, do_fs_data, reference_element_data,
                                   local_space_data, test, result)

        result.global_matrix_a = result.global_matrix_a_data.create(do_fs_data.number_do_fs, do_fs_data.number_do_fs)
        result.dirichlet_matrix_a = result.dirichlet_matrix_a_data.create(do_fs_data.number_do_fs, do_fs_data.number_strongs)

        return result


    def post_process_solution(self,
                              geometry_utilities_confi: gedim.GeometryUtilitiesConfig,
                              mesh: gedim.MeshMatricesDAO,
                              mesh_geometric_data: gedim.MeshUtilities.MeshGeometricData2D,
                              do_fs_data: polydim.pde_tools.do_fs.DOFsManager.DOFsData,
                              do_fs_data_indices: polydim.pde_tools.do_fs.DOFsManager.CellsDOFsIndicesData,
                              count_do_fs_data: assembler_utilities.CountDOFsData,
                              reference_element_data: polydim.pde_tools.local_space_pcc_2_d.ReferenceElement_Data,
                              assembler_data: EllipticPCC2DProblemData,
                              test: ITest) -> PostProcessData:

        result = self.PostProcessData()

        result.residual_norm = 0.0
        if do_fs_data.number_do_fs > 0:
            residual = assembler_data.global_matrix_a @ assembler_data.solution - assembler_data.right_hand_side
            result.residual_norm = np.linalg.norm(residual)


        result.cell0_ds_numeric = np.zeros([mesh.cell0_d_total_number()])
        result.cell0_ds_exact = np.zeros([mesh.cell0_d_total_number()])

        for p in range(mesh.cell0_d_total_number()):

            result.cell0_ds_exact[p] = test.exact_solution(np.expand_dims(mesh.cell0_d_coordinates(p), axis=1))[0]

            local_do_fs = do_fs_data.cells_do_fs[0][p]

            for loc_i in range(len(local_do_fs)):

                local_dof_i = local_do_fs[loc_i]

                match local_dof_i.type:
                    case polydim.pde_tools.do_fs.DOFsManager.DOFsData.DOF.Types.strong:
                        result.cell0_ds_numeric[p] = assembler_data.solution_dirichlet[local_dof_i.global_index]
                        pass
                    case polydim.pde_tools.do_fs.DOFsManager.DOFsData.DOF.Types.dof:
                        result.cell0_ds_numeric[p] = assembler_data.solution[local_dof_i.global_index]
                    case _:
                        raise ValueError("Unknown DOF Type")


        result.cell2_ds_error_l2 = np.zeros([mesh.cell2_d_total_number()])
        result.cell2_ds_norm_l2 = np.zeros([mesh.cell2_d_total_number()])
        result.cell2_ds_error_h1 = np.zeros([mesh.cell2_d_total_number()])
        result.cell2_ds_norm_h1 = np.zeros([mesh.cell2_d_total_number()])
        result.error_L2 = 0.0
        result.norm_L2 = 0.0
        result.error_H1 = 0.0
        result.norm_H1 = 0.0
        result.mesh_size = 0.0

        for c in range(mesh.cell2_d_total_number()):

            local_space_data=polydim.pde_tools.local_space_pcc_2_d.create_local_space(
                geometry_utilities_confi.tolerance1_d,
                geometry_utilities_confi.tolerance2_d,
                mesh_geometric_data,
                c,
                reference_element_data)

            basis_functions_values = polydim.pde_tools.local_space_pcc_2_d.basis_functions_values(
                reference_element_data, local_space_data, polydim.vem.pcc.ProjectionTypes.pi0k)
            basis_functions_derivative_values = polydim.pde_tools.local_space_pcc_2_d.basis_functions_derivative_values(
                reference_element_data, local_space_data)

            cell2_d_internal_quadrature = polydim.pde_tools.local_space_pcc_2_d.internal_quadrature(reference_element_data, local_space_data)


            exact_solution_values = test.exact_solution(cell2_d_internal_quadrature.points)
            exact_derivative_solution_values = test.exact_derivative_solution(cell2_d_internal_quadrature.points)

            assembler_utilities_obj = assembler_utilities()
            local_count_do_fs = assembler_utilities_obj.local_count_do_fs(2, c, [do_fs_data])

            do_fs_values = assembler_utilities_obj.global_solution_to_local_solution(2,
                                                                                     c,
                                                                                     [do_fs_data_indices],
                                                                                     count_do_fs_data,
                                                                                     local_count_do_fs,
                                                                                     assembler_data.solution,
                                                                                     assembler_data.solution_dirichlet)

            local_error_l2 = (basis_functions_values @ do_fs_values - exact_solution_values)**2
            local_norm_l2 = (basis_functions_values @ do_fs_values)**2

            result.cell2_ds_error_l2[c] = np.sum(cell2_d_internal_quadrature.weights * local_error_l2)
            result.cell2_ds_norm_l2[c] = np.sum(cell2_d_internal_quadrature.weights * local_norm_l2)

            local_error_h1 = ((basis_functions_derivative_values[0] @ do_fs_values - exact_derivative_solution_values[0])**2
                              + (basis_functions_derivative_values[1] @ do_fs_values - exact_derivative_solution_values[1])**2)

            local_norm_h1 = ((basis_functions_derivative_values[0] @ do_fs_values)**2 +
                             (basis_functions_derivative_values[1] @ do_fs_values)**2)

            result.cell2_ds_error_h1[c] = np.sum(cell2_d_internal_quadrature.weights * local_error_h1)
            result.cell2_ds_norm_h1[c] = np.sum(cell2_d_internal_quadrature.weights * local_norm_h1)

            if mesh_geometric_data.cell2_ds_diameters[c] > result.mesh_size:
                result.mesh_size = mesh_geometric_data.cell2_ds_diameters[c]


        result.error_l2 = np.sqrt(np.sum(result.cell2_ds_error_l2))
        result.norm_l2 = np.sqrt(np.sum(result.cell2_ds_norm_l2))
        result.error_h1 = np.sqrt(np.sum(result.cell2_ds_error_h1))
        result.norm_h1 = np.sqrt(np.sum(result.cell2_ds_norm_h1))

        return result

