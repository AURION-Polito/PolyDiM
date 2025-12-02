import numpy as np
from pypolydim import polydim, gedim
from Elliptic_MCC_2D.test_definition import ITest
import scipy.sparse.linalg as sla
from pypolydim.assembler_utilities import assembler_utilities
from scipy.sparse import coo_array
from typing import List


class Assembler:

    class EllipticMCC2DProblemData:

        def __init__(self):
            self.global_matrix_a_data = assembler_utilities.SparseMatrix()
            self.neumann_matrix_a_data = assembler_utilities.SparseMatrix()
            self.global_matrix_a: coo_array
            self.neumann_matrix_a: coo_array
            self.right_hand_side: np.ndarray = np.ndarray(0)
            self.solution: np.ndarray = np.ndarray(0)
            self.solution_neumann: np.ndarray = np.ndarray(0)


    class PostProcessData:

        def __init__(self):
            self.cell2_ds_numeric_pressure: np.ndarray = np.ndarray(0)
            self.cell2_ds_exact_pressure: np.ndarray = np.ndarray(0)

            self.cell2_ds_error_l2_pressure: np.ndarray = np.ndarray(0)
            self.cell2_ds_norm_l2_pressure: np.ndarray = np.ndarray(0)
            self.error_l2_pressure: float = 0.0
            self.norm_l2_pressure: float = 0.0

            self.cell2_ds_super_error_l2_pressure: np.ndarray = np.ndarray(0)
            self.super_error_l2_pressure: float = 0.0

            self.cell2_ds_error_l2_velocity: np.ndarray = np.ndarray(0)
            self.cell2_ds_norm_l2_velocity: np.ndarray = np.ndarray(0)
            self.error_l2_velocity: float = 0.0
            self.norm_l2_velocity: float = 0.0

            self.mesh_size: float = 0.0
            self.residual_norm: float = 0.0

    @staticmethod
    def compute_strong_term(cell2_d_index: int,
                            mesh: gedim.MeshMatricesDAO,
                            mesh_do_fs_info: polydim.pde_tools.do_fs.DOFsManager.MeshDOFsInfo,
                            do_fs_data: polydim.pde_tools.do_fs.DOFsManager.DOFsData,
                            reference_element_data: polydim.pde_tools.local_space_mcc_2_d.ReferenceElement_Data,
                            local_space_data: polydim.pde_tools.local_space_mcc_2_d.LocalSpace_Data,
                            test: ITest,
                            assembler_data: EllipticMCC2DProblemData) -> None:

        # Assemble strong boundary condition on Cell1Ds
        for ed in range(mesh.cell2_d_number_edges(cell2_d_index)):

            cell1_d_index = mesh.cell2_d_edge(cell2_d_index, ed)

            boundary_info = mesh_do_fs_info.cells_boundary_info[1][cell1_d_index]
            local_do_fs = do_fs_data.cells_do_fs[1][cell1_d_index]

            if (boundary_info.type != polydim.pde_tools.do_fs.DOFsManager.MeshDOFsInfo.BoundaryInfo.BoundaryTypes.strong
                    or len(local_do_fs) == 0):
                continue

            edge_do_fs_coordinates = polydim.pde_tools.local_space_mcc_2_d.edge_dofs_coordinates(reference_element_data, local_space_data, ed)

            neumann_values = test.strong_boundary_condition(boundary_info.Marker, edge_do_fs_coordinates.points)
            strong_boundary_values = polydim.pde_tools.local_space_mcc_2_d.edge_dofs(reference_element_data, local_space_data, ed, edge_do_fs_coordinates, neumann_values)


            assert len(local_do_fs) == len(strong_boundary_values)

            for loc_i in range(len(local_do_fs)):

                local_dof_i = local_do_fs[loc_i]

                match local_dof_i.type:
                    case polydim.pde_tools.do_fs.DOFsManager.DOFsData.DOF.Types.strong:
                        assembler_data.solution_neumann[local_dof_i.global_index] = strong_boundary_values[loc_i]
                        pass
                    case polydim.pde_tools.do_fs.DOFsManager.DOFsData.DOF.Types.dof:
                        pass
                    case _:
                        raise ValueError("Unknown DOF Type")

    def compute_weak_term(self,
                          cell2_d_index: int,
                          mesh: gedim.MeshMatricesDAO,
                          mesh_do_fs_info: polydim.pde_tools.do_fs.DOFsManager.MeshDOFsInfo,
                          do_fs_data: polydim.pde_tools.do_fs.DOFsManager.DOFsData,
                          reference_element_data: polydim.pde_tools.local_space_mcc_2_d.ReferenceElement_Data,
                          local_space_data: polydim.pde_tools.local_space_mcc_2_d.LocalSpace_Data,
                          test: ITest,
                          assembler_data: EllipticMCC2DProblemData) -> None:

        num_vertices = mesh.cell2_d_number_vertices(cell2_d_index)

        for ed in range(num_vertices):
            cell1_d_index = mesh.cell2_d_edge(cell2_d_index, ed)

            boundary_info = mesh_do_fs_info.cells_boundary_info[1][cell1_d_index]

            if boundary_info.type != polydim.pde_tools.do_fs.DOFsManager.MeshDOFsInfo.BoundaryInfo.BoundaryTypes.weak:
                continue

            # compute values of Neumann condition
            quadrature = polydim.pde_tools.local_space_mcc_2_d.edge_quadrature(reference_element_data, local_space_data, ed)
            dirichlet_values = test.weak_boundary_condition(boundary_info.marker, quadrature.points)
            basis_function_values = polydim.pde_tools.local_space_mcc_2_d.velocity_basis_functions_values_on_edges(ed, reference_element_data, local_space_data, quadrature.points)

            weak_boundary_values = -basis_function_values.T @ np.diag(quadrature.weights) @ dirichlet_values

            local_do_fs = do_fs_data.cells_do_fs[1][cell1_d_index]
            for loc_i in range(len(local_do_fs)):
                local_dof_i = local_do_fs[loc_i]

                match local_dof_i.type:
                    case polydim.pde_tools.do_fs.DOFsManager.DOFsData.DOF.Types.strong:
                        continue
                    case polydim.pde_tools.do_fs.DOFsManager.DOFsData.DOF.Types.dof:
                        assembler_data.right_hand_side[local_dof_i.global_index] += weak_boundary_values[loc_i]
                        pass
                    case _:
                        raise ValueError("Unknown DOF Type")

    @staticmethod
    def solve(count_do_fs_data: assembler_utilities.CountDOFsData,
              assembler_data: EllipticMCC2DProblemData) -> None:

        if count_do_fs_data.num_total_strongs > 0:
            assembler_data.right_hand_side -= assembler_data.neumann_matrix_a @ assembler_data.solution_neumann

        assembler_data.solution = sla.spsolve(assembler_data.global_matrix_a.tocsc(), assembler_data.right_hand_side)


    def assemble(self,
                 geometry_utilities_confi: gedim.GeometryUtilitiesConfig,
                 mesh: gedim.MeshMatricesDAO,
                 mesh_geometric_data: gedim.MeshUtilities.MeshGeometricData2D,
                 mesh_do_fs_info: List[polydim.pde_tools.do_fs.DOFsManager.MeshDOFsInfo],
                 do_fs_data: List[polydim.pde_tools.do_fs.DOFsManager.DOFsData],
                 do_fs_data_indices: List[polydim.pde_tools.do_fs.DOFsManager.CellsDOFsIndicesData],
                 count_do_fs_data: assembler_utilities.CountDOFsData,
                 reference_element_data: polydim.pde_tools.local_space_mcc_2_d.ReferenceElement_Data,
                 test: ITest) -> EllipticMCC2DProblemData:

        result = self.EllipticMCC2DProblemData()

        result.right_hand_side = np.zeros([count_do_fs_data.num_total_do_fs])
        result.solution_neumann = np.zeros([count_do_fs_data.num_total_strongs])

        assembler_utils_obj = assembler_utilities()

        for c in range(mesh.cell2_d_total_number()):

            local_space_data = polydim.pde_tools.local_space_mcc_2_d.create_local_space(geometry_utilities_confi.tolerance1_d,
                                                                                        geometry_utilities_confi.tolerance2_d,
                                                                                        mesh_geometric_data,
                                                                                        c,
                                                                                        reference_element_data)

            velocity_basis_functions_values = polydim.pde_tools.local_space_mcc_2_d.velocity_basis_functions_values(reference_element_data, local_space_data)
            velocity_basis_functions_divergence_values = polydim.pde_tools.local_space_mcc_2_d.velocity_basis_functions_divergence_values(reference_element_data, local_space_data)
            pressure_basis_functions_values = polydim.pde_tools.local_space_mcc_2_d.pressure_basis_functions_values(reference_element_data, local_space_data)

            cell2_d_internal_quadrature = polydim.pde_tools.local_space_mcc_2_d.internal_quadrature(reference_element_data, local_space_data)

            inverse_diffusion_term_values = test.inverse_diffusion_term(cell2_d_internal_quadrature.points)
            mixed_advection_term_values = test.mixed_advection_term(cell2_d_internal_quadrature.points)
            reaction_term_values = test.reaction_term(cell2_d_internal_quadrature.points)
            source_term_values = test.source_term(cell2_d_internal_quadrature.points)
            weights = cell2_d_internal_quadrature.weights

            local_a = np.zeros([velocity_basis_functions_values[0].shape[1], velocity_basis_functions_values[0].shape[1]])

            # coefficient for vem stabilizing form
            k_max = 0.0

            for d1 in range(2):
                for d2 in range(2):
                    local_a += (velocity_basis_functions_values[d1].T
                                @ np.diag(weights * inverse_diffusion_term_values[d1 + 3 * d2])
                                @ velocity_basis_functions_values[d2])

                    k_values_max = np.max(inverse_diffusion_term_values[d1 + 3 * d2])
                    if k_values_max > k_max:
                        k_max = k_values_max

            # stiffness matrix
            local_a += k_max * polydim.pde_tools.local_space_mcc_2_d.stabilization_matrix(reference_element_data, local_space_data)

            # reaction
            local_m =  pressure_basis_functions_values.T @ np.diag(weights * reaction_term_values) @ pressure_basis_functions_values

            # advection
            local_t = pressure_basis_functions_values.T @ np.diag(weights * mixed_advection_term_values[0]) @ velocity_basis_functions_values[0]

            for d in np.arange(1, 2):
                local_t += pressure_basis_functions_values.T @ np.diag(weights * mixed_advection_term_values[d]) @ velocity_basis_functions_values[d]

            # divergence
            local_b = pressure_basis_functions_values.T @ np.diag(weights) @ velocity_basis_functions_divergence_values

            local_rhs = np.concatenate([np.zeros([velocity_basis_functions_values[0].shape[1]]),
                                        pressure_basis_functions_values.T @ np.diag(weights) @ source_term_values], axis=0)


            local_lhs = np.block([[local_a, -(local_t + local_b).T], [local_b, local_m]])

            global_do_fs = do_fs_data[0].cells_global_do_fs[2][c]
            assert polydim.pde_tools.local_space_mcc_2_d.velocity_size(reference_element_data, local_space_data) == len(global_do_fs)


            local_count_do_fs = assembler_utils_obj.local_count_do_fs(2, c, do_fs_data)

            local_to_global_data = assembler_utilities.LocalMatrixToGlobalMatrixDOFsData()
            local_to_global_data.do_fs_data_indices = do_fs_data_indices
            local_to_global_data.local_offsets = local_count_do_fs.offset_do_fs
            local_to_global_data.global_offsets_do_fs = count_do_fs_data.offset_do_fs
            local_to_global_data.global_offsets_strongs = count_do_fs_data.offset_strongs

            assembler_utilities.assemble_local_matrix_to_global_matrix(2,
                                                                       c,
                                                                       local_to_global_data,
                                                                       local_to_global_data,
                                                                       local_lhs,
                                                                       result.global_matrix_a_data,
                                                                       result.neumann_matrix_a_data,
                                                                       local_rhs,
                                                                       result.right_hand_side)


            self.compute_strong_term(c, mesh, mesh_do_fs_info[0], do_fs_data[0], reference_element_data, local_space_data,
                                     test, result)
            self.compute_weak_term(c, mesh, mesh_do_fs_info[0], do_fs_data[0], reference_element_data,
                                   local_space_data, test, result)

        result.global_matrix_a = result.global_matrix_a_data.create(count_do_fs_data.num_total_do_fs, count_do_fs_data.num_total_do_fs)
        result.dirichlet_matrix_a = result.neumann_matrix_a_data.create(count_do_fs_data.num_total_do_fs, count_do_fs_data.num_total_do_fs)

        return result


    def post_process_solution(self,
                              geometry_utilities_confi: gedim.GeometryUtilitiesConfig,
                              mesh: gedim.MeshMatricesDAO,
                              mesh_geometric_data: gedim.MeshUtilities.MeshGeometricData2D,
                              do_fs_data: List[polydim.pde_tools.do_fs.DOFsManager.DOFsData],
                              do_fs_data_indices: List[polydim.pde_tools.do_fs.DOFsManager.CellsDOFsIndicesData],
                              count_do_fs_data: assembler_utilities.CountDOFsData,
                              reference_element_data: polydim.pde_tools.local_space_mcc_2_d.ReferenceElement_Data,
                              assembler_data: EllipticMCC2DProblemData,
                              test: ITest) -> PostProcessData:

        result = self.PostProcessData()

        result.residual_norm = 0.0
        if count_do_fs_data.num_total_do_fs > 0:
            residual = assembler_data.global_matrix_a @ assembler_data.solution - assembler_data.right_hand_side
            result.residual_norm = np.linalg.norm(residual)


        result.cell2_ds_numeric_pressure = np.zeros([mesh.cell2_d_total_number()])
        result.cell2_ds_exact_pressure = np.zeros([mesh.cell2_d_total_number()])
        result.cell2_ds_error_l2_pressure = np.zeros([mesh.cell2_d_total_number()])
        result.cell2_ds_super_error_l2_pressure = np.zeros([mesh.cell2_d_total_number()])
        result.cell2_ds_norm_l2_pressure = np.zeros([mesh.cell2_d_total_number()])
        result.cell2_ds_error_l2_velocity = np.zeros([mesh.cell2_d_total_number()])
        result.cell2_ds_norm_l2_velocity = np.zeros([mesh.cell2_d_total_number()])

        assembler_utilities_obj = assembler_utilities()
        for c in range(mesh.cell2_d_total_number()):

            local_space_data=polydim.pde_tools.local_space_mcc_2_d.create_local_space(
                geometry_utilities_confi.tolerance1_d,
                geometry_utilities_confi.tolerance2_d,
                mesh_geometric_data,
                c,
                reference_element_data)

            pressure_basis_functions_values = polydim.pde_tools.local_space_mcc_2_d.pressure_basis_functions_values(
                reference_element_data, local_space_data)
            velocity_basis_functions_derivative_values = polydim.pde_tools.local_space_mcc_2_d.velocity_basis_functions_values(
                reference_element_data, local_space_data)

            cell2_d_internal_quadrature = polydim.pde_tools.local_space_mcc_2_d.internal_quadrature(reference_element_data, local_space_data)


            local_count_do_fs = assembler_utilities_obj.local_count_do_fs(2, c, do_fs_data)
            do_fs_values = assembler_utilities_obj.global_solution_to_local_solution(2,
                                                                                     c,
                                                                                     do_fs_data_indices,
                                                                                     count_do_fs_data,
                                                                                     local_count_do_fs,
                                                                                     assembler_data.solution,
                                                                                     assembler_data.solution_neumann)

            do_fs_velocity_values = do_fs_values[0:local_count_do_fs.offset_do_fs[1]]
            do_fs_pressure_values = do_fs_values[local_count_do_fs.offset_do_fs[1]:]

            exact_pressure_values = test.exact_pressure(cell2_d_internal_quadrature.points)
            exact_velocity_values = test.exact_velocity(cell2_d_internal_quadrature.points)

            # Export solution
            result.cell2_ds_numeric_pressure[c] = np.sum(cell2_d_internal_quadrature.weights * (pressure_basis_functions_values @ do_fs_pressure_values))
            result.cell2_ds_exact_pressure[c] = np.sum(cell2_d_internal_quadrature.weights * exact_pressure_values)

            # Error l2 pressure
            local_error_l2_pressure = (pressure_basis_functions_values @ do_fs_pressure_values - exact_pressure_values)**2
            local_norm_l2_pressure = (pressure_basis_functions_values @ do_fs_pressure_values)**2

            result.cell2_ds_error_l2_pressure[c] = np.sum(cell2_d_internal_quadrature.weights * local_error_l2_pressure)
            result.cell2_ds_norm_l2_pressure[c] = np.sum(cell2_d_internal_quadrature.weights * local_norm_l2_pressure)

            # Interpolate Exact Pressure - Super Convergence Error
            rhs = pressure_basis_functions_values.T @ np.diag(cell2_d_internal_quadrature.weights) @ exact_pressure_values
            lhs = pressure_basis_functions_values.T @ np.diag(cell2_d_internal_quadrature.weights) @ pressure_basis_functions_values
            interpolant_coefficients = np.linalg.solve(lhs, rhs)

            local_super_error_l2_pressure = (pressure_basis_functions_values @ do_fs_pressure_values - pressure_basis_functions_values @ interpolant_coefficients)**2
            result.cell2_ds_super_error_l2_pressure[c] = np.sum(cell2_d_internal_quadrature.weights * local_super_error_l2_pressure)

            local_error_l2_velocity = ((velocity_basis_functions_derivative_values[0] @ do_fs_velocity_values - exact_velocity_values[0])**2
                              + (velocity_basis_functions_derivative_values[1] @ do_fs_velocity_values - exact_velocity_values[1])**2)

            local_norm_l2_velocity = ((velocity_basis_functions_derivative_values[0] @ do_fs_velocity_values)**2 +
                             (velocity_basis_functions_derivative_values[1] @ do_fs_velocity_values)**2)

            result.cell2_ds_error_l2_velocity[c] = np.sum(cell2_d_internal_quadrature.weights * local_error_l2_velocity)
            result.cell2_ds_norm_l2_velocity[c] = np.sum(cell2_d_internal_quadrature.weights * local_norm_l2_velocity)

            if mesh_geometric_data.cell2_ds_diameters[c] > result.mesh_size:
                result.mesh_size = mesh_geometric_data.cell2_ds_diameters[c]


        result.error_l2_pressure = np.sqrt(np.sum(result.cell2_ds_error_l2_pressure ))
        result.super_error_l2_pressure = np.sqrt(np.sum(result.cell2_ds_super_error_l2_pressure))
        result.norm_l2_pressure  = np.sqrt(np.sum(result.cell2_ds_norm_l2_pressure ))
        result.error_l2_velocity = np.sqrt(np.sum(result.cell2_ds_error_l2_velocity))
        result.norm_l2_velocity = np.sqrt(np.sum(result.cell2_ds_norm_l2_velocity))

        return result

