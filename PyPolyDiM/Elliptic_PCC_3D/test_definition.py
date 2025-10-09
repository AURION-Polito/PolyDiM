from abc import ABC, abstractmethod
import numpy as np
from pypolydim import polydim

class ITest(ABC):

    @staticmethod
    @abstractmethod
    def domain() -> polydim.pde_tools.mesh.pde_mesh_utilities.PDE_Domain_3D:
        pass

    @staticmethod
    @abstractmethod
    def boundary_info():
        pass

    @abstractmethod
    def diffusion_term(self, points: np.ndarray) -> np.ndarray:
        pass

    @abstractmethod
    def source_term(self, points: np.ndarray) -> np.ndarray:
        pass

    @abstractmethod
    def strong_boundary_condition(self, marker: int, points: np.ndarray):
        pass

    @abstractmethod
    def weak_boundary_condition(self, marker: int, points: np.ndarray):
        pass

    @abstractmethod
    def exact_solution(self, points: np.ndarray):
        pass

    @abstractmethod
    def exact_derivative_solution(self, points: np.ndarray):
        pass

class EllipticPolynomialProblem(ITest):

    @staticmethod
    def domain() -> polydim.pde_tools.mesh.pde_mesh_utilities.PDE_Domain_3D:
        pde_domain = polydim.pde_tools.mesh.pde_mesh_utilities.PDE_Domain_3D()

        pde_domain.vertices = np.array([[0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0],
                                       [0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0],
                                       [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]], dtype=np.float64)

        # create edges
        pde_domain.edges = np.array([[0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3],
                                       [1, 2, 3, 0, 5, 6, 7, 4, 4, 5, 6, 7]], dtype=np.int64)

        # create faces
        pde_domain.faces = [np.array([[0, 1, 2, 3], [0, 1, 2, 3]], dtype=np.int64),
                            np.array([[4, 5, 6, 7], [4, 5, 6, 7]], dtype=np.int64),
                            np.array([[0, 3, 7, 4], [3, 11, 7, 8]], dtype=np.int64),
                            np.array([[1, 2, 6, 5], [1, 10, 5, 9]], dtype=np.int64),
                            np.array([[0, 1, 5, 4], [0, 9, 4, 8]], dtype=np.int64),
                            np.array([[3, 2, 6, 7], [2, 10, 6, 11]], dtype=np.int64)]

        pde_domain.volume = 1.0
        pde_domain.shape_type = polydim.pde_tools.mesh.pde_mesh_utilities.PDE_Domain_3D.Domain_Shape_Types.parallelepiped

        return pde_domain

    @staticmethod
    def boundary_info():

        info_internal = polydim.pde_tools.do_fs.DOFsManager.MeshDOFsInfo.BoundaryInfo(polydim.pde_tools.do_fs.DOFsManager.MeshDOFsInfo.BoundaryInfo.BoundaryTypes.none)
        info_internal.marker = 0

        info_dirichlet = polydim.pde_tools.do_fs.DOFsManager.MeshDOFsInfo.BoundaryInfo(
            polydim.pde_tools.do_fs.DOFsManager.MeshDOFsInfo.BoundaryInfo.BoundaryTypes.strong)
        info_dirichlet.marker = 1

        return {
            0: info_internal,
            1: info_dirichlet,
            2: info_dirichlet,
            3: info_dirichlet,
            4: info_dirichlet,
            5: info_dirichlet,
            6: info_dirichlet,
            7: info_dirichlet,
            8: info_dirichlet,
            9: info_dirichlet,
            10: info_dirichlet,
            11: info_dirichlet,
            12: info_dirichlet,
            13: info_dirichlet,
            14: info_dirichlet,
            15: info_dirichlet,
            16: info_dirichlet,
            17: info_dirichlet,
            18: info_dirichlet,
            19: info_dirichlet,
            20: info_dirichlet,
            21: info_dirichlet,
            22: info_dirichlet,
            23: info_dirichlet,
            24: info_dirichlet,
            25: info_dirichlet,
            26: info_dirichlet
        }

    def diffusion_term(self, points: np.ndarray) -> np.ndarray:
        return np.ones(points.shape[1])

    def source_term(self, points: np.ndarray):
        return 32.0 * (points[0, :] * (1.0 - points[0, :]) + points[1, :] * (1.0 - points[1, :]))

    def strong_boundary_condition(self, marker: int, points: np.ndarray):

        if marker != 1:
            raise ValueError("not valid marker")

        return self.exact_solution(points)

    def weak_boundary_condition(self, marker: int, points: np.ndarray):
        raise ValueError("not valid marker")

    def exact_solution(self, points: np.ndarray):
        return 16.0 * points[0, :] * (1.0 - points[0, :]) * points[1, :] * (1.0 - points[1, :]) + 1.1

    def exact_derivative_solution(self, points: np.ndarray):
        return [16.0 * (1.0 - 2.0 * points[0, :]) * points[1, :] * (1.0 - points[1, :]),
                16.0 * points[0, :] * (1.0 - points[0, :]) * (1.0 - 2.0 * points[1, :]),
                np.zeros(points.shape[1])]