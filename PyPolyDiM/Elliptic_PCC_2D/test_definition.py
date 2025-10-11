from abc import ABC, abstractmethod
import numpy as np
from pypolydim import polydim

class ITest(ABC):

    @staticmethod
    @abstractmethod
    def domain() -> polydim.pde_tools.mesh.pde_mesh_utilities.PDE_Domain_2D:
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
    def domain() -> polydim.pde_tools.mesh.pde_mesh_utilities.PDE_Domain_2D:
        pde_domain = polydim.pde_tools.mesh.pde_mesh_utilities.PDE_Domain_2D()

        pde_domain.vertices = np.array([[0.0, 1.0, 1.0, 0.0],
                                       [0.0, 0.0, 1.0, 1.0],
                                       [0.0, 0.0, 0.0, 0.0]])

        pde_domain.area = 1.0
        pde_domain.shape_type = polydim.pde_tools.mesh.pde_mesh_utilities.PDE_Domain_2D.Domain_Shape_Types.parallelogram

        return pde_domain

    @staticmethod
    def boundary_info():

        info_internal = polydim.pde_tools.do_fs.DOFsManager.MeshDOFsInfo.BoundaryInfo(polydim.pde_tools.do_fs.DOFsManager.MeshDOFsInfo.BoundaryInfo.BoundaryTypes.none)
        info_internal.marker = 0

        info_dirichlet = polydim.pde_tools.do_fs.DOFsManager.MeshDOFsInfo.BoundaryInfo(
            polydim.pde_tools.do_fs.DOFsManager.MeshDOFsInfo.BoundaryInfo.BoundaryTypes.strong)
        info_dirichlet.marker = 1

        info_neumann_top = polydim.pde_tools.do_fs.DOFsManager.MeshDOFsInfo.BoundaryInfo(
            polydim.pde_tools.do_fs.DOFsManager.MeshDOFsInfo.BoundaryInfo.BoundaryTypes.weak)
        info_neumann_top.marker = 4

        info_neumann_right = polydim.pde_tools.do_fs.DOFsManager.MeshDOFsInfo.BoundaryInfo(
            polydim.pde_tools.do_fs.DOFsManager.MeshDOFsInfo.BoundaryInfo.BoundaryTypes.weak)
        info_neumann_right.marker = 2

        info_neumann_none = polydim.pde_tools.do_fs.DOFsManager.MeshDOFsInfo.BoundaryInfo(
            polydim.pde_tools.do_fs.DOFsManager.MeshDOFsInfo.BoundaryInfo.BoundaryTypes.none)
        info_neumann_none.marker = 0

        return {
            0: info_internal,
            1: info_dirichlet,
            2: info_dirichlet,
            3: info_neumann_none,
            4: info_dirichlet,
            5: info_dirichlet,
            6: info_neumann_right,
            7: info_neumann_top,
            8: info_dirichlet
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

        derivatives = self.exact_derivative_solution(points)

        match marker:
            case 2:
                return derivatives[0]
            case 4:
                return derivatives[1]
            case _:
                raise ValueError("unknown marker")

    def exact_solution(self, points: np.ndarray):
        return 16.0 * points[0, :] * (1.0 - points[0, :]) * points[1, :] * (1.0 - points[1, :]) + 1.1

    def exact_derivative_solution(self, points: np.ndarray):
        return [16.0 * (1.0 - 2.0 * points[0, :]) * points[1, :] * (1.0 - points[1, :]),
                16.0 * points[0, :] * (1.0 - points[0, :]) * (1.0 - 2.0 * points[1, :])]