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

    # K
    @abstractmethod
    def diffusion_term(self, points: np.ndarray) -> np.ndarray:
        pass

    # K^{-1}
    @abstractmethod
    def inverse_diffusion_term(self, points: np.ndarray) -> np.ndarray:
        pass


    @abstractmethod
    def reaction_term(self, points: np.ndarray) -> np.ndarray:
        pass

    # beta
    @abstractmethod
    def advection_term(self, points: np.ndarray) -> np.ndarray:
        pass

    # k^{-1} beta
    @abstractmethod
    def mixed_advection_term(self, points: np.ndarray) -> np.ndarray:
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

    # p
    @abstractmethod
    def exact_pressure(self, points: np.ndarray):
        pass

    # u = - K \nabla p + beta u
    @abstractmethod
    def exact_velocity(self, points: np.ndarray):
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
            polydim.pde_tools.do_fs.DOFsManager.MeshDOFsInfo.BoundaryInfo.BoundaryTypes.weak)
        info_dirichlet.marker = 2

        return {
            0: info_internal,
            1: info_internal,
            2: info_internal,
            3: info_internal,
            4: info_internal,
            5: info_dirichlet,
            6: info_dirichlet,
            7: info_dirichlet,
            8: info_dirichlet
        }

    def diffusion_term(self, points: np.ndarray) -> np.ndarray:

        num_points = points.shape[1]
        diffusion_term_values = np.zeros([9, num_points])
        diffusion_term_values[0, :] = np.ones(num_points)
        diffusion_term_values[4, :] = np.ones(num_points)

        return diffusion_term_values

    def inverse_diffusion_term(self, points: np.ndarray) -> np.ndarray:

        num_points = points.shape[1]
        diffusion_term_values = np.zeros([9, num_points])
        diffusion_term_values[0, :] = np.ones(num_points)
        diffusion_term_values[4, :] = np.ones(num_points)

        return diffusion_term_values

    def reaction_term(self, points: np.ndarray) -> np.ndarray:

        num_points = points.shape[1]

        return np.zeros([num_points])

    def advection_term(self, points: np.ndarray) -> np.ndarray:

        num_points = points.shape[1]
        advection_term_values = np.zeros([3, num_points])

        return advection_term_values

    def mixed_advection_term(self, points: np.ndarray) -> np.ndarray:

        num_points = points.shape[1]
        advection_term_values = np.zeros([3, num_points])

        return advection_term_values

    def source_term(self, points: np.ndarray):
        return 2.0 * np.pi * np.pi * np.sin(np.pi * points[0, :]) * np.sin(np.pi * points[1, :])

    def strong_boundary_condition(self, marker: int, points: np.ndarray):

        match marker:
            case _:
                raise ValueError("unknown marker")

    def weak_boundary_condition(self, marker: int, points: np.ndarray):

        if marker != 2:
            raise ValueError("unknown marker")

        return self.exact_pressure(points)

    def exact_pressure(self, points: np.ndarray):
        return np.sin(np.pi * points[0, :]) * np.sin(np.pi * points[1, :])

    def exact_velocity(self, points: np.ndarray):
        return [-np.pi * np.cos(np.pi * points[0, :]) * np.sin(np.pi * points[1, :]),
                -np.pi * np.sin(np.pi * points[0, :]) * np.cos(np.pi * points[1, :])]