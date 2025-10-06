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
    def diffusion_term(self, points: np.array):
        pass

    @abstractmethod
    def source_term(self, points: np.array):
        pass

    @abstractmethod
    def strong_boundary_condition(self, marker: int, points: np.array):
        pass

    @abstractmethod
    def weak_boundary_condition(self, marker: int, points: np.array):
        pass

    @abstractmethod
    def exact_solution(self, points: np.array):
        pass

    @abstractmethod
    def exact_derivatives_solution(self, points: np.array):
        pass

class EllipticPolynomialProblem(ITest):

    @staticmethod
    def domain() -> polydim.pde_tools.mesh.pde_mesh_utilities.PDE_Domain_2D:
        pde_domain = polydim.pde_tools.mesh.pde_mesh_utilities.PDE_Domain_2D()

        pde_domain.vertices = np.array([[0.0, 1.0, 1.0, 0.0],
                                       [0.0, 0.0, 1.0, 1.0],
                                       [0.0, 0.0, 0.0, 0.0]])

        pde_domain.area = 1.0

        return pde_domain

    @staticmethod
    def boundary_info():
        return

    def diffusion_term(self, points: np.array):
        return np.ones(points.shape[1])

    def source_term(self, points: np.array):
        return 32.0 * points[0, :] * (1.0 - points[0, :]) * points[1, :] * (1.0 - points[1, :])

    def strong_boundary_condition(self, marker: int, points: np.array):

        if marker != 1:
            raise ValueError("not valid marker")

        return self.exact_solution(points)

    def weak_boundary_condition(self, marker: int, points: np.array):
        raise ValueError("not valid marker")

    def exact_solution(self, points: np.array):
        return 16.0 * points[0, :] * (1.0 - points[0, :]) * points[1, :] * (1.0 - points[1, :]) + 1.1

    def exact_derivatives_solution(self, points: np.array):
        return [16.0 * (1.0 - 2.0 * points[0, :]) * points[1, :] * (1.0 - points[1, :]),
                16.0 * points[0, :] * (1.0 - points[0, :]) * (1.0 - 2.0 * points[1, :]) + 1.1]