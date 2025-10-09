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
                                       [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]])

        # create edges
        pde_domain.edges = np.zeros([2, 12])
        pde_domain.edges[:, 0] = [0, 1]
        pde_domain.edges[:, 1] = [1, 2]
        pde_domain.edges[:, 2] = [2, 3]
        pde_domain.edges[:, 3] = [3, 0]
        pde_domain.edges[:, 4] = [4, 5]
        pde_domain.edges[:, 5] = [5, 6]
        pde_domain.edges[:, 6] = [6, 7]
        pde_domain.edges[:, 7] = [7, 4]
        pde_domain.edges[:, 8] = [0, 4]
        pde_domain.edges[:, 9] = [1, 5]
        pde_domain.edges[:, 10] = [2, 6]
        pde_domain.edges[:, 11] = [3, 7]

        # create faces
        pde_domain.faces = []
        face = np.zeros([2, 4])
        pde_domain.faces[0].row(0) << 0, 1, 2, 3;
        pde_domain.faces[0].row(1) << 0, 1, 2, 3;

        pde_domain.faces[1].row(0) << 4, 5, 6, 7;
        pde_domain.faces[1].row(1) << 4, 5, 6, 7;

        pde_domain.faces[2].row(0) << 0, 3, 7, 4;
        pde_domain.faces[2].row(1) << 3, 11, 7, 8;

        pde_domain.faces[3].row(0) << 1, 2, 6, 5;
        pde_domain.faces[3].row(1) << 1, 10, 5, 9;

        pde_domain.faces[4].row(0) << 0, 1, 5, 4;
        pde_domain.faces[4].row(1) << 0, 9, 4, 8;

        pde_domain.faces[5].row(0) << 3, 2, 6, 7;
        pde_domain.faces[5].row(1) << 2, 10, 6, 11;

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
        raise ValueError("not valid marker")

    def exact_solution(self, points: np.ndarray):
        return 16.0 * points[0, :] * (1.0 - points[0, :]) * points[1, :] * (1.0 - points[1, :]) + 1.1

    def exact_derivative_solution(self, points: np.ndarray):
        return [16.0 * (1.0 - 2.0 * points[0, :]) * points[1, :] * (1.0 - points[1, :]),
                16.0 * points[0, :] * (1.0 - points[0, :]) * (1.0 - 2.0 * points[1, :])]