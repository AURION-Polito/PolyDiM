from Elliptic_PCC_2D.test_definition import ITest, EllipticPolynomialProblem
from pypolydim import gedim, polydim
def create_test(test_id: int) -> ITest:

    match test_id:
        case 1:
            return EllipticPolynomialProblem()
        case _:
            raise ValueError("not valid test id")

def create_mesh(geometry_utilities: gedim.GeometryUtilities, mesh_utilities: gedim.MeshUtilities,
                mesh_type: polydim.pde_tools.mesh.pde_mesh_utilities.MeshGenerator_Types_2D,
                mesh_max_relative_area: float, import_path: str,
                pde_domain: polydim.pde_tools.mesh.pde_mesh_utilities.PDE_Domain_2D,
                mesh: gedim.MeshMatricesDAO):

    if (polydim.pde_tools.mesh.pde_mesh_utilities.MeshGenerator_Types_2D.triangular == mesh_type or
            polydim.pde_tools.mesh.pde_mesh_utilities.MeshGenerator_Types_2D.minimal == mesh_type or
            polydim.pde_tools.mesh.pde_mesh_utilities.MeshGenerator_Types_2D.polygonal == mesh_type or
            polydim.pde_tools.mesh.pde_mesh_utilities.MeshGenerator_Types_2D.squared == mesh_type or
        polydim.pde_tools.mesh.pde_mesh_utilities.MeshGenerator_Types_2D.random_distorted == mesh_type):

        polydim.pde_tools.mesh.pde_mesh_utilities.create_mesh_2_d(geometry_utilities,
                                                                  mesh_utilities,
                                                                  mesh_type,
                                                                  pde_domain,
                                                                  mesh_max_relative_area,
                                                                  mesh)

    elif (polydim.pde_tools.mesh.pde_mesh_utilities.MeshGenerator_Types_2D.csv_importer == mesh_type or
            polydim.pde_tools.mesh.pde_mesh_utilities.MeshGenerator_Types_2D.off_importer == mesh_type):

        polydim.pde_tools.mesh.pde_mesh_utilities.import_mesh_2_d(mesh_utilities,
                                                                  mesh_type,
                                                                  import_path,
                                                                  mesh)
    else:
        raise ValueError("MeshGenerator " + str(mesh_type) + " not supported")
