import argparse
import numpy
from pypolydim import polydim, gedim
from Elliptic_PCC_2D.program_utilities import create_test, create_mesh

if __name__=='__main__':

    parser =argparse.ArgumentParser()
    parser.add_argument('-order','--method-order',dest='method_order', default=1, type=int, help="Method order")
    parser.add_argument('-method','--method-type',dest='method_type', default=1, type=int, help="Method type")
    parser.add_argument('-test', '--test-id', dest='test_id', default=1, type=int, help="Test type")
    parser.add_argument('-mesh', '--mesh-type', dest='mesh_type', default=5, type=int, help="Mesh type")
    parser.add_argument('-area', '--mesh-max-relative-area', dest='max_relative_area', default=0.1, type=float, help="Mesh max relative area")
    parser.add_argument('-import', '--import-path', dest='import_path', default='./', type=str, help="Mesh Import Path")
    args = parser.parse_args()

    mesh_type = polydim.pde_tools.mesh.pde_mesh_utilities.MeshGenerator_Types_2D(args.mesh_type)

    # Set problem
    test = create_test(args.test_id)
    pde_domain = test.domain()

    # Create mesh
    geometry_utilities_config = gedim.GeometryUtilitiesConfig()
    geometry_utilities = gedim.GeometryUtilities(geometry_utilities_config)
    mesh_utilities = gedim.MeshUtilities()

    mesh_data = gedim.MeshMatrices()
    mesh = gedim.MeshMatricesDAO(mesh_data)

    create_mesh(geometry_utilities, mesh_utilities, mesh_type, args.max_relative_area, args.import_path, pde_domain, mesh)




