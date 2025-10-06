import argparse
import numpy
from pypolydim import polydim, gedim
from Elliptic_PCC_2D.program_utilities import create_test

if __name__=='__main__':

    parser =argparse.ArgumentParser()
    parser.add_argument('-order','--method-order',dest='method_order', default=1, type=int, help="Method order")
    parser.add_argument('-method','--method-type',dest='method_type', default=1, type=int, help="Method type")
    parser.add_argument('-test', '--test-id', dest='test_id', default=1, type=int, help="Test type")
    config = parser.parse_args()

    # Set problem
    test = create_test(config.test_id)

    # Create mesh
    geometry_utilities_config = gedim.GeometryUtilitiesConfig()
    geometry_utilities = gedim.GeometryUtilities(geometry_utilities_config)

