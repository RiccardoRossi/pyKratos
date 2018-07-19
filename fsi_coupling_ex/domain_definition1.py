from __future__ import print_function, absolute_import, division 
import sys
sys.path.append("..")
print(sys.path)

from numpy import *
from pyKratos import *
import domain_solver
        

#######################################################################################
#######################################################################################
# DOMAIN 1
#######################################################################################
#######################################################################################
# 9 - 8
# | / |
# 7 - 6
# | / |
# 5 - 4
# | / |
# 3 - 2
# |  /|
# | / |
# 0 - 1

node_list = {
    0: array([0.0, 0.0]),
    1: array([0.1, 0.0]),
    2: array([0.1, 0.2]),
    3: array([0.0, 0.2]),
    4: array([0.1, 0.3]),
    5: array([0.0, 0.3]),
    6: array([0.1, 0.4]),
    7: array([0.0, 0.4]),
    8: array([0.1, 0.5]),
    9: array([0.0, 0.5]),
}

property_list = {
    0: {YOUNG_MODULUS: 200e4,
        POISSON_RATIO: 0.3,
        DENSITY: 100.0,
        BODY_FORCE_X: 0.0,
        BODY_FORCE_Y: 0.0,
    }
}


element_connectivities = {
    1: [0, [0, 1, 2]],
    2: [0, [0, 2, 3]],
    3: [0, [2, 4, 3]],
    4: [0, [3, 4, 5]],
    5: [0, [5, 4, 6]],
    6: [0, [5, 6, 7]],
    7: [0, [7, 6, 8]],
    8: [0, [7, 8, 9]],
}

condition_connectivities = {
    1: [0, [0]],
    2: [0, [1]],
    3: [0, [2]],
    4: [0, [3]],
    5: [0, [4]],
    6: [0, [5]],
    7: [0, [6]],
    8: [0, [7]],
    9: [0, [8]],
    10: [0, [9]]
}


t = []
d = []

#import object Domain1
domain1_solver = domain_solver.DomainSolver(node_list, property_list, element_connectivities, condition_connectivities)

domain1_solver.ApplyDirichletCondition(
        {
            VELOCITY_X: [
                    [0, True, 0.0],  # first column is Id of node, second col if fixity, third is imposed value
                    [1, True, 0.0],
                ],
            VELOCITY_Y: [
                [0, True, 0.0],
                [1, True, 0.0],
            ],
            FSI_INTERFACE: [
                [0, False, 1.0],
                [3, False, 1.0],
                [5, False, 1.0],
                [7, False, 1.0],
                [9, False, 1.0],
            ]
        }
    )

domain1_solver.ApplyNodalForces(
        {
            8 : [1.0, 0.0]
        }
    )

