from __future__ import print_function, absolute_import, division 
import sys
sys.path.append("..")
print(sys.path)

from numpy import *
from pyKratos import *

class DomainSolver:

    def __init__(self, node_list, property_list, element_connectivities, condition_connectivities):
        # add variables to be allocated from the list in variables.py
        solution_step_variables = [
            EXTERNAL_FORCE_X,
            EXTERNAL_FORCE_Y,
            DISPLACEMENT_X,
            DISPLACEMENT_Y,
            VELOCITY_X,
            VELOCITY_Y,
            IS_LAGRANGIAN,
            REACTION_X,
            REACTION_Y,
            FSI_INTERFACE
        ]



        buffer_size = 3  # store current step and 2 in the past
        self.model_part = ModelPart(buffer_size, solution_step_variables)
        self.model_part.AddNodes(node_list)
        self.model_part.AddProperties(property_list)
        self.model_part.AddElements("triangular_elastic_element", element_connectivities)
        self.model_part.AddConditions("point_condition_2d", condition_connectivities)

        for node in self.model_part.NodeIterators():
            node.SetSolutionStepValue(IS_LAGRANGIAN,0,True)
        
        import gear_scheme
        self.time_scheme = gear_scheme.GearScheme(self.model_part)

        #import pyKratos.builder_and_solver  

        self.builder = builder_and_solver.BuilderAndSolver(
            self.model_part, self.time_scheme)

        
        self.strategy = solving_strategy.SolvingStrategy(
            self.model_part, self.time_scheme, self.builder)

        self.strategy.Initialize() #we will redo this when we Apply Dirichlet


    def ComputeReactions(self):
        b = self.strategy.ComputeReactionsVector()

        counter = 0
        for node in self.model_part.NodeIterators():
            rx = b[node.EquationId(VELOCITY_X)]
            ry = b[node.EquationId(VELOCITY_Y)]
            node.SetSolutionStepValue(REACTION_X, 0, -rx)
            node.SetSolutionStepValue(REACTION_Y, 0, -ry)
            counter += 2

    def ApplyDirichletCondition(self, nodal_values):
        self.model_part.AddNodalValues(nodal_values)
        self.strategy.Initialize() #here we mount the matrix

    def ApplyNodalForces(self, nodal_forces):
        for key,value in nodal_forces.items():
            node = self.model_part.Nodes[key]

            print(key,value)
            node.SetSolutionStepValue(EXTERNAL_FORCE_X, 0, value[0])
            node.SetSolutionStepValue(EXTERNAL_FORCE_Y, 0, value[1])


    def GetReactions(self, node_list):
        reactions = {}
        for index in node_list:
            node = self.model_part.Nodes[index]
            reactions.update({index, [
                node.GetSolutionStepValue(REACTION_X), 
                node.GetSolutionStepValue(REACTION_Y)
                ]})
        return reactions
        

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
        DENSITY: 1000.0,
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
domain1 = DomainSolver(node_list, property_list, element_connectivities, condition_connectivities)

domain1.ApplyDirichletCondition(
        {
            VELOCITY_X: [
                    [0, True, 0.0],  # first column is Id of node, second col if fixity, third is imposed value
                    [1, True, 0.0],
                ],
            VELOCITY_Y: [
                [0, True, 0.0],
                [1, True, 0.0],
            ],
        }
    )

domain1.ApplyNodalForces(
        {
            8 : [1.0, 0.0]
        }
    )
