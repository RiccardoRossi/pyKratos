from __future__ import print_function, absolute_import, division 
import sys
sys.path.append("..")
print(sys.path)

from numpy import *
from variables import *

# add variables to be allocated from the list in variables.py
solution_step_variables = [
    DISPLACEMENT_X,
    DISPLACEMENT_Y,
    VELOCITY_X,
    VELOCITY_Y,
    IS_LAGRANGIAN,
]

#     2
#    /|
#   / |
# 0 - 1

node_list = {
    0: array([0.0, 0.0]),
    1: array([0.1, 0.0]),
    2: array([0.1, 0.2]),
}

property_list = {
    0: {YOUNG_MODULUS: 200e6,
        POISSON_RATIO: 0.3,
        DENSITY: 1000.0,
        BODY_FORCE_X: 10.0,
        BODY_FORCE_Y: 0.0,
       }
}


element_connectivities = {
    1: [0, [0, 1, 2]],
}

nodal_values = {
    VELOCITY_X: [
        [0, True, 0.0],  # first column is Id of node, second col if fixity, third is imposed value
        [1, True, 0.0],
    ],
VELOCITY_Y: [
    [0, True, 0.0],
    [1, True, 0.0],
],
}


from model_part import *
from gid_io import *

buffer_size = 3  # store current step and 2 in the past
model_part = ModelPart(buffer_size, solution_step_variables)
#model_part.AddNodes(node_list)
#model_part.AddProperties(property_list)
#model_part.AddElements("triangular_elastic_element", element_connectivities)
#model_part.AddNodalValues(nodal_values)
gid_io = GidIO("test.mdpa", "test")
gid_io.ReadModelPart(model_part)
print (model_part)

import solving_strategy
import gear_scheme
time_scheme = gear_scheme.GearScheme(model_part)

import builder_and_solver
builder_and_solver = builder_and_solver.BuilderAndSolver(
    model_part, time_scheme)

strategy = solving_strategy.SolvingStrategy(
    model_part, time_scheme, builder_and_solver)
strategy.Initialize()


t = []
d = []

import plot_contour
dt = 0.000025
nsteps = 10

mesh_name = "test"
gid_io.WriteMesh(model_part, mesh_name)
for i in range(1,nsteps):
    time = i*dt
    model_part.CloneTimeStep(time)
    print("time = ", time)
    strategy.Solve()

    gid_io.WriteNodalResults(DISPLACEMENT_X, model_part.NodeIterators(), time)
        
    #outname = "disp"+str(time)+".png"
    #plot_contour.PlotContour(model_part.NodeIterators(), DISPLACEMENT_X, outname )

    #outname = "vx"+str(time)+".png"
    #plot_contour.PlotContour(model_part.NodeIterators(), VELOCITY_X, outname )

    t.append(time)
    d.append(model_part.Nodes[2].GetSolutionStepValue(DISPLACEMENT_X,0))
    
print(d)

import matplotlib.pyplot as plt
plt.plot(t,d)
plt.show()