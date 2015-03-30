from __future__ import print_function, absolute_import, division 
import sys
sys.path.append("..")
print(sys.path)

from numpy import *
from pyKratos import *

#example = "cavity" 
#example = "gravity"
#example = "shear_x"
example = "inlet"

# add variables to be allocated from the list in variables.py
solution_step_variables = [
    VELOCITY_X,
    VELOCITY_Y,
    PRESSURE,
    IS_LAGRANGIAN,
    EXTERNAL_FORCE_X,
    EXTERNAL_FORCE_Y
]

if example == "gravity":
    property_list = {
        0: {VISCOSITY: 1.0,
            DENSITY: 1.0,
            BODY_FORCE_X: 0.0,
            BODY_FORCE_Y: -10.0,
            }
        }
else:
    property_list = {
        0: {VISCOSITY: 1.0,
            DENSITY: 1.0,
            BODY_FORCE_X: 0.0,
            BODY_FORCE_Y: 0.0,
            }
        }
        
#defining a 1 by 1 square
nx = 11
dx = 1.0/(nx-1)
ny = 11
dy = 1.0/(ny-1)
import generate_square
node_list, element_connectivities,face_connectivities = generate_square.GenerateSquare(nx, dx, ny, dy)



#import py_kratos
buffer_size = 3  # store current step and 2 in the past
model_part = ModelPart(buffer_size, solution_step_variables)
model_part.AddNodes(node_list)
model_part.AddProperties(property_list)
model_part.AddElements("stokes_element_2d", element_connectivities)
model_part.AddConditions("neumann_face_condition_2d", face_connectivities)
gid_io = GidIO("cavity.mdpa", "cavity")

bottom_nodes = []
left_nodes = []
right_nodes = []
top_nodes = []

for node in model_part.NodeIterators():
    if(node.coordinates[0] <= 0.000001):
        left_nodes.append(node)
    elif  node.coordinates[0] >(nx-1)*dx-0.00001:
        right_nodes.append(node)
    elif(node.coordinates[1] < 0.00001):
        bottom_nodes.append(node)
    elif(node.coordinates[1] >(ny-1)*dy-0.0001):
        top_nodes.append(node)
        

if example == "cavity":
    for node in bottom_nodes:
        node.Fix(VELOCITY_X)
        node.Fix(VELOCITY_Y)
        
    for node in left_nodes:
        node.Fix(VELOCITY_X)
        node.Fix(VELOCITY_Y)
    for node in right_nodes:
        node.Fix(VELOCITY_X)
        node.Fix(VELOCITY_Y)
    for node in top_nodes:
        node.Fix(VELOCITY_X)
        node.Fix(VELOCITY_Y)
        node.SetSolutionStepValue(VELOCITY_X,0,1.0)
    #fixing the node at the center of the bottom face
    model_part.Nodes[int(nx/2)+1].Fix(PRESSURE)
elif example=="gravity":
    for node in bottom_nodes:
        node.Fix(VELOCITY_X)
        node.Fix(VELOCITY_Y)
              
    for node in left_nodes:
        node.Fix(VELOCITY_X)
        node.Fix(VELOCITY_Y)
    for node in right_nodes:
        node.Fix(VELOCITY_X)
        node.Fix(VELOCITY_Y)
        
    for node in top_nodes:
        pass
elif example == "shear_x":
    for node in model_part.NodeIterators():
        if(node.coordinates[1] == 0.0):
            node.Fix(VELOCITY_X)
            node.Fix(VELOCITY_Y)
        elif(node.coordinates[1] == 1.0):
            node.Fix(VELOCITY_X)
            node.Fix(VELOCITY_Y)
            node.SetSolutionStepValue(VELOCITY_X,0,1.0)
elif example == "inlet":
    for node in model_part.NodeIterators():
        if(node.coordinates[0] == 0.0):
            node.Fix(VELOCITY_X)
            node.Fix(VELOCITY_Y)
            node.SetSolutionStepValue(VELOCITY_X,0,1.0)
        
        if(node.coordinates[1] == (ny-1)*dy or node.coordinates[1] == 0.0):
            node.Fix(VELOCITY_X)
            node.Fix(VELOCITY_Y)
            node.SetSolutionStepValue(VELOCITY_X,0,0.0) ##this overwrites the one before on corner nodes
            
        for node in right_nodes:
            node.SetSolutionStepValue(EXTERNAL_FORCE_X,0,-100.0)
            node.SetSolutionStepValue(EXTERNAL_FORCE_Y,0,0.0)
            
         


import gear_scheme
time_scheme = gear_scheme.GearScheme(model_part)

builder_and_solver = builder_and_solver.BuilderAndSolver(
    model_part, time_scheme)

strategy = solving_strategy.SolvingStrategy(
    model_part, time_scheme, builder_and_solver)

strategy.Initialize()

mesh_name = "test"
gid_io.WriteMesh(model_part, mesh_name)
dt = 0.1
nsteps = 10
for i in range(1,nsteps):
    time = i*dt
    model_part.CloneTimeStep(time)
    print("time = ", time)
    strategy.Solve()
    gid_io.WriteNodalResults(PRESSURE, model_part.NodeIterators(), time)
    gid_io.WriteNodalResults(VELOCITY, model_part.NodeIterators(), time)
    #plot_contour.PlotContour(model_part.NodeIterators(), VELOCITY_Y, outname )
    
    #for node in model_part.NodeIterators():
        #if(node.coordinates[1] < 0.00001):
            #print(node.GetSolutionStepValue(TEMPERATURE,0))

# strategy.SpyMatrix()
#import plot_contour
#plot_contour.PlotContour(model_part.NodeIterators(), TEMPERATURE)