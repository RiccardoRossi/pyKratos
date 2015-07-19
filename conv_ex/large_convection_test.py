from __future__ import print_function, absolute_import, division 
import sys
sys.path.append("..")
from numpy import *
from pyKratos import *
# add variables to be allocated from the list in variables.py

# add variables to be allocated from the list in variables.py
solution_step_variables = [
    TEMPERATURE,
    VELOCITY_X,
    VELOCITY_Y,
    VOLUME_HEAT_SOURCE
]

property_list = {
    0: {CONDUCTIVITY: 0.001,
        }
}


nx = 40
dx = 1.0/(nx-1)
ny = 40
dy = 1.0/(ny-1)
import generate_square
node_list, element_connectivities, face_connectivities = generate_square.GenerateSquare(nx, dx, ny, dy)

#import py_kratos
#from model_part import *

buffer_size = 3  # store current step and 2 in the past
model_part = ModelPart(buffer_size, solution_step_variables)
model_part.AddNodes(node_list)
model_part.AddProperties(property_list)
model_part.AddElements(
    "conv_diff_element2D", element_connectivities)
# model_part.AddNodalValues(nodal_values)

# apply a velocity_y


v_y = 20.0
for node in model_part.NodeIterators():
    if(node.Id < nx):
        node.Fix(TEMPERATURE)
    if(node.Id < nx / 2):
        node.SetSolutionStepValue(TEMPERATURE, 0, 1.0)
    node.SetSolutionStepValue(VELOCITY_Y, 0, v_y)
    node.SetSolutionStepValue(VELOCITY_X, 0, 0.3 * v_y)
    
    if(node.coordinates[0] == 0.0):
        node.Fix(TEMPERATURE)
        node.SetSolutionStepValue(TEMPERATURE, 0, 1.0)

    if(node.coordinates[1] == 1.0):
        node.Fix(TEMPERATURE)
        node.SetSolutionStepValue(TEMPERATURE, 0, 0.0)

#import solving_strategy
#import static_scheme
time_scheme = static_scheme.StaticScheme(model_part)

#import builder_and_solver
builder_and_solver = builder_and_solver.BuilderAndSolver(
    model_part, time_scheme)

strategy = solving_strategy.SolvingStrategy(
    model_part, time_scheme, builder_and_solver)

strategy.Initialize()

mesh_name = "ConvectionDiffusion"
gid_io = GidIO("convection_diffusion.mdpa", "convection_diffusion")

gid_io.WriteMesh(model_part, mesh_name)


#import plot_contour
dt = 0.1
nsteps = 20
for i in range(1,nsteps):
    time = i*dt
    model_part.CloneTimeStep(time)
    print("before  solve")
    strategy.Solve()
    print("after solve")
    
    gid_io.WriteNodalResults(TEMPERATURE, model_part.NodeIterators(), time)
    #outname = "fig"+str(time)+".png"
    #plot_contour.PlotContour(model_part.NodeIterators(), TEMPERATURE, outname )
    
    #for node in model_part.NodeIterators():
        #if(node.coordinates[1] < 0.00001):
            #print(node.GetSolutionStepValue(TEMPERATURE,0))

# strategy.SpyMatrix()
#import plot_contour
#plot_contour.PlotContour(model_part.NodeIterators(), TEMPERATURE)