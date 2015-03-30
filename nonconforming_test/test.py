from __future__ import print_function, absolute_import, division 
import sys
sys.path.append("..")
from numpy import *
from pyKratos import *
# add variables to be allocated from the list in variables.py
solution_step_variables = [
    TEMPERATURE,
    VELOCITY_X,
    VELOCITY_Y,
    VOLUME_HEAT_SOURCE,
    FACE_HEAT_FLUX
]

property_list = {
    0: {CONDUCTIVITY: 1.0,
        DENSITY: 1.0e-3}
}

# 5  -  6  - 7
# |  /  |  / |
# 2  -  3  - 4
# |     /    |
# 0 -  -  -  1

node_list = {
    0: array([0.0, 0.0]),
    1: array([0.2, 0.0]),
    2: array([0.0, 0.1]),
    3: array([0.1, 0.1]),
    4: array([0.2, 0.1]),
    5: array([0.0, 0.2]),
    6: array([0.1, 0.2]),
    7: array([0.2, 0.2]),
    }

element_connectivities = {
    1: [0, [0, 1, 4]],
    2: [0, [0, 4, 2]],
    3: [0, [2, 3, 6]],
    4: [0, [2,6,5]],
    5: [0, [3,4,7]],
    6: [0, [3,7,6]],
}
condition_connectivities = {
    1: [0, [7,6]],
    2: [0, [6,5]],
}

nodal_values = {
    VELOCITY_X: [
        [0, True, 0.0],  # first column is Id of node, second col if fixity, third is imposed value
        [1, True, 0.0],
        [2, True, 0.0],
        [3, True, 0.0],
        [4, True, 0.0],
        [5, True, 0.0],
        [6, True, 0.0],
        [7, True, 0.0],
    ],
VELOCITY_Y: [
        [0, True, 0.0],  # first column is Id of node, second col if fixity, third is imposed value
        [1, True, 0.0],
        [2, True, 0.0],
        [3, True, 0.0],
        [4, True, 0.0],
        [5, True, 0.0],
        [6, True, 0.0],
        [7, True, 0.0],
],
#TEMPERATURE: [
        #[5, True, 10.0],  # first column is Id of node, second col if fixity, third is imposed value
        #[6, True, 10.0],
        #[7, True, 10.0],
#],
TEMPERATURE: [
        [0, True, 10.0],  # first column is Id of node, second col if fixity, third is imposed value
        [1, True, 10.0],
],
FACE_HEAT_FLUX: [
        [5, True, 10.0],  # first column is Id of node, second col if fixity, third is imposed value
        [6, True, 10.0],
        [7, True, 10.0],
],
}

#import py_kratos

buffer_size = 3  # store current step and 2 in the past
model_part = ModelPart(buffer_size, solution_step_variables)
model_part.AddNodes(node_list)
model_part.AddProperties(property_list)
model_part.AddElements(
    "conv_diff_element2D", element_connectivities)
model_part.AddConditions(
    "temperature_neumann_face_condition_2d", condition_connectivities)
model_part.AddNodalValues(nodal_values)
gid_io = GidIO("convection_diffusion.mdpa", "convection_diffusion")

# apply a velocity_y


time_scheme = static_scheme.StaticScheme(model_part)

builder_and_solver = builder_and_solver.BuilderAndSolver(
    model_part, time_scheme)

strategy = solving_strategy.SolvingStrategy(
    model_part, time_scheme, builder_and_solver)
strategy.builder_and_solver.use_sparse_matrices = False

strategy.Initialize()
mesh_name = "ConvectionDiffusion"
gid_io.WriteMesh(model_part, mesh_name)


#constraint_type
#0 - replace row
#1 - non symmetric
#2 - symmetric
#
constraint_type = 2

dt = 100000.0 #0.1
nsteps = 20
for i in range(1,nsteps):
    time = i*dt
    print("time: {:.6}".format(time))
    model_part.CloneTimeStep(time)
    #strategy.Solve()
    #strategy.SpyMatrix()
    
    A, dx, b = strategy.builder_and_solver.Build(strategy.A, strategy.dx, strategy.b)
    
    ##substituting line 3 by a constraint
    if(constraint_type == 0): #replace row
        A[3,:] = 0.0
        A[3,2]=-0.5
        A[3,3]=+1.0
        A[3,4]=-0.5
        b[3] = 0.0
    elif(constraint_type == 1): #adding a NON-SYMM constraint operator
        aux =  zeros(A.shape)
        for i in range(0, A.shape[0]):
            aux[i,i] = 1.0
        aux[3,2] = 0.5
        aux[3,3] = 0.0   
        aux[3,4] = 0.5

        tmp = A.dot(aux)    
        A = tmp.copy()
        A[3,:] = 0.0
        A[3,3] = 1.0
    elif(constraint_type == 2): #adding a NON-SYMM constraint operator    
        aux =  zeros(A.shape)
        for i in range(0, A.shape[0]):
            aux[i,i] = 1.0
        aux[3,3] = 0.0
        aux[3,2] = 0.5
        aux[3,4] = 0.5
        print(aux)    
        print(A)
        tmp = A.dot(aux)    
        A = (aux.transpose()).dot(tmp)
        b = (aux.transpose()).dot(b)    
        A[3,3] = 1.0
    
    #print(A)
    
    A, dx, b = strategy.builder_and_solver.ApplyDirichlet(A, dx, b)

    dx = linalg.solve(A, b)
    
    if(constraint_type != 0):
        dx = aux.dot(dx.copy())
            
    strategy.scheme.Update(strategy.builder_and_solver.dofset, dx)  
   
    gid_io.WriteNodalResults(TEMPERATURE, model_part.NodeIterators(), time)
    
    #for node in model_part.NodeIterators():
        #if(node.coordinates[1] < 0.00001):
            #print(node.GetSolutionStepValue(TEMPERATURE,0))

# strategy.SpyMatrix()
#import plot_contour
#plot_contour.PlotContour(model_part.NodeIterators(), TEMPERATURE)