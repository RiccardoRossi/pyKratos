from __future__ import print_function, absolute_import, division 
import sys
sys.path.append("..")
print(sys.path)

from numpy import *
from pyKratos import *

#example = "closed" 
example = "gravity"

# add variables to be allocated from the list in variables.py
solution_step_variables = [
    VELOCITY_X,
    VELOCITY_Y,
    DISPLACEMENT_X,
    DISPLACEMENT_Y,
    PRESSURE,
    IS_LAGRANGIAN,
    EXTERNAL_FORCE_X,
    EXTERNAL_FORCE_Y

]

if example == "gravity":
    property_list = {
        0: {VISCOSITY: 1.0, #will use property 0 for the fluid
            DENSITY: 1.0,
            BODY_FORCE_X: 0.0,
            BODY_FORCE_Y: -10.0,
            },
        1: {YOUNG_MODULUS: 200e6, #will use property 1 for the structure
            POISSON_RATIO: 0.3,
            DENSITY: 1.0,
            BODY_FORCE_X: -10.0,
            BODY_FORCE_Y: 0.0,
            }
        }
else:
    property_list = {
        0: {VISCOSITY: 1.0, #will use property 0 for the fluid
            DENSITY: 1.0,
            BODY_FORCE_X: 0.0,
            BODY_FORCE_Y: 0.0,
            },
        1: {YOUNG_MODULUS: 200e1, #will use property 1 for the structure
            POISSON_RATIO: 0.3,
            DENSITY: 1.0,
            BODY_FORCE_X: 0.0,
            BODY_FORCE_Y: 0.0,
            }
        }
        
#######################################################
#### first we construct a model part made of just fluid Elements 
#defining a 1 by 1 square
nx = 11
dx = 1.0/(nx-1)
ny = 11
dy = 1.0/(ny-1)
import generate_square
node_list, element_connectivities, face_connectivities = generate_square.GenerateSquare(nx, dx, ny, dy)

#import py_kratos
buffer_size = 3  # store current step and 2 in the past
model_part = ModelPart(buffer_size, solution_step_variables)
model_part.AddNodes(node_list)
model_part.AddProperties(property_list)
model_part.AddElements("stokes_element_2d", element_connectivities)

n_elements_original = len(model_part.Elements)
print("original elements",n_elements_original)

#######################################################
#### now we generate a new model part, in which we exclude 
####the last column of elements to the righ
fluid_model_part = ModelPart(buffer_size, solution_step_variables)
fluid_model_part.Properties = model_part.Properties #copy properties from the original model part

#do not include the rightmost nodes
right_node_coordinates = (nx-1)*dx-0.00001

for node in model_part.NodeIterators():
    if node.coordinates[0] < right_node_coordinates:
        fluid_model_part.Nodes.update( {node.Id: node} )

#do not include elements which have a node on the rightmost boundary
for elem in model_part.ElementIterators():
    avoid_including = False
    for node in elem.geometry:
        if node.coordinates[0] > right_node_coordinates:
            avoid_including = True
    if avoid_including == False:
        fluid_model_part.Elements.update( {elem.Id: elem})

n_fluid_elements = len(fluid_model_part.Elements)
print("fluid elements",n_fluid_elements)

#######################################################
#generate a "structure_model_part" only including elements
#in the rightmost column
structure_model_part = ModelPart(buffer_size, solution_step_variables)
structure_model_part.Properties = model_part.Properties #copy properties from the original model part

#do not include elements which have a node on the rightmost boundary
for elem in model_part.ElementIterators():
    structural = False
    for node in elem.geometry:
        if node.coordinates[0] > right_node_coordinates:
            structural = True
    if structural == True:
        #here we need to create a "structural" element using the geometry of the existing one
        #note that we replace the property by propert 1
        import triangular_elastic_element
        structure_model_part.Elements.update( {elem.Id: triangular_elastic_element.ElasticityElement(elem.Id, structure_model_part.Properties[1], elem.geometry) } )
                 
##now add all of the elements taking them out of the elements
for elem in structure_model_part.ElementIterators():
    for node in elem.geometry:
        structure_model_part.Nodes.update( { node.Id: node } )

n_structure_elements = len(structure_model_part.Elements)
print("structure elements",n_structure_elements)
     
        
########################################################
#erase the original model part
model_part = 0

########################################################
#finally create a model part that contains everything 
fsi_model_part = ModelPart(buffer_size, solution_step_variables)
fsi_model_part.Properties = fluid_model_part.Properties #copy properties from the original model part

for elem in fluid_model_part.ElementIterators():
    fsi_model_part.Elements.update( {elem.Id:  elem} )
for elem in structure_model_part.ElementIterators():
    fsi_model_part.Elements.update( {elem.Id:  elem} )
    

for node in fluid_model_part.NodeIterators():
    fsi_model_part.Nodes.update( {node.Id:  node} )
for node in structure_model_part.NodeIterators():
    fsi_model_part.Nodes.update( {node.Id:  node} )
    
n_fsi_elements = len(fsi_model_part.Elements)
print("fsi elements",n_fsi_elements)





###now apply some boundary conditions ON THE FSI_MODEL_PART



bottom_nodes = []
left_side_nodes = []
right_side_nodes = []
top_nodes = []

for node in fsi_model_part.NodeIterators():
    if(node.coordinates[0] <= 0.000001):
        left_side_nodes.append(node)
    elif( node.coordinates[0] >(nx-1)*dx-0.00001):
        right_side_nodes.append(node)
    elif(node.coordinates[1] < 0.00001):
        bottom_nodes.append(node)
    elif(node.coordinates[1] >(ny-1)*dy-0.0001):
        top_nodes.append(node)
        

if example == "gravity":
    for node in bottom_nodes:
        node.Fix(VELOCITY_X)
        node.Fix(VELOCITY_Y)
        
    for node in left_side_nodes:
        node.Fix(VELOCITY_X)
        node.Fix(VELOCITY_Y)
        
    for node in fsi_model_part.NodeIterators():
        if(node.coordinates[0] == 1.0 and node.coordinates[1] == 1.0):
            top_right_node = node
elif example=="closed":
    for node in fsi_model_part.NodeIterators():
        if(node.coordinates[0] == 0.0):
            node.Fix(VELOCITY_X)
            node.Fix(VELOCITY_Y)
            node.SetSolutionStepValue(VELOCITY_X,0,1.0)
                
        if(node.coordinates[1] == (ny-1)*dy or node.coordinates[1] == 0.0):
            node.Fix(VELOCITY_X)
            node.Fix(VELOCITY_Y)
            node.SetSolutionStepValue(VELOCITY_X,0,0.0) ##this overwrites the one before on corner nodes


import gear_scheme
time_scheme = gear_scheme.GearScheme(fsi_model_part)

builder_and_solver = builder_and_solver.BuilderAndSolver(
    fsi_model_part, time_scheme)

strategy = solving_strategy.SolvingStrategy(
    fsi_model_part, time_scheme, builder_and_solver)

strategy.Initialize()

#import plot_contour

#outname = "structure"+".png"
#plot_contour.PlotContour(fsi_model_part.NodeIterators(), IS_LAGRANGIAN, outname )

GiDIO = gid_io.GidIO("inputfile.mdpa","gid_out")
GiDIO.WriteMesh(fluid_model_part,"outmesh")
t = []
disp_hist = []
dt = 0.1
nsteps = 100
for i in range(1,nsteps):
    time = i*dt
    t.append(time)
    fsi_model_part.CloneTimeStep(time)
    print("time: {:.6}".format(time))
    strategy.Solve()
    GiDIO.WriteNodalResults(PRESSURE,fluid_model_part.NodeIterators(), time)
    GiDIO.WriteNodalResults(DISPLACEMENT,fluid_model_part.NodeIterators(), time)
    GiDIO.WriteNodalResults(VELOCITY,fluid_model_part.NodeIterators(), time)
        
#    outname = "press"+str(time)+".png"
#    plot_contour.PlotContour(fluid_model_part.NodeIterators(), PRESSURE, outname )

 #   outname = "vx"+str(time)+".png"
#    plot_contour.PlotContour(fsi_model_part.NodeIterators(), VELOCITY_X, outname )

#    outname = "vy"+str(time)+".png"
#    plot_contour.PlotContour(fsi_model_part.NodeIterators(), VELOCITY_Y, outname )
    
#    if example=="gravity":
#        disp_hist = top_right_node.GetSolutionStepValue(DISPLACEMENT_X,0)
        
#import matplotlib.pyplot as plt
#plt.plot(t,d)
#plt.show()       
    
    #for node in model_part.NodeIterators():
        #if(node.coordinates[1] < 0.00001):
            #print(node.GetSolutionStepValue(TEMPERATURE,0))

# strategy.SpyMatrix()
#import plot_contour
#plot_contour.PlotContour(model_part.NodeIterators(), TEMPERATURE)