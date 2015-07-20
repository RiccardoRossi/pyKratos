from __future__ import print_function, absolute_import, division 
import sys
sys.path.append("..")
from numpy import *
from pyKratos import *

# add variables to be allocated from the list in variables.py
solution_step_variables = [
    TEMPERATURE,
    VOLUME_HEAT_SOURCE
]

property_list = {
    0: {CONDUCTIVITY: 0.001,
        DENSITY: 1000.0}
}

#generate a 5 by 5 square to play with
nx = 5
dx = 1.0/(nx-1)
ny = 5
dy = 1.0/(ny-1)
import generate_square
node_list, element_connectivities, face_connectivities = generate_square.GenerateSquare(nx, dx, ny, dy)

#import py_kratos

buffer_size = 3  # store current step and 2 in the past
model_part = ModelPart(buffer_size, solution_step_variables)
model_part.AddNodes(node_list)
model_part.AddProperties(property_list)
model_part.AddElements(
    "conv_diff_element2D", element_connectivities)
model_part.AddConditions(
    "temperature_neumann_face_condition_2d", face_connectivities)
gid_io = GidIO("convection_diffusion.mdpa", "convection_diffusion")

#model_part is a container of containers
model_part.Nodes #this is the list of nodes in the model
model_part.Elements #this is the list of Elements in the volume
model_part.Conditions #this is the list of "Conditions" used to apply Boundary conditions which require integration on the boundary
model_part.ProcessInfo #a dictionary used to store data and pass them to control elemental behaviour

#make a loop on all the nodes in the model and set the TEMPERATURE to a value
for node in model_part.NodeIterators():
    
    #node.coordinates contains the nodal position in space
    d = node.coordinates[0]**2 + node.coordinates[1]**2 - 0.75
    
    #here we assign the value of d to the node
    node.SetSolutionStepValue(TEMPERATURE, 0, d    )
    
    print("we are now considering the node with Id: ",node.Id)
    
    #we can verify the content by doing
    print("node now contains the following value of TEMPERATURE: ",node.GetSolutionStepValue(TEMPERATURE,0) )
    
    #at the previous time step it was zero...
    print("node at the previous step conainted the following value of TEMPERATURE: ",node.GetSolutionStepValue(TEMPERATURE,1 ) ) #note the ,1 to go 1 step in the past
    
    #i can check if a DOF of the node is fixed by doing
    print("is TEMPERATURE fixed for the node? ",node.IsFixed(TEMPERATURE) )
    
    #if i want to fix it i could do
    node.Fix(TEMPERATURE)
    
    #and to free it
    node.Free(TEMPERATURE)    
    
#i can also iterate over all the elements in the model
for elem in model_part.ElementIterators():
    #every element has an unique Id
    print("element Id", elem.Id)
    
    #an element contains a geometry which in turn can be iterated for its nodes
    geom = elem.geometry
    for node in geom:
        print("    Id = ",node.Id)
        
    #it also has a "property" which is in common between a group of elements
    #prop can be used to store variables which are common to a groupd o f nodes
    prop = elem.prop
    print("my conductivity is", prop[CONDUCTIVITY] ) #note that i am using directly prop as a python dictionary (which is actually what it is!!)
    #note that if you change the value of 

    
#i can also iterate over all the elements in the model
for cond in model_part.ConditionIterators():
    #every element has an unique Id
    print("Condition Id", elem.Id)
    
    #an element contains a geometry which in turn can be iterated for its nodes
    geom = cond.geometry
    for node in geom:
        print("    Id = ",node.Id)
        
    #it also has a "property" which is in common between a group of elements
    #prop can be used to store variables which are common to a groupd o f nodes
    prop = cond.prop
    print("my conductivity is", prop[CONDUCTIVITY] ) #note that i am using directly prop as a python dictionary (which is actually what it is!!)
    #note that if you change the value of 
    
#NOTE: essentially all of the databases in pyKRatos are implemented on the top of python dictionaries. 
#if you try to extract from the dictionary a variable that is not there an error will be thrown
#to correct it you need first to add the variable to the dictionary

#model_part.Properties[0].update( {TEMPERATURE: 123.456} ) ##uncomment this line and the next one will work
model_part.Properties[0][TEMPERATURE] #this throws an error since the TEMPERATURE is not in PROPERTIES