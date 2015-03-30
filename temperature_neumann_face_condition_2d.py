from __future__ import print_function, absolute_import, division 
import math
from pyKratos import *
from numpy import *
from scipy import linalg


def Create(Id, prop, list_of_nodes):
    geom = line2d.Line2D(list_of_nodes)
    return TemperatureNeumannFaceCondition(Id, prop, geom)


class TemperatureNeumannFaceCondition:
    
    integration_order = 2  # this is like a c "static variable" one for all of the objects of this type
    include_dynamics = True

    def __init__(self, Id, prop, geometry):
        self.Id = Id
        self.prop = prop
        self.geometry = geometry

    def GetDofsPerNode(self):
        return 1

    def GetScalarValueOnGauss(self, var, N,step=0):
        value = 0.0
        for i in range(0, self.geometry.GetNumberOfNodes() ):
            value += N[i] * self.geometry[i].GetSolutionStepValue(var, step)
        return value
                
    def CalculateLocalSystem(self,ProcessInfo):
        order = self.integration_order
        nnodes = self.geometry.GetNumberOfNodes()
        dofs_per_node = self.GetDofsPerNode()
        mat_size = nnodes*dofs_per_node
        
        [Ns, derivatives, weights] = self.geometry.ShapeFunctions(order)
        number_of_gauss = len(Ns)
        
        RHS = zeros(mat_size)  # no external forces so far
        LHS = zeros((mat_size,mat_size))

        #integrate external forces to the RHS
        for gauss in range(0, number_of_gauss):
            weight = weights[gauss]
            N = Ns[gauss]
            
            #NOTE: external force is given "per unit area"
            fext = self.GetScalarValueOnGauss(FACE_HEAT_FLUX,N)
            #print(fext)
            
            for i in range(0,nnodes):               
                RHS[i] += weight*N[i]*fext

        return [LHS, RHS]


    # this function returns a list with the node and unkowns to be solved for
    def GetDofList(self):
        unknowns = []
        unknowns.append(Dof(self.geometry[0], TEMPERATURE))
        unknowns.append(Dof(self.geometry[1], TEMPERATURE))

        return unknowns

    def EquationId(self):
        equation_ids = []
        equation_ids.append(self.geometry[0].EquationId(TEMPERATURE))
        equation_ids.append(self.geometry[1].EquationId(TEMPERATURE))
        return equation_ids

    def GetValues(self, step=0):
        values = zeros(self.GetDofsPerNode()*self.geometry.GetNumberOfNodes())
        values[0] = self.geometry[0].GetSolutionStepValue(TEMPERATURE, step)
        values[1] = self.geometry[1].GetSolutionStepValue(TEMPERATURE, step)

        return values


