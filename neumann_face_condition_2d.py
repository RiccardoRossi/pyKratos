from __future__ import print_function, absolute_import, division 
import math
from pyKratos import *
from numpy import *
from scipy import linalg


def Create(Id, prop, list_of_nodes):
    geom = line2d.Line2D(list_of_nodes)
    return NeumannFaceCondition(Id, prop, geom)


class NeumannFaceCondition:
    #this elements constructs a stiffness matrix which mixes velocities and pressures
    #for each pair of nodes I and J this matrix has a 3*3 subblock ordered as
    #  | Kvv Kvp | 
    #  | Kpv Kpp |3x3
    
    integration_order = 2  # this is like a c "static variable" one for all of the objects of this type
    include_dynamics = True

    def __init__(self, Id, prop, geometry):
        self.Id = Id
        self.prop = prop
        self.geometry = geometry

    def GetDofsPerNode(self):
        return 2

    def GetVectorValueOnGauss(self, var_x, var_y, N,step=0):
        value = zeros(2)
        for i in range(0, self.geometry.GetNumberOfNodes() ):
            value[0] += N[i] * self.geometry[i].GetSolutionStepValue(var_x, step)
            value[1] += N[i] * self.geometry[i].GetSolutionStepValue(var_y, step)
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
            fext = self.GetVectorValueOnGauss(EXTERNAL_FORCE_X,EXTERNAL_FORCE_Y,N)
            #print(fext)
            
            for i in range(0,nnodes):               
                for k in range(0,2):
                    RHS[i*2+k] += weight*N[i]*fext[k]
            
            
            
            if(self.include_dynamics == True):
                coeffs =  ProcessInfo[BDF_COEFFICIENTS]
                c0 = coeffs[0]
                c1 = coeffs[1]
                c2 = coeffs[2]
                
                #if an external density is given add it to the sys 
                if(EXTERNAL_UNIT_DENSITY in self.prop):
                    ext_density = self.prop[EXTERNAL_UNIT_DENSITY]
                    
                    #part of "acc" to the RHS
                    v1gauss = self.GetVectorValueOnGauss(VELOCITY_X,VELOCITY_Y,N,1) #old step
                    v2gauss = self.GetVectorValueOnGauss(VELOCITY_X,VELOCITY_Y,N,2) #two steps ago
                    arhs = c1*v1gauss + c2*v2gauss 
                    #print(arhs)
                    for i in range(0,nnodes):
                        RHS[i*2   ] -= (weight*ext_density*N[i])*arhs[0]
                        RHS[i*2 +1] -= (weight*ext_density*N[i])*arhs[1]
        
                    #part of "acc" to the LHS
                    for i in range(0, 3):
                        for j in range(0, 3):
                            tmp = (weight*ext_density*c0)*N[i]*N[j] 
                            LHS[i*2  , j*2  ] += tmp
                            LHS[i*2+1, j*2+1] += tmp

                    
        # compute RESIDUAL for the RHS
        # since the problem is LINEAR this can be done as
        # RHS = fext - LHS*values
        # note that this is done out of the integration loop!
        values = self.GetValues()  # get values of unknown at the nodes
        RHS -= dot(LHS, values)
        return [LHS, RHS]

        return C

    # this function returns a list with the node and unkowns to be solved for
    def GetDofList(self):
        unknowns = []
        unknowns.append(Dof(self.geometry[0], VELOCITY_X))
        unknowns.append(Dof(self.geometry[0], VELOCITY_Y))
        unknowns.append(Dof(self.geometry[1], VELOCITY_X))
        unknowns.append(Dof(self.geometry[1], VELOCITY_Y))
        return unknowns

    def EquationId(self):
        equation_ids = []
        equation_ids.append(self.geometry[0].EquationId(VELOCITY_X))
        equation_ids.append(self.geometry[0].EquationId(VELOCITY_Y))
        equation_ids.append(self.geometry[1].EquationId(VELOCITY_X))
        equation_ids.append(self.geometry[1].EquationId(VELOCITY_Y))
        return equation_ids

    def GetValues(self, step=0):
        values = zeros(self.GetDofsPerNode()*self.geometry.GetNumberOfNodes())
        values[0] = self.geometry[0].GetSolutionStepValue(VELOCITY_X, step)
        values[1] = self.geometry[0].GetSolutionStepValue(VELOCITY_Y, step)
        values[2] = self.geometry[1].GetSolutionStepValue(VELOCITY_X, step)
        values[3] = self.geometry[1].GetSolutionStepValue(VELOCITY_Y, step)
        return values


