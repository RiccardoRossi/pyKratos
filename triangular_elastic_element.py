from __future__ import print_function, absolute_import, division 
import math
from pyKratos import *
from numpy import *
from scipy import linalg


def Create(Id, prop, list_of_nodes):
    geom = triangle.Triangle(list_of_nodes)
    return ElasticityElement(Id, prop, geom)


class ElasticityElement:
    integration_order = 2  # this is like a c "static variable" one for all of the objects of this type

    def __init__(self, Id, prop, geometry):
        self.Id = Id
        self.prop = prop
        self.geometry = geometry
        
        #mark as lagrangian all of the nodes in the element
        for node in self.geometry:
            node.SetSolutionStepValue(IS_LAGRANGIAN,0,True)

    def GetDofsPerNode(self):
        return 2



    

    #this _auxiliary function computes the stiffness contribution
    def _ComputeStiffnessContribution(self,ProcessInfo):
        order = self.integration_order
        [Ns, derivatives, weights] = self.geometry.ShapeFunctions(order)

        RHSstatic = zeros(6)  # no external forces so far
        K = zeros((6, 6))
        
        body_force = zeros(2)
        body_force[0] = self.prop[BODY_FORCE_X]
        body_force[1] = self.prop[BODY_FORCE_Y]
        density = self.prop[DENSITY]
        
        nnodes = self.geometry.GetNumberOfNodes()

        number_of_gauss = len(Ns)
        for gauss in range(0, number_of_gauss):
            N = Ns[gauss]
            DN_DX = derivatives[gauss]
            
            #computation of external forces
            for i in range(0,nnodes):
                for k in range(0,2):
                    RHSstatic[i*2+k] += weights[gauss]*density*N[i]*body_force[k]
            
            ##computation of internal forces
            C = self.ComputeElasticityTensor()
            B = self.ComputeB(DN_DX)

            tmp = dot(C, B)
            K += weights[gauss] * dot(B.transpose(), tmp)

        # compute RESIDUAL for the RHS
        # since the problem is linear this can be done as
        # RHS = fext - LHS*values
        # note that this is done out of the integration loop!
        values = self._GetAllValuesVector(DISPLACEMENT_X,DISPLACEMENT_Y)  # get values of unknown at the nodes
        RHSstatic -= dot(K, values)

        return [K, RHSstatic]
    
    #this _auxiliary function computes the mass contribution
    def _ComputeMassContribution(self,ProcessInfo):
        density = self.prop[DENSITY]

        RHSinertia = zeros(6)
        M = zeros((6, 6))
        
        order = self.integration_order
        if(order == 1):
            raise Exception("exact mass matrix can not be computer with an integration order of 1. Either switch to order 2 or do mass lumping")
        [Ns, derivatives, weights] = self.geometry.ShapeFunctions(order)

        ##compute external forces
        

        ##compute stiffness
        number_of_gauss = len(Ns)
        for gauss in range(0, number_of_gauss):

            N = Ns[gauss]
            weight = weights[gauss]

            for i in range(0, 3):
                for j in range(0, 3):
                    M[i * 2, j * 2] += N[i] * N[j] * density * weight
                    M[i * 2 + 1, j * 2 + 1] += N[i] * N[j] * density * weight

        acc = self._GetAccelerationsVector(ProcessInfo)
        RHSinertia -= dot(M,acc)

        return [M, RHSinertia]
    
    def CalculateLocalSystem(self,ProcessInfo):
        K, RHSstatic = self._ComputeStiffnessContribution(ProcessInfo)
        M, RHSinertia = self._ComputeMassContribution(ProcessInfo)
        
        #get BDF coefficients to correctly compute the LHS
        coeffs =  ProcessInfo[BDF_COEFFICIENTS]
        c0 = coeffs[0]
        c1 = coeffs[1]
        c2 = coeffs[2]

        
        RHS = RHSstatic + RHSinertia
        LHS = c0*M + 1.0/c0*K
        
        return [LHS,RHS]

    def ComputeB(self, DN_DX):

        B = zeros((3, 6))

        for i in range(0, 3):
            index = i * 2
            B[0, index + 0] = DN_DX[i, 0]
            B[0, index + 1] = 0
            B[1, index + 0] = 0
            B[1, index + 1] = DN_DX[i, 1]
            B[2, index + 0] = DN_DX[i, 1]
            B[2, index + 1] = DN_DX[i, 0]

        return B

    def ComputeElasticityTensor(self):
        C = zeros((3, 3))
        E = self.prop[YOUNG_MODULUS]
        NU = self.prop[POISSON_RATIO]

        c1 = E / (1.00 - NU * NU)

        c2 = NU * c1

        c3 = E / 2 / (1 + NU)

        C = [[c1, c2, 0.0], [c2, c1, 0.0], [0.0, 0.0, c3]]

        return C

    # this function returns a list with the node and unkowns to be solved for
    def GetDofList(self):
        unknowns = []
        unknowns.append(Dof(self.geometry[0], VELOCITY_X))
        unknowns.append(Dof(self.geometry[0], VELOCITY_Y))
        unknowns.append(Dof(self.geometry[1], VELOCITY_X))
        unknowns.append(Dof(self.geometry[1], VELOCITY_Y))
        unknowns.append(Dof(self.geometry[2], VELOCITY_X))
        unknowns.append(Dof(self.geometry[2], VELOCITY_Y))
        return unknowns

    def EquationId(self):
        equation_ids = []
        equation_ids.append(self.geometry[0].EquationId(VELOCITY_X))
        equation_ids.append(self.geometry[0].EquationId(VELOCITY_Y))
        equation_ids.append(self.geometry[1].EquationId(VELOCITY_X))
        equation_ids.append(self.geometry[1].EquationId(VELOCITY_Y))
        equation_ids.append(self.geometry[2].EquationId(VELOCITY_X))
        equation_ids.append(self.geometry[2].EquationId(VELOCITY_Y))
        return equation_ids

    def _GetAllValuesVector(self, var_x,var_y, step=0):
        values = zeros(6)
        values[0] = self.geometry[0].GetSolutionStepValue(var_x, step)
        values[1] = self.geometry[0].GetSolutionStepValue(var_y, step)
        values[2] = self.geometry[1].GetSolutionStepValue(var_x, step)
        values[3] = self.geometry[1].GetSolutionStepValue(var_y, step)
        values[4] = self.geometry[2].GetSolutionStepValue(var_x, step)
        values[5] = self.geometry[2].GetSolutionStepValue(var_y, step)
        return values


    def _GetAccelerationsVector(self,ProcessInfo):
        values = zeros(6)
        coeffs =  ProcessInfo[BDF_COEFFICIENTS]
        c0 = coeffs[0]
        c1 = coeffs[1]
        c2 = coeffs[2]
                
        v0 = self._GetAllValuesVector(VELOCITY_X,VELOCITY_Y,0)
        v1 = self._GetAllValuesVector(VELOCITY_X,VELOCITY_Y,1)
        v2 = self._GetAllValuesVector(VELOCITY_X,VELOCITY_Y,2)
               
        return c0*v0+c1*v1+c2*v2
