from __future__ import print_function, absolute_import, division 
from numpy import *
from pyKratos import *


class GearScheme:
    '''this is an INCREMENTAL scheme, designed to handle a correction of a predicted value
    it is designed so to handle both the linear and non-linear case
    the fundamental idea is that a second order backward differenci formula is used both for the velocity and 
    for the displacements (this last one, if needed)
    If we call 
    vn0 the velocity at the current step
    vn1 the velocity one step in the past
    vn2 the velocity two steps in the past
    
    dn0 the displacemnt at the current step
    dn1 the displacemnt one step in the past
    dn2 the displacemnt two steps in the past
    
    than we assume
    d(vn0)/dt |tn0 = c0*vn0 + c1*vn1 + c2*vn2 
    d(dn0)/dt |tn0 = c0*dn0 + c1*dn1 + c2*dn2 
    
    with
    c0 = 1.5/dt #3/(2dt)
    c1 = -2.0/dt
    c2 = 0.5/dt
    
    this scheme is designed to solve a system of the type
    M*d(vn0)/dt + D*vn0 + K*dn0 = fext
    
    the scheme provides a "guess" (prediction) of the values of vn0,dn0
    on the basis of this elements are expected to give back
    
    RHS = fext - M*d(vn0)/dt - D*vn0 - K*dn0
    
    and 
    
    LHS = d(RHS)/d(vn0) = c0*M + D + (1/c0)*K
    NOTE: this implies that elements are expected to be written in terms 
    of unknown VELOCITIES
    
    NOTE 2: only nodes marked as "IS_LAGRANGIAN" (node.GetSolutionStepValue(IS_LAGRANGIAN)==True)
    ##will be moved
    '''

    def __init__(self, model_part):
        self.model_part = model_part

    def CalculateLocalSystem(self, elem):
        return elem.CalculateLocalSystem(self.model_part.ProcessInfo)

    def EquationId(self, elem):
        return elem.EquationId()
    
    #compute the BDF coefficients c0, c1, c2 to be used in the computations
    #such coefficients are stored in a vector BDFvector = (c0,c1,c2)
    def ComputeBDFCoefficients(self):
        dt = self.model_part.ProcessInfo[DELTA_TIME]
        if(dt == 0):
            raise Exception("Dt can not be zero!!")
        c0 = 1.5/dt #3/(2dt)
        c1 = -2.0/dt
        c2 = 0.5/dt
        self.model_part.ProcessInfo.update( {BDF_COEFFICIENTS: array( [c0,c1,c2] ) } )

    def Predict(self):
        self.ComputeBDFCoefficients()
        
        #do a prediction of the velocity and of the displacement at  at the next step
        for node in self.model_part.NodeIterators():
            if node.IsFixed(VELOCITY_X)==False:
                #predict component x of velocity
                vn1 = node.GetSolutionStepValue(VELOCITY_X,1)
                vn2 = node.GetSolutionStepValue(VELOCITY_X,2)
                vn0 = 2.0*vn1-vn2
                node.SetSolutionStepValue(VELOCITY_X,0,vn0)
            
            if node.IsFixed(VELOCITY_Y)==False:
                #predict component y of velocity
                vn1 = node.GetSolutionStepValue(VELOCITY_Y,1)
                vn2 = node.GetSolutionStepValue(VELOCITY_Y,2)
                vn0 = 2.0*vn1-vn2
                node.SetSolutionStepValue(VELOCITY_Y,0,vn0)
            
        self.UpdateDisplacements(self.model_part.NodeIterators())
        
    def UpdateDisplacements(self,nodes):
        coeffs =  self.model_part.ProcessInfo[BDF_COEFFICIENTS]
        c0 = coeffs[0]
        c1 = coeffs[1]
        c2 = coeffs[2]
        
        for node in nodes:
            if(node.GetSolutionStepValue(IS_LAGRANGIAN,0) == True):
                #update component x of displacement
                dn1 = node.GetSolutionStepValue(DISPLACEMENT_X,1)
                dn2 = node.GetSolutionStepValue(DISPLACEMENT_X,2)
                vn0 = node.GetSolutionStepValue(VELOCITY_X,0)
                dn0 = 1.0/c0*(vn0 - c1*dn1-c2*dn2)
                node.SetSolutionStepValue(DISPLACEMENT_X,0,dn0)
                
                #update component y of displacement
                dn1 = node.GetSolutionStepValue(DISPLACEMENT_Y,1)
                dn2 = node.GetSolutionStepValue(DISPLACEMENT_Y,2)
                vn0 = node.GetSolutionStepValue(VELOCITY_Y,0)
                dn0 = 1.0/c0*(vn0 - c1*dn1-c2*dn2)
                node.SetSolutionStepValue(DISPLACEMENT_Y,0,dn0)
        

    def Update(self, dofset, dx):

        # here we do x = old_value + dx
        #this is expected to update VELOCITY (and eventually other quantities)
        #so after this loop the value VELOCITY has the "final" value (of the iteration)
        for dof, value in zip(dofset, dx):
            if dof.IsFixed() == False:
                old_value = dof.GetValue()
                dof.SetValue(old_value + value)
        
        #update displacements following BDF2 formula
        self.UpdateDisplacements(self.model_part.NodeIterators())