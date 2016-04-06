# ----------------------------------------------------------------------
# author  : Martin Ruchti
# contact : martin.ruchti@tum.de
# ----------------------------------------------------------------------
from __future__ import print_function, absolute_import, division 
from numpy import *
from pyKratos import *

'''
    The implementation of this scheme follows the implementation of the residual based predictor corrector Bossak scheme 
    from Kratos Multiphysics.
    
    The Bossak scheme is a special version of the generalized alpha scheme (alpha_f = 0)
'''
class BossakScheme:

    def __init__(self, model_part, alpha):
        
        self.model_part         = model_part
        self.alphaBossak        = alpha
        self.gammaNewmark       = 0.5 - self.alphaBossak
        self.betaNewmark        = 0.25 * ( ( 1.0 - self.alphaBossak )**2 )
        
        #init parameters with 0.0
        self.ma0 = 0.0
        self.ma1 = 0.0
        self.ma2 = 0.0
        self.ma3 = 0.0
        self.ma4 = 0.0
        self.ma5 = 0.0
        self.mam = 0.0

    def CalculateLocalSystem(self, elem):
        return elem.CalculateLocalSystem(self.model_part.ProcessInfo)
        
    def CalculateLocalAdjointSystem(self, elem, dragForceDirection, stepRelative):
        if(stepRelative == 0):
            return elem.CalculateInitialLocalAdjointSystem(dragForceDirection, self.model_part.ProcessInfo)
        else:
            return elem.CalculateLocalAdjointSystem(dragForceDirection, self.model_part.ProcessInfo, stepRelative)
    
    def CalculateLocalAdjointDOTSystem(self, elem, dragForceDirection, stepRelative):
        if(stepRelative == 0):
            return elem.CalculateInitialLocalAdjointDOTSystem(dragForceDirection, self.model_part.ProcessInfo)
        else:
            return elem.CalculateLocalAdjointDOTSystem(dragForceDirection, self.model_part.ProcessInfo, stepRelative)

    def EquationId(self, elem):
        return elem.EquationId()
    
    #compute the integration coefficients
    def ComputeCoefficients(self):
        
        #retrieve current time step from process info
        dt = self.model_part.ProcessInfo[DELTA_TIME]
        
        #check that dt is not 0, because this would cause the inertia terms to
        #be inf.
        if(dt == 0):
            raise Exception("Dt can not be zero!!")
        
        #calculate coefficients for the newmark/alpha scheme
        self.ma0 = 1.0 / (self.gammaNewmark * dt);
        #self.ma1 = dt * self.betaNewmark / self.gammaNewmark;
        self.ma2 = (-1.0 + self.gammaNewmark) / self.gammaNewmark;
        self.ma3 = dt;
        self.ma4 = dt**2 *(-2.0 * self.betaNewmark + 1.0) / 2.0;
        self.ma5 = dt**2 * self.betaNewmark;
        #self.mam = (1.0 - self.alphaBossak) / (self.gammaNewmark * dt);
        
        #write coefficients to the process info to be accessed within the elements
        self.model_part.ProcessInfo.update( {BOSSAK_COEFFICIENTS: array( [self.alphaBossak,self.gammaNewmark] ) } )


    def Predict(self):
        
        self.ComputeCoefficients()
        #if nodes are not fixed, set the new value to the old one
        #this is the starting point for the iteration process
        for node in self.model_part.NodeIterators():
            if node.IsFixed(VELOCITY_X)==False:
                node.SetSolutionStepValue(VELOCITY_X,0,node.GetSolutionStepValue(VELOCITY_X,1))
            
            if node.IsFixed(VELOCITY_Y)==False:
                node.SetSolutionStepValue(VELOCITY_Y,0,node.GetSolutionStepValue(VELOCITY_Y,1))
                
            if node.IsFixed(PRESSURE)==False:
                node.SetSolutionStepValue(PRESSURE,0,node.GetSolutionStepValue(PRESSURE,1))
            
        #update displacement and accceleration
        self.UpdateAcceleration(self.model_part.NodeIterators())
        self.UpdateDisplacements(self.model_part.NodeIterators())
        
    def UpdateDisplacements(self,nodes):
        
        for node in nodes:
            if(node.GetSolutionStepValue(IS_LAGRANGIAN,0) == True):
                #update component x of displacement
                dn1 = node.GetSolutionStepValue(DISPLACEMENT_X,1)
                an0 = node.GetSolutionStepValue(ACCELERATION_X,0)
                an1 = node.GetSolutionStepValue(ACCELERATION_X,1)
                vn1 = node.GetSolutionStepValue(VELOCITY_X,1)
                dn0 = dn1 + self.ma3 * vn1 + self.ma4 * an1 + self.ma5 * an0
                node.SetSolutionStepValue(DISPLACEMENT_X,0,dn0)
                
                #update component y of displacement
                dn1 = node.GetSolutionStepValue(DISPLACEMENT_Y,1)
                an0 = node.GetSolutionStepValue(ACCELERATION_Y,0)
                an1 = node.GetSolutionStepValue(ACCELERATION_Y,1)
                vn1 = node.GetSolutionStepValue(VELOCITY_Y,1)
                dn0 = dn1 + self.ma3 * vn1 + self.ma4 * an1 + self.ma5 * an0
                node.SetSolutionStepValue(DISPLACEMENT_Y,0,dn0)
        
    def UpdateAcceleration(self,nodes):
        
        for node in nodes:
            #predict acceleration_x
            if node.IsFixed(ACCELERATION_X)==False:
                aOld            = node.GetSolutionStepValue(ACCELERATION_X,1)
                deltaVel        = node.GetSolutionStepValue(VELOCITY_X,0) - node.GetSolutionStepValue(VELOCITY_X,1)
                aNewX           = self.ma0 * deltaVel + self.ma2*aOld
                node.SetSolutionStepValue(ACCELERATION_X, 0, aNewX)
            #predict acceleration_y
            if node.IsFixed(ACCELERATION_Y)==False:
                aOld            = node.GetSolutionStepValue(ACCELERATION_Y,1)
                deltaVel        = node.GetSolutionStepValue(VELOCITY_Y,0) - node.GetSolutionStepValue(VELOCITY_Y,1)
                aNewY           = self.ma0 * deltaVel + self.ma2*aOld
                node.SetSolutionStepValue(ACCELERATION_Y, 0, aNewY)

    def Update(self, dofset, dx):
        # here we do x = old_value + dx
        #this is expected to update VELOCITY (and eventually other quantities)
        for dof, value in zip(dofset, dx):
            if dof.IsFixed() == False:
                old_value = dof.GetValue()
                dof.SetValue(old_value + value)
        
        #update the dependent values
        self.UpdateDisplacements(self.model_part.NodeIterators())
        self.UpdateAcceleration(self.model_part.NodeIterators())
