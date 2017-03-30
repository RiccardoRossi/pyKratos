from __future__ import print_function, absolute_import, division 
import math
from pyKratos import *
from numpy import *
from scipy import linalg


def Create(Id, prop, list_of_nodes):
    geom = point2d.Point2D(list_of_nodes)
    return Point2DCondition(Id, prop, geom)

class Point2DCondition:

    def __init__(self, Id, prop, geometry):
        self.Id = Id
        self.prop = prop
        self.geometry = geometry

    def GetDofsPerNode(self):
        return 2

                
    def CalculateLocalSystem(self,ProcessInfo):
        fx = self.geometry[0].GetSolutionStepValue(EXTERNAL_FORCE_X,0) 
        fy = self.geometry[0].GetSolutionStepValue(EXTERNAL_FORCE_Y,0)   
        
        RHS = zeros(2)  
        LHS = zeros((2,2))
        
        RHS[0] = fx
        RHS[1] = fy

        return [LHS, RHS]

    # this function returns a list with the node and unkowns to be solved for
    def GetDofList(self):
        unknowns = []
        unknowns.append(Dof(self.geometry[0], VELOCITY_X))
        unknowns.append(Dof(self.geometry[0], VELOCITY_Y))
        return unknowns

    def EquationId(self):
        equation_ids = []
        equation_ids.append(self.geometry[0].EquationId(VELOCITY_X))
        equation_ids.append(self.geometry[0].EquationId(VELOCITY_Y))
        return equation_ids

    def GetValues(self, step=0):
        values = zeros(self.GetDofsPerNode()*self.geometry.GetNumberOfNodes())
        values[0] = self.geometry[0].GetSolutionStepValue(VELOCITY_X, step)
        values[1] = self.geometry[0].GetSolutionStepValue(VELOCITY_Y, step)
        return values


