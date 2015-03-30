from __future__ import print_function, absolute_import, division 
from numpy import *
from .variables import *


class StaticScheme:

    '''this is an INCREMENTAL scheme, designed to handle a correction of a predicted value
    it is designed so to handle both the linear and non-linear case'''

    def __init__(self, model_part):
        self.model_part = model_part

    def CalculateLocalSystem(self, elem):
        return elem.CalculateLocalSystem(self.model_part.ProcessInfo)

    def EquationId(self, elem):
        return elem.EquationId()

    def Predict(self):
        # eventually do a prediction of the solution. here do nothing
        pass

    def Update(self, dofset, dx):

        # here we do x = old_value + dx
        for dof, value in zip(dofset, dx):
            if dof.IsFixed() == False:
                old_value = dof.GetValue()
                dof.SetValue(old_value + value)
