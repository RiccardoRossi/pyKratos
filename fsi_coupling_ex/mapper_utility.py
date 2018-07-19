from numpy import *
from pyKratos import *

class Mapper:
    def __init__(self, mp1, mp2):

        self.interface1 = []
        self.interface2 = []

        for node in mp1.NodeIterators():
            if(node.GetSolutionStepValue(FSI_INTERFACE,0) == 1.0):
                self.interface1.append(node)

        for node in self.interface1:
            for other in mp2.NodeIterators():
                d = node.coordinates - other.coordinates
                if( linalg.norm(d) < 1e-9 ):
                    self.interface2.append(other)
                    break

        if(len(self.interface1) != len(self.interface2)):
            raise Exception("different number of nodes in interface 1 and 2")

    def Map_1to2(self, origin_var, destination_var, factor=1.0):
        for i1,i2 in zip(self.interface1, self.interface2):
            i2.SetSolutionStepValue(destination_var,0, factor*i1.GetSolutionStepValue(origin_var,0))

    def Map_2to1(self, origin_var, destination_var, factor=1.0):
        for i1,i2 in zip(self.interface1, self.interface2):
            i1.SetSolutionStepValue(destination_var,0, factor*i2.GetSolutionStepValue(origin_var,0))

