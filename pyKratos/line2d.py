from __future__ import print_function, absolute_import, division 
import math
from numpy import *


class Line2D:

    def __init__(self, node_list):
        if(len(node_list) != 2):
            raise Exception("wrong number of nodes! should be 2!!")
        self.nodes = node_list

        for node in self.nodes:
            if(node.Id < 0):
                raise Exception("node with Id lesser than 0 found")

    # def Nodes(self):
        # return self.nodes

    def __getitem__(self, key):
        return self.nodes[key]
    
    def GetNumberOfNodes(self):
        return 2

    def ShapeFunctions(self, order=1):
        '''this function provides the shape function values, derivatives and integration_weight'''
        '''at the location of the gauss points. Order of integration is controlled'''
        '''by the optional parameter "order".'''
        '''N[gauss][i] contains the shape function of node i computed at the position of "gauss" '''
        '''derivatives[gauss][i,k] contains the derivative of node i, component k at the position of gauss '''
        '''weights[gauss] includes the integration weights, including the det of the jacobian, to be used '''
        '''at the gauss point'''
        derivatives = []
        weights = []
        Ncontainer = []

        x10 = self.nodes[1].coordinates[0] - self.nodes[0].coordinates[0]
        y10 = self.nodes[1].coordinates[1] - self.nodes[0].coordinates[1]

        lenght = math.sqrt(x10**2 + y10**2)
        
        if(order == 1):  # give back 1 single integration point
            one_half = 1.0 / 2.0
            Ncontainer = [array([one_half, one_half])]

            weights = [lenght]
            derivatives = [] #it is a 1d line in 2d space ... global derivatives are not defined

        elif(order == 2):  # gives back 3 integration points
            aux = math.sqrt(1.0 / 3.0) / 2.0 + 0.5

            Ncontainer.append(array([aux, 1.0-aux]))
            Ncontainer.append(array([1.0-aux, aux]))

            weights = [0.5*lenght, 0.5*lenght]

            derivatives = [] #it is a 1d line in 2d space ... global derivatives are not defined

        else:
            raise Exception("integration order not implemented")

        return [Ncontainer, derivatives, weights]
