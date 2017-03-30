from __future__ import print_function, absolute_import, division 
import math
from numpy import *


class Point2D:

    def __init__(self, node_list):
        if(len(node_list) != 1):
            raise Exception("wrong number of nodes! should be 1!!")
        self.nodes = node_list

        for node in self.nodes:
            if(node.Id < 0):
                raise Exception("node with Id lesser than 0 found")

    # def Nodes(self):
        # return self.nodes

    def __getitem__(self, key):
        return self.nodes[key]
    
    def GetNumberOfNodes(self):
        return 1

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
        
        Ncontainer.append(array([1.0])) #does not make sense
        weights = [1.0]
        derivatives = [] #does not make sense

        return [Ncontainer, derivatives, weights]
