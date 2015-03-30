from __future__ import print_function, absolute_import, division 
from .node import *
from .variables import *


class ModelPart:

    def __init__(self, buffer_size, solution_step_variables):
        self.Nodes = {}  # empty dictionary
        self.Properties = {}  # empty dictionary
        self.Elements = {}  # empty dictionary
        self.Conditions = {}  # empty dictionary
        self.buffer_size = buffer_size
        self.solution_step_variables = solution_step_variables
        
        self.ProcessInfo = {TIME: 0.0, DELTA_TIME: 0.0}  # empty dictionary


    # function to access nodes as a list
    def NodeIterators(self):
        return list(self.Nodes.values())

    def ElementIterators(self):
        return list(self.Elements.values())
    
    def ConditionIterators(self):
        return list(self.Conditions.values())
    
    def PropertyIterators(self):
        return list(self.Properties.values())

    def CloneTimeStep(self, time):
        for node in self.NodeIterators():
            node.AdvanceInTime()
        old_time = self.ProcessInfo[TIME]
        self.ProcessInfo[TIME] = time
        self.ProcessInfo[DELTA_TIME] = time-old_time

    # function to create a list of nodes and give it to the model part
    def AddNodes(self, dict_of_nodes):
        for node_id, coords in list(dict_of_nodes.items()):
            if node_id in list(self.Nodes.keys()):
                error_string = "trying to add a node already existing with id =" + \
                    str(node_id)
                raise Exception(error_string)
            else:
                node = Node(node_id, coords)
                node.SetBufferSize(self.buffer_size)
                for var in self.solution_step_variables:
                    node.AddVariable(var)

                self.Nodes.update({node_id: node})

    def AddNode(self, node_id, coordinates):
            if node_id in list(self.Nodes.keys()):
                error_string = "trying to add a node already existing with id =" + \
                    str(node_id)
                raise Exception(error_string)
            else:
                node = Node(node_id, coordinates)
                node.SetBufferSize(self.buffer_size)
                for var in self.solution_step_variables:
                    node.AddVariable(var)

                self.Nodes.update({node_id: node})

    # add properties
    def AddProperties(self, dict_of_properties):
        self.Properties.update(dict_of_properties)

    def AddElements(self, element_name, dict_of_elements):
        el = __import__(element_name)

        for el_id, tmp in list(dict_of_elements.items()):
            # get property we assign to the element
            prop = self.Properties[tmp[0]]  # obtain the property pointer

            # get nodes that form the element
            connectivity = tmp[1]
            el_nodes = []
            for i in connectivity:
                el_nodes.append(self.Nodes[i])

            # create the element
            elem = el.Create(el_id, prop, el_nodes)

            self.Elements.update({el_id: elem})

    def AddConditions(self, condition_name, dict_of_conditions):
        cond = __import__(condition_name)

        for el_id, tmp in list(dict_of_conditions.items()):
            # get property we assign to the element
            prop = self.Properties[tmp[0]]  # obtain the property pointer

            # get nodes that form the element
            connectivity = tmp[1]
            el_nodes = []
            
            for i in connectivity:
                el_nodes.append(self.Nodes[i])

            # create the element
            condition = cond.Create(el_id, prop, el_nodes)

            self.Conditions.update({el_id: condition})


    def AddNodalValues(self, nodal_values):
        for var, values in list(nodal_values.items()):
            for entry in values:
                node = self.Nodes[entry[0]]
                # apply values
                node.SetSolutionStepValue(var, 0, entry[2])

                # fix if required
                if entry[1] == True:
                    node.Fix(var)

    def WriteMesh(self):
        pass

    def Check(self):
        pass

    def __str__(self):
        return "ModelPart:\n    Number of Nodes: {0}\n    Nunber of Elements: {1}".format(len(self.Nodes), len(self.Elements))

