from __future__ import print_function, absolute_import, division 
class Node:

    def __init__(self, Id, coordinates):
        self.variables = []
        self.var_is_fixed = dict()
        self.equation_id = dict()
        self.coordinates = coordinates
        self.Id = Id

    def SetBufferSize(self, buffer_size):
        for i in range(0, buffer_size):
            self.variables.append(dict())

    def AddVariable(self, variable_name):
        for i in range(0, len(self.variables)):
            self.variables[i][variable_name] = 0
        self.var_is_fixed[variable_name] = False
        self.equation_id[variable_name] = -1

    def AdvanceInTime(self):
        for i in range(len(self.variables)-1,0,-1):
            for key in list(self.variables[i].keys()):
                self.variables[i][key] = self.variables[i - 1][key]

    def Fix(self, variable_name):
        self.var_is_fixed[variable_name] = True

    def Free(self, variable_name):
        self.var_is_fixed[variable_name] = False

    def IsFixed(self, variable_name):
        return self.var_is_fixed[variable_name]

    def EquationId(self, variable_name):
        return self.equation_id[variable_name]

    def GetSolutionStepValue(self, variable_name, step):
        return self.variables[step][variable_name]

    def SetSolutionStepValue(self, variable_name, step, value):
        if variable_name in list(self.variables[step].keys()):
            self.variables[step][variable_name] = value
        else:
            raise Exception(
                "trying to set an inexisting variable with name ",
                variable_name,
                " on node ",
                self.Id)
    def __str__(self):
        return  "Node #{0} with {1}".format(self.Id, self.variables)
        return s

class Dof:

    def __init__(self, node, variable):
        self.node = node
        self.variable = variable

    def IsFixed(self):
        return self.node.IsFixed(self.variable)

    def GetEquationId(self):
        return self.node.equation_id[self.variable]

    def SetEquationId(self, new_id):
        self.node.equation_id[self.variable] = new_id

    def GetValue(self, step=0):
        return self.node.GetSolutionStepValue(self.variable, step)

    def SetValue(self, value, step=0):
        self.node.SetSolutionStepValue(self.variable, step, value)

    def __eq__(self, other):
        if isinstance(other, Dof):
            return ((self.node.Id == other.node.Id) and (self.variable == other.variable))
        else:
            return False

    def __lt__(self, other):
        if(self.node.Id < other.node.Id):
            return True
        elif(self.node.Id == other.node.Id):
            return self.variable < other.variable
        else:
            return False

    def __hash__(self):
        return hash((self.node.Id, self.variable))
