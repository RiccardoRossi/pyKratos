from __future__ import print_function, absolute_import, division 
from numpy import *
from .variables import *


class SolvingStrategy:

    def __init__(self, model_part, scheme, builder_and_solver):
        self.model_part = model_part
        self.scheme = scheme
        self.builder_and_solver = builder_and_solver
        self.A = zeros((0, 0))
        self.b = zeros((0))
        self.dx = zeros((0))

    def Initialize(self):
        # find list of dofs
        self.builder_and_solver.SetupDofSet()

        # allocate memory for the system
        self.A, self.x, self.b = self.builder_and_solver.SetupSystem(
            self.A, self.dx, self.b)

    def Solve(self):
        # eventually do prediction
        self.scheme.Predict()

        # do build and solve
        self.A, self.dx, self.b = self.builder_and_solver.BuildAndSolve(
            self.A, self.x, self.b)

        # call scheme to do update
        self.scheme.Update(self.builder_and_solver.dofset, self.dx)

    def SpyMatrix(self):
        try:
            import matplotlib.pylab as pl
            pl.spy(self.A)
            pl.show()
        except:
            raise Exception(
                "error in function Spy. Probably matplotlib not installed")
