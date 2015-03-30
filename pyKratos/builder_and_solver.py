from __future__ import print_function, absolute_import, division 
from numpy import *
from scipy import linalg
from scipy import sparse


class BuilderAndSolver:
    use_sparse_matrices = True

    '''ATTENTION!!
    this builder and solver assumes elements to be written IN RESIDUAL FORM and hence
    solves FOR A CORRECTION Dx'''

    def __init__(self, model_part, scheme):
        self.scheme = scheme
        self.model_part = model_part
        self.dofset = set()
        self.dirichlet_dof = set()

    def SetupDofSet(self):
        '''this function shapes the system to be built'''

        # start by iterating over all the elements and obtaining the list of
        # dofs
        aux = set()
        for elem in self.model_part.ElementIterators():
            unknowns = elem.GetDofList()

            for aaa in unknowns:
                aux.add(aaa)

        self.dofset = sorted(aux)

        # for dof in self.dofset:
            # print dof.node.Id, " ",dof.variable

        # assign an equation id
        counter = 0
        for dof in self.dofset:
            dof.SetEquationId(counter)
            counter += 1

            if(dof.IsFixed()):
                self.dirichlet_dof.add(dof)

    def SetupSystem(self, A, dx, b):
        ndof = len(self.dofset)

        # allocate systme vectors
        b = zeros(ndof)
        dx = zeros(ndof)

        # allocate system matrix
        if(self.use_sparse_matrices == False):  # dense case
            A = zeros((ndof, ndof))
        else:  # allocate non zeros and transofrm to csr
            A = sparse.dok_matrix((ndof, ndof))
            for elem in self.model_part.ElementIterators():
                # get non zero positions
                equation_id = self.scheme.EquationId(elem)
                for i in range(0, len(equation_id)):
                    eq_i = equation_id[i]
                    for j in range(0, len(equation_id)):
                        eq_j = equation_id[j]
                        A[eq_i,
                            eq_j] = 1.0  # set it to 1 to ensure it is converted well
                        # problem here is that in converting zero entries are
                        # discarded
            A = A.tocsr()

        return [A, dx, b]

    # this function sets to
    def SetToZero(self, A, dx, b):
        ndof = len(self.dofset)

        if(self.use_sparse_matrices == False):
            # allocating a dense matrix. This should be definitely improved
            A = zeros((ndof, ndof))
            b = zeros(ndof)
            dx = zeros(ndof)
        else:
            # print A.todense()
            A = A.multiply(
                0.0)  # only way i found to set to zero is to multiply by zero
            b = zeros(ndof)
            dx = zeros(ndof)
        return [A, dx, b]

    def ApplyDirichlet(self, A, dx, b):
        ndof = A.shape[0]
        if(self.use_sparse_matrices == False):
            for dof in self.dirichlet_dof:
                fixed_eqn = dof.GetEquationId()
                for i in range(0, ndof):
                    A[fixed_eqn, i] = 0.0
                    A[i, fixed_eqn] = 0.0
                A[fixed_eqn, fixed_eqn] = 1.0
                b[fixed_eqn] = 0.0  # note that this is zero since we assume residual form!
        else:
            # expensive loop: exactly set to 1 the diagonal
            # could be done cheaper, but i want to guarantee accuracy
            aux = ones(ndof, dtype=bool)
            for dof in self.dirichlet_dof:
                eq_id = dof.GetEquationId()
                aux[eq_id] = False


            ij = A.nonzero()
            for i, j in zip(ij[0], ij[1]):
                if(aux[i] == False or aux[j] == False):
                    A[i, j] = 0.0

            for dof in self.dirichlet_dof:
                eq_id = dof.GetEquationId()
                A[eq_id, eq_id] = 1.0
                b[eq_id] = 0.0

        return [A, dx, b]
    
    def Build(self, A, dx, b):
        A, dx, b = self.SetToZero(A, dx, b)

        for elem in self.model_part.ElementIterators():
            # find where to assemble
            equation_id = self.scheme.EquationId(elem)

            # compute local contribution to the stiffness matrix
            [lhs, rhs] = self.scheme.CalculateLocalSystem(elem)

            # assembly to the matrix
            for i in range(0, len(equation_id)):
                eq_i = equation_id[i]
                b[eq_i] += rhs[i]
                for j in range(0, len(equation_id)):
                    eq_j = equation_id[j]
                    A[eq_i, eq_j] += lhs[i, j]

        for cond in self.model_part.ConditionIterators():
            # find where to assemble
            equation_id = self.scheme.EquationId(cond)

            # compute local contribution to the stiffness matrix
            [lhs, rhs] = self.scheme.CalculateLocalSystem(cond)

            # assembly to the matrix
            for i in range(0, len(equation_id)):
                eq_i = equation_id[i]
                b[eq_i] += rhs[i]
                for j in range(0, len(equation_id)):
                    eq_j = equation_id[j]
                    A[eq_i, eq_j] += lhs[i, j]
        return [A, dx, b]

    def BuildAndSolve(self, A, dx, b):
        A, dx, b = self.Build(A, dx, b)

        A, dx, b = self.ApplyDirichlet(A, dx, b)

        # print A
        if(self.use_sparse_matrices == False):
            dx = linalg.solve(A, b)
        else:

            from scipy.sparse.linalg import spsolve
            dx = sparse.linalg.spsolve(A, b)

        return [A, dx, b]
