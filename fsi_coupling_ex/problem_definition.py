from __future__ import print_function, absolute_import, division 
import sys
sys.path.append("..")
print(sys.path)

from numpy import *
from pyKratos import *



import domain_definition1
import domain_definition2
solver1 = domain_definition1.domain1_solver #mp with dirichlet conditions on interface
solver2 = domain_definition2.domain2_solver

##we apply dirichlet condition to solver1 fsi nodes
for node in solver1.model_part.NodeIterators():
    if(node.GetSolutionStepValue(FSI_INTERFACE,0) == 1):
        node.Fix(VELOCITY_X) #note that we fix VELOCITY and note displacements since that is the unknown    
        node.Fix(VELOCITY_Y)

import mapper_utility
mapper = mapper_utility.Mapper(solver1.model_part, solver2.model_part)




#preparing for output
zero_based_indices_for_nodes = True
GiDIO1 = gid_io.GidIO("inputfile.mdpa","domain1",zero_based_indices_for_nodes)
GiDIO1.WriteMesh(solver1.model_part,"outmesh1")
GiDIO2 = gid_io.GidIO("inputfile.mdpa","domain2",zero_based_indices_for_nodes)
GiDIO2.WriteMesh(solver2.model_part,"outmesh2")

def CoupledSolve(solver1, solver2, mapper, w, max_iter=10, tol=1e-6):


    # mapper.Map_2to1(DISPLACEMENT_X, DISPLACEMENT_X)
    # mapper.Map_2to1(DISPLACEMENT_Y, DISPLACEMENT_Y)
    mapper.Map_2to1(VELOCITY_X, VELOCITY_X)
    mapper.Map_2to1(VELOCITY_Y, VELOCITY_Y)

    norm_rd = 1e9
    it = 0
    while norm_rd > tol and it < max_iter:

        tmp = []
        for node in mapper.interface1:
            tmp.append(node.GetSolutionStepValue(VELOCITY_X,0))
            tmp.append(node.GetSolutionStepValue(VELOCITY_Y,0))
        dbefore = array(tmp)

        solver1.strategy.Solve()
        solver1.ComputeReactions()  

        mapper.Map_1to2(REACTION_X, EXTERNAL_FORCE_X, -1.0) #note that we change sign to the reaction
        mapper.Map_1to2(REACTION_Y, EXTERNAL_FORCE_Y, -1.0)

        solver2.strategy.Solve()

        # mapper.Map_2to1(DISPLACEMENT_X, DISPLACEMENT_X)
        # mapper.Map_2to1(DISPLACEMENT_Y, DISPLACEMENT_Y)
        mapper.Map_2to1(VELOCITY_X, VELOCITY_X)
        mapper.Map_2to1(VELOCITY_Y, VELOCITY_Y)

        tmp = []
        for node in mapper.interface1:
            tmp.append(node.GetSolutionStepValue(VELOCITY_X,0))
            tmp.append(node.GetSolutionStepValue(VELOCITY_Y,0))
        dafter = array(tmp)
        # print("debefore = ", dbefore)
        # print("dafter = ", dafter)

        rd = dafter - dbefore
        norm_rd = linalg.norm(rd)
        print("iteration = ", it, " norm_rd = ", norm_rd)

        drelaxed = dbefore + w*rd
        counter = 0
        for node in mapper.interface1:
            node.SetSolutionStepValue(VELOCITY_X,0,drelaxed[counter])
            node.SetSolutionStepValue(VELOCITY_Y,0,drelaxed[counter+1])
            counter += 2

        it += 1
        
        
relaxation_factor = 0.1    
max_iter = 300
tol = 1e-6
dt = 0.01
nsteps = 100
for i in range(1,nsteps):
    time = i*dt
    solver1.model_part.CloneTimeStep(time)
    solver2.model_part.CloneTimeStep(time)

    if(i > 2):
        print("time: {:.6}".format(time))
        CoupledSolve(solver1,solver2,mapper, relaxation_factor, max_iter, tol)
        
        GiDIO1.WriteNodalResults(DISPLACEMENT,solver1.model_part.NodeIterators(), time)
        GiDIO1.WriteNodalResults(VELOCITY,solver1.model_part.NodeIterators(), time)
        GiDIO1.WriteNodalResults(REACTION,solver1.model_part.NodeIterators(), time)
        GiDIO1.WriteNodalResults(FSI_INTERFACE,solver1.model_part.NodeIterators(), time)
        GiDIO1.WriteNodalResults(EXTERNAL_FORCE,solver1.model_part.NodeIterators(), time)
            
        GiDIO2.WriteNodalResults(DISPLACEMENT,solver2.model_part.NodeIterators(), time)
        GiDIO2.WriteNodalResults(VELOCITY,solver2.model_part.NodeIterators(), time)
        GiDIO2.WriteNodalResults(REACTION,solver2.model_part.NodeIterators(), time)
        GiDIO2.WriteNodalResults(FSI_INTERFACE,solver2.model_part.NodeIterators(), time)
        GiDIO2.WriteNodalResults(EXTERNAL_FORCE,solver2.model_part.NodeIterators(), time)
        
    #print(d)

    #import matplotlib.pyplot as plt
    #plt.plot(t,d)
    #plt.show()