from __future__ import print_function, absolute_import, division 
import math
from pyKratos import *
from numpy import *
from scipy import linalg


def Create(Id, prop, list_of_nodes):
    geom = triangle.Triangle(list_of_nodes)
    return StokesElement(Id, prop, geom)


class StokesElement:
    #this elements constructs a stiffness matrix which mixes velocities and pressures
    #for each pair of nodes I and J this matrix has a 3*3 subblock ordered as
    #  | Kvv Kvp | 
    #  | Kpv Kpp |3x3
    
    integration_order = 2  # this is like a c "static variable" one for all of the objects of this type
    include_dynamics = True

    def __init__(self, Id, prop, geometry):
        self.Id = Id
        self.prop = prop
        self.geometry = geometry

    def GetDofsPerNode(self):
        return 3

    def GetVectorValueOnGauss(self, var_x, var_y, N,step=0):
        value = zeros(2)
        for i in range(0, 3):
            value[0] += N[i] * self.geometry[i].GetSolutionStepValue(var_x, step)
            value[1] += N[i] * self.geometry[i].GetSolutionStepValue(var_y, step)
        return value
    
    def Assemble_vv_part(self,Kvv, LHS):
        nnodes = self.geometry.GetNumberOfNodes()
        for i in range(0,nnodes):
            for j in range(0,nnodes):
                for k in range(0,2):
                    for l in range(0,2):
                        LHS[i*3+k,j*3+l] += Kvv[i*2+k,j*2+l]

        return LHS

    def Assemble_vp_part(self,Kvp, LHS):
        nnodes = self.geometry.GetNumberOfNodes()
        for i in range(0,nnodes):
            for j in range(0,nnodes):
                for k in range(0,2):
                    LHS[i*3+k,j*3+2] += Kvp[i*2+k,j]
        return LHS
    
    def Assemble_pv_part(self,Kpv, LHS):
        nnodes = self.geometry.GetNumberOfNodes()
        for i in range(0,nnodes):
            for j in range(0,nnodes):
                for k in range(0,2):
                    LHS[i*3+2,j*3+k] += Kpv[i,j*2+k]
        return LHS
    
    def Assemble_pp_part(self,Kpp, LHS):
        nnodes = self.geometry.GetNumberOfNodes()
        for i in range(0,nnodes):
            for j in range(0,nnodes):
                LHS[i*3+2,j*3+2] += Kpp[i,j]
        return LHS
    
    def AssembleResiduals(self,rv,rp,RHS):
        nnodes = self.geometry.GetNumberOfNodes()
        for i in range(0,nnodes):
            RHS[i*3    ] += rv[i*2]
            RHS[i*3 +1 ] += rv[i*2+1]
            RHS[i*3 +2 ] += rp[i]
        return RHS    
            
    def CalculateLocalSystem(self,ProcessInfo):
        order = self.integration_order
        nnodes = self.geometry.GetNumberOfNodes()
        dofs_per_node = self.GetDofsPerNode()
        mat_size = nnodes*dofs_per_node
        
        [Ns, derivatives, weights] = self.geometry.ShapeFunctions(order)

        RHS = zeros(mat_size)  # no external forces so far
        LHS = zeros((mat_size,mat_size))

        number_of_gauss = len(Ns)
        
        Kvv = zeros((nnodes*2,nnodes*2)) 
        Kvp = zeros((nnodes*2,nnodes)) 
        Kpv = zeros((nnodes,nnodes*2)) 
        Kpp = zeros((nnodes,nnodes)) 
        
        rv = zeros(nnodes*2)
        rp = zeros(nnodes)
        
        nu = self.prop[VISCOSITY]
        density = self.prop[DENSITY]
        body_force = zeros(2)
        body_force[0] = self.prop[BODY_FORCE_X]
        body_force[1] = self.prop[BODY_FORCE_Y]
        
        # estimate element lenght
        Area = 0.0
        for gauss in range(0, number_of_gauss):
            Area += weights[gauss]
        h = sqrt(2.0 * Area)
        
        
        for gauss in range(0, number_of_gauss):
            weight = weights[gauss]
            N = Ns[gauss]
            DN_DX = derivatives[gauss]
            C = self.ComputeConstitutiveMatrix()
            B = self.ComputeB(DN_DX)
            
            ############################################
            ## GALERKIN TERMS
            ############################################
            #viscous compute v-v contribution Kvp
            Kvv += weight*density*dot( B.transpose() , dot(C,B) )
            
            #compute Kvp term (grad_q, p)
            for i in range(0,nnodes):
                for j in range(0,nnodes):
                    for k in range(0,2): #loop on dimension
                        Kvp[i*2+k,j] += -weight*DN_DX[i,k]*N[j]
                        
            #compute Kpv term - will write the p equation divided by the density for better scaling
            for i in range(0,nnodes):
                for j in range(0,nnodes):
                    for k in range(0,2): #loop on dimension
                        Kpv[i,j*2+k] += -weight*N[i]*DN_DX[j,k] #changed sign! to have symmetry
                        
            #compute Kpp term - no contripution of pp from galerking terms
            
            #contribution to residual 
            for i in range(0,nnodes):
                for k in range(0,2):
                    rv[i*2+k] += weight*density*N[i]*body_force[k]
                    
            #rp term has no galerkin contribution
            
            ############################################
            ## STABILIZATION TERMS
            ############################################
            tau = 1.0/(density*nu/(h*h) )
            
            #contribution to Kvv
            fvector = zeros(nnodes*2)
            for i in range(0,nnodes):
                for k in range(0,2):
                    fvector[i*2+k] = body_force[k]
            Bf = weight*density*nu*tau*dot(B,fvector)
            rv += dot(B.transpose(), Bf)
            
            #contribution to Kpp
            Kpp += -weight*tau*dot(DN_DX,DN_DX.transpose()) #changed sign to have symmetry
            
            #contribution to rp
            for i in range(0,nnodes):
                for k in range(0,2):
                    rp[i] += -weight*tau*DN_DX[i,k]*body_force[k] #changed sign to have symmetry
                    
            
            #tau2 contribution - divdiv contribution to Kvv
            tau2 = density*nu
            div_vec = zeros(nnodes*2)
            for i in range(0,nnodes):
                for k in range(0,2):
                    div_vec[2*i+k] = DN_DX[i,k]
            Kvv += weight*tau2*outer(div_vec,div_vec)
            
            ##ADD DYNAMICS IF NEEDED - using BDF2
            if(self.include_dynamics == True):
                #dt = ProcessInfo[DELTA_TIME]
                #if(dt == 0):
                    #raise Exception("Dt can not be zero!!")
                coeffs =  ProcessInfo[BDF_COEFFICIENTS]
                c0 = coeffs[0]
                c1 = coeffs[1]
                c2 = coeffs[2]
                                
                #part of "acc" to the RHS
                ##0vec = self.GetValues(0) #current step
                v1gauss = self.GetVectorValueOnGauss(VELOCITY_X,VELOCITY_Y,N,1) #old step
                v2gauss = self.GetVectorValueOnGauss(VELOCITY_X,VELOCITY_Y,N,2) #two steps ago
                arhs = c1*v1gauss + c2*v2gauss 
                #print(arhs)
                for i in range(0,nnodes):
                    rv[i*2   ] -= (weight*density*N[i])*arhs[0]
                    rv[i*2 +1] -= (weight*density*N[i])*arhs[1]
      
                #part of "acc" to the LHS
                for i in range(0, 3):
                    for j in range(0, 3):
                        tmp = (weight*density*c0)*N[i]*N[j] 
                        Kvv[i*2  , j*2  ] += tmp
                        Kvv[i*2+1, j*2+1] += tmp
                        
                #STABILIZATION - contribution to rp
                for i in range(0,nnodes):
                    for k in range(0,2):
                        rp[i] += weight*tau*DN_DX[i,k]*arhs[k] #changed sign to have symmetry
                
                #STABILIZATION - contribution to Kpv
                for i in range(0,nnodes):
                    for j in range(0,nnodes):
                        for k in range(0,2):
                            Kpv[i,2*j+k] += weight*tau*DN_DX[i,k]*N[j]

                    
        LHS = self.Assemble_vv_part(Kvv,LHS)
        LHS = self.Assemble_vp_part(Kvp,LHS)
        LHS = self.Assemble_pv_part(Kpv,LHS)
        LHS = self.Assemble_pp_part(Kpp,LHS)
        RHS = self.AssembleResiduals(rv,rp,RHS)
            
        #print(LHS)
        #err
        # compute RESIDUAL for the RHS
        # since the problem is LINEAR this can be done as
        # RHS = fext - LHS*values
        # note that this is done out of the integration loop!
        values = self.GetValues()  # get values of unknown at the nodes
        RHS -= dot(LHS, values)

        return [LHS, RHS]

    def ComputeB(self, DN_DX):

        B = zeros((3, 6))

        for i in range(0, 3):
            index = i * 2
            B[0, index + 0] = DN_DX[i, 0]
            B[0, index + 1] = 0
            B[1, index + 0] = 0
            B[1, index + 1] = DN_DX[i, 1]
            B[2, index + 0] = DN_DX[i, 1]
            B[2, index + 1] = DN_DX[i, 0]

        return B

    def ComputeConstitutiveMatrix(self):
        C = zeros((3, 3))
        nu = self.prop[VISCOSITY]
        density = self.prop[DENSITY]

        C[0,0] = 2.0*nu*density #TODO: check this!!!
        C[1,1] = 2.0*nu*density #TODO: check this!!!
        C[2,2] = nu*density #TODO: check this!!!
        return C

    # this function returns a list with the node and unkowns to be solved for
    def GetDofList(self):
        unknowns = []
        unknowns.append(Dof(self.geometry[0], VELOCITY_X))
        unknowns.append(Dof(self.geometry[0], VELOCITY_Y))
        unknowns.append(Dof(self.geometry[0], PRESSURE))
        unknowns.append(Dof(self.geometry[1], VELOCITY_X))
        unknowns.append(Dof(self.geometry[1], VELOCITY_Y))
        unknowns.append(Dof(self.geometry[1], PRESSURE))
        unknowns.append(Dof(self.geometry[2], VELOCITY_X))
        unknowns.append(Dof(self.geometry[2], VELOCITY_Y))
        unknowns.append(Dof(self.geometry[2], PRESSURE))
        return unknowns

    def EquationId(self):
        equation_ids = []
        equation_ids.append(self.geometry[0].EquationId(VELOCITY_X))
        equation_ids.append(self.geometry[0].EquationId(VELOCITY_Y))
        equation_ids.append(self.geometry[0].EquationId(PRESSURE))
        equation_ids.append(self.geometry[1].EquationId(VELOCITY_X))
        equation_ids.append(self.geometry[1].EquationId(VELOCITY_Y))
        equation_ids.append(self.geometry[1].EquationId(PRESSURE))
        equation_ids.append(self.geometry[2].EquationId(VELOCITY_X))
        equation_ids.append(self.geometry[2].EquationId(VELOCITY_Y))
        equation_ids.append(self.geometry[2].EquationId(PRESSURE))
        return equation_ids

    def GetValues(self, step=0):
        values = zeros(3*self.geometry.GetNumberOfNodes())
        values[0] = self.geometry[0].GetSolutionStepValue(VELOCITY_X, step)
        values[1] = self.geometry[0].GetSolutionStepValue(VELOCITY_Y, step)
        values[2] = self.geometry[0].GetSolutionStepValue(PRESSURE, step)
        values[3] = self.geometry[1].GetSolutionStepValue(VELOCITY_X, step)
        values[4] = self.geometry[1].GetSolutionStepValue(VELOCITY_Y, step)
        values[5] = self.geometry[1].GetSolutionStepValue(PRESSURE, step)
        values[6] = self.geometry[2].GetSolutionStepValue(VELOCITY_X, step)
        values[7] = self.geometry[2].GetSolutionStepValue(VELOCITY_Y, step)
        values[8] = self.geometry[2].GetSolutionStepValue(PRESSURE, step)
        return values


