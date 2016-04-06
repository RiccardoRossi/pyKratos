# ----------------------------------------------------------------------
# author  : Martin Ruchti
# contact : martin.ruchti@tum.de
# ----------------------------------------------------------------------
from __future__ import print_function, absolute_import, division 
import math
from pyKratos import *
from numpy import *
from scipy import linalg

'''
    Navier-Stokes element - Variational Multiscale (VMS)
    
    This element contains the formulation of the Navier Stokes VMS element
    with the bossak time scheme included. 
    
    ***NOTE*** this element is formulated in the residual form,
               solving the element system, gives the correction steps
               dui, dvi, dpi
    
    Stabilization terms come from the ASGS formulation
    
    ***NOTE***
    use this element with:
    bossak_scheme.py
    
'''


def Create(Id, prop, list_of_nodes):
    geom = triangle.Triangle(list_of_nodes)
    return NavierStokesElement(Id, prop, geom)


class NavierStokesElement:
    #this elements constructs a stiffness matrix which mixes velocities and pressures
    #for each pair of nodes I and J this matrix has a 3*3 subblock ordered as
    #  | Kvv Kvp | 
    #  | Kpv Kpp |3x3
    
    integration_order = 1  # this is like a c "static variable" one for all of the objects of this type
    include_dynamics = False

    def __init__(self, Id, prop, geometry):
        self.Id         = Id
        self.prop       = prop
        self.geometry   = geometry
        self.nnodes     = self.geometry.GetNumberOfNodes()
        dofs_per_node   = self.GetDofsPerNode()
        self.mat_size   = self.nnodes*dofs_per_node
        #self.tauDyn     = 0.01# from kratos
        self.tauDyn     = 0.05# channel example
        #self.tauDyn     = 0.8 # results in a strong incompressibility constraint
        self.integration_order = 1

    def GetDofsPerNode(self):
        return 3

    def GetVectorValueOnGauss(self, var_x, var_y, N,step=0):
        value = zeros(2)
        for i in range(0, 3):
            value[0] += N[i] * self.geometry[i].GetSolutionStepValue(var_x, step)
            value[1] += N[i] * self.geometry[i].GetSolutionStepValue(var_y, step)
        return value
    
    def Assemble_vv_part(self,Kvv, LHS):
        for i in range(0,self.nnodes):
            for j in range(0,self.nnodes):
                for k in range(0,2):
                    for l in range(0,2):
                        LHS[i*3+k,j*3+l] += Kvv[i*2+k,j*2+l]

        return LHS

    def Assemble_vp_part(self,Kvp, LHS):
        for i in range(0,self.nnodes):
            for j in range(0,self.nnodes):
                for k in range(0,2):
                    LHS[i*3+k,j*3+2] += Kvp[i*2+k,j]
        return LHS
    
    def Assemble_pv_part(self,Kpv, LHS):
        for i in range(0,self.nnodes):
            for j in range(0,self.nnodes):
                for k in range(0,2):
                    LHS[i*3+2,j*3+k] += Kpv[i,j*2+k]
        return LHS
    
    def Assemble_pp_part(self,Kpp, LHS):
        for i in range(0,self.nnodes):
            for j in range(0,self.nnodes):
                LHS[i*3+2,j*3+2] += Kpp[i,j]
        return LHS
    
    def AssembleResiduals(self,rv,rp,RHS):
        for i in range(0,self.nnodes):
            RHS[i*3    ] += rv[i*2]
            RHS[i*3 +1 ] += rv[i*2+1]
            RHS[i*3 +2 ] += rp[i]
        return RHS    
            
            
    def CalculateLocalSystem(self,ProcessInfo, stepRelative = 0):
        
        order           = self.integration_order
        [Ns, derivatives, weights] = self.geometry.ShapeFunctions(order)
        number_of_gauss = len(Ns)
        dt              = ProcessInfo[DELTA_TIME]
        coeffs          = ProcessInfo[BOSSAK_COEFFICIENTS]
        alpha           = coeffs[0]
        gamma           = coeffs[1]
        
        #create empty matrices/vectors of the right size
        RHS = zeros(self.mat_size)
        LHS = zeros((self.mat_size,self.mat_size))
        #mass stabilization matrix for rhs
        Mp  = zeros((self.nnodes, self.nnodes * 2))

        #create empty matrices for the different contributions of the LHS matrix
        Kvv = zeros((self.nnodes*2,self.nnodes*2))
        Kvp = zeros((self.nnodes*2,self.nnodes)) 
        Kpv = zeros((self.nnodes,self.nnodes*2))
        Kpp = zeros((self.nnodes,self.nnodes)) 

        #create empty vectors for the different contributions of the RHS vector
        rv = zeros(self.nnodes*2)
        rp = zeros(self.nnodes)
        
        #compute advection velocity and taus for every gauss point
        adVel = self.computeAdvectionVelocity(Ns, derivatives, weights, stepRelative)
        [tau1, tau2] = self.computeTauOneAndTwo(number_of_gauss, adVel, dt)
        
        lumpedM = self.computeLumpedMassMatrix(Ns, weights)
        #lumpedM = self.computeConsistentMassMatrixContribution(Ns, weights)
        
        #mass stabilization term
        self.computeMassStabilizationTermVelocity(lumpedM, Ns, derivatives, weights, tau1, number_of_gauss, adVel)
        #Kvv
        #inertia term**
        self.computeLHSInertia(Kvv, lumpedM, alpha, gamma, dt)
        #convection term
        self.computeConvectionTerm(Kvv, Ns, derivatives, weights, adVel)
        #convection stab. term
        self.computeConvectionStabilization(Kvv, adVel, tau1, Ns, derivatives, weights)
        #velocity stabilization
        self.computeVelocityStabilization(Kvv, derivatives, weights, number_of_gauss, tau2)
        #viscous term
        self.computeViscousTerm(Kvv, derivatives, weights, number_of_gauss)
        #Kvp
        #pressure term pv
        #pressure stab. term pv
        self.computeKvp(Kvp, adVel, Ns, derivatives, weights, tau1, number_of_gauss)
        #Kpv
        #pressure term qu
        #pressure stab. term qu
        self.computeKpv(Kpv, adVel, Ns, derivatives, weights, tau1, number_of_gauss)
        #Kpp
        #pressure stab. term
        self.computeKpp(Kpp, derivatives, weights, tau1, number_of_gauss)
            
        #assemble LHS
        LHS     = self.Assemble_vv_part(Kvv,LHS)
        LHS     = self.Assemble_vp_part(Kvp,LHS)
        LHS     = self.Assemble_pv_part(Kpv,LHS)
        LHS     = self.Assemble_pp_part(Kpp,LHS)

        #RHS
        #compute RHS mass stabilization terms for the pressure block
        self.computeMassStabilizationTermPressure(Mp , Ns, derivatives, weights, tau1, number_of_gauss, adVel)
        LHS     = self.Assemble_pv_part(((1.0-alpha)/(gamma*dt))*Mp,LHS)
        #intertial term**
        self.computeRHSInertia(rv, rp, lumpedM, Mp, alpha, gamma, dt)
        #convection stab. term vf
        self.computeConvectionStabilizationRHS(rv, derivatives, weights, number_of_gauss, adVel, tau1)
        #external forces
        self.computeExternalForce(rv, Ns, weights)
        #projection terms ***
        #pressure stab. term qf
        self.computePressureStabilizationRHS(rp, derivatives, weights, number_of_gauss, tau1)
        
        
        
        #assemble RHS
        RHS     = self.AssembleResiduals(rv, rp, RHS)
        
        # compute RESIDUAL for the RHS
        # the problem is NON-LINEAR however, this is done as
        # RHS = fext - LHS*values
        # For a solution, the local system has to be re-calculated iteratively,
        # to approach the equilibrium point for this time step
        
        values = self.GetValues(stepRelative)  # get values of unknown at the nodes
        RHS -= dot(LHS, values)

        #print("Shape RHS = ", RHS.shape)
        #print("Shape LHS = ", LHS.shape)
        #input("Printed shapes")

        #import matplotlib.pylab as pl
        #pl.spy(LHS)
        #pl.show()
        #input("Spyed matrix")

        return [LHS, RHS]
        
        
    def computeConvectionStabilization(self, Kvv, adVel, tauOne, Ns, derivatives, weights):

        density = self.prop[DENSITY]
        
        for gauss in range(0, len(Ns)):

            weight = weights[gauss]
            N      = Ns[gauss]
            DN_DX  = derivatives[gauss]
            AGradN = zeros(self.nnodes)
            AdvVel = adVel[gauss]
            tau1   = tauOne[gauss]

            for i in range(0,self.nnodes):
                AGradN[i] = AdvVel[0] * DN_DX[i,0] + AdvVel[1] * DN_DX[i,1]

            for i in range(0,self.nnodes):
                for j in range(0,self.nnodes):
                    Kvv[i*2,j*2]            += weight * density**2 * tau1 * AGradN[i] * AGradN[j]
                    Kvv[i*2+1,j*2+1]        += weight * density**2 * tau1 * AGradN[i] * AGradN[j]


    def computeConvectionTerm(self, Kvv, Ns, derivatives, weights, adVels):

        density  = self.prop[DENSITY]
        
        for gauss in range(0, len(Ns)):
            weight      = weights[gauss]
            N           = Ns[gauss]
            DN_DX       = derivatives[gauss]
            adVel       = adVels[gauss]

            #compute scalar product a*grad(N)
            AGradN = zeros(self.nnodes)
            for i in range(0,self.nnodes):
                AGradN[i] = adVel[0] * DN_DX[i,0] + adVel[1] * DN_DX[i,1]
                
            for i in range(0,self.nnodes):
                for j in range(0,self.nnodes):
                    Kvv[i*2,j*2]            += weight * density * N[i] * AGradN[j]
                    Kvv[i*2+1,j*2+1]        += weight * density * N[i] * AGradN[j]

    def computeVelocityStabilization(self, Kvv, derivatives, weights, number_of_gauss, tauTwo):

        for gauss in range(0, number_of_gauss):
            tau2  = tauTwo[gauss]
            DN_DX = derivatives[gauss]
            weight = weights[gauss]
            
            for i in range(0,self.nnodes):
                for j in range(0,self.nnodes):
                    Kvv[i*2,j*2]     += tau2 * DN_DX[i,0] * DN_DX[j,0] * weight
                    Kvv[i*2,j*2+1]   += tau2 * DN_DX[i,0] * DN_DX[j,1] * weight
                    Kvv[i*2+1,j*2]   += tau2 * DN_DX[i,1] * DN_DX[j,0] * weight
                    Kvv[i*2+1,j*2+1] += tau2 * DN_DX[i,1] * DN_DX[j,1] * weight

    def computeViscousTerm(self, Kvv, derivatives, weights, number_of_gauss):
    
        nu          = self.prop[VISCOSITY]
        density     = self.prop[DENSITY]
        mu          = nu*density
        fourThirds  = 4.0/3.0
        twoThirds   = 2.0/3.0
        
        for gauss in range(0, number_of_gauss):
            DN_DX = derivatives[gauss]
            weight = weights[gauss]
            
            for i in range(0,self.nnodes):
                for j in range(0,self.nnodes):
                
                    #diagonal terms
                    Kvv[i*2,j*2]     += mu * weight * ( fourThirds * DN_DX[i,0] * DN_DX[j,0] +  DN_DX[i,1] * DN_DX[j,1] )
                    Kvv[i*2+1,j*2+1] += mu * weight * ( DN_DX[i,0] * DN_DX[j,0] +  fourThirds * DN_DX[i,1] * DN_DX[j,1] )
                    #off-diagonal terms
                    Kvv[i*2,j*2+1]   += mu * weight * ( -twoThirds * DN_DX[i,0] * DN_DX[j,1] +  DN_DX[i,1] * DN_DX[j,0] )
                    Kvv[i*2+1,j*2]   += mu * weight * ( -twoThirds * DN_DX[i,1] * DN_DX[j,0] +  DN_DX[i,0] * DN_DX[j,1] )

    
    def computeLHSInertia(self, Kvv, lumpedM, alpha, gamma, dt):
        
        coeff   = ( 1.0 - alpha ) / ( gamma * dt )
        
        Kvv += coeff * lumpedM
     
    
    def computeMassStabilizationTermVelocity(self, lumpedM, Ns, derivatives, weights, tau1, number_of_gauss, adVel):
        
        density     = self.prop[DENSITY]
        
        for gauss in range(0, number_of_gauss):
            N       = Ns[gauss]
            DN_DX   = derivatives[gauss]
            weight  = weights[gauss]
            tauOne  = tau1[gauss]
            AdvVel  = adVel[gauss]
            
            coeff   = weight * tauOne * density**2
            
            #compute scalar product a*grad(N)
            AGradN = zeros(self.nnodes)
            for i in range(0,self.nnodes):
                AGradN[i] = AdvVel[0] * DN_DX[i,0] + AdvVel[1] * DN_DX[i,1]
            
            for i in range(0,self.nnodes):
                for j in range(0,self.nnodes):
                    
                    lumpedM[2*i, 2*j]       += coeff * AGradN[i] * N[j]
                    lumpedM[2*i+1, 2*j+1]   += coeff * AGradN[i] * N[j]
                    
        
    def computeKvp(self, Kvp, adVel, Ns, derivatives, weights, tauOne, number_of_gauss):
    
        density     = self.prop[DENSITY]
    
        for gauss in range(0,number_of_gauss):
            N       = Ns[gauss]
            DN_DX   = derivatives[gauss]
            weight  = weights[gauss]
            tau1    = tauOne[gauss]
            AdvVel  = adVel[gauss]
            
            #compute scalar product a*grad(N)
            AGradN = zeros(self.nnodes)
            for i in range(0,self.nnodes):
                AGradN[i] = AdvVel[0] * DN_DX[i,0] + AdvVel[1] * DN_DX[i,1]
                
            for i in range(0,self.nnodes):
                for j in range(0,self.nnodes):
                    Kvp[i*2,j]    += weight * ( - DN_DX[i,0] * N[j] + density * tau1 * AGradN[i] * DN_DX[j,0] )   
                    Kvp[i*2+1,j]  += weight * ( - DN_DX[i,1] * N[j] + density * tau1 * AGradN[i] * DN_DX[j,1] )
            
            
    def computeKpv(self, Kpv, adVel, Ns, derivatives, weights, tau1, number_of_gauss):
    
        density     = self.prop[DENSITY]
        
        for gauss in range(0,number_of_gauss):
            N       = Ns[gauss]
            DN_DX   = derivatives[gauss]
            weight  = weights[gauss]
            tauOne  = tau1[gauss]
            AdvVel = adVel[gauss]
            
            #compute scalar product a*grad(N)
            AGradN = zeros(self.nnodes)
            for i in range(0,self.nnodes):
                AGradN[i] = AdvVel[0] * DN_DX[i,0] + AdvVel[1] * DN_DX[i,1]
                
            for i in range(0,self.nnodes):
                for j in range(0,self.nnodes):
                    Kpv[i,j*2]    += weight * ( + N[i] * DN_DX[j,0] + density * tauOne * DN_DX[i,0] * AGradN[j] )
                    Kpv[i,j*2+1]  += weight * ( + N[i] * DN_DX[j,1] + density * tauOne * DN_DX[i,1] * AGradN[j] )
                
                
    def computeKpp(self, Kpp, derivatives, weights, tauOne, number_of_gauss):
            
        for gauss in range(0,number_of_gauss):
            DN_DX   = derivatives[gauss]
            weight  = weights[gauss]
            tau1    = tauOne[gauss]
            
            for i in range(0,self.nnodes):
                for j in range(0,self.nnodes):
                    Kpp[i,j] += weight * tau1 * ( DN_DX[i,0] * DN_DX[j,0] + DN_DX[i,1] * DN_DX[j,1] )
                    
                    
    def computeConvectionStabilizationRHS(self, rv, derivatives, weights, number_of_gauss, adVels, tauOne):
        
        nu              = self.prop[VISCOSITY]
        density         = self.prop[DENSITY]
        f               = zeros(2)
        f[0]            = self.prop[BODY_FORCE_X]
        f[1]            = self.prop[BODY_FORCE_Y]
        
        for gauss in range(0,number_of_gauss):
            DN_DX       = derivatives[gauss]
            weight      = weights[gauss]
            tau1        = tauOne[gauss]
            advVel      = adVels[gauss]
            
            #compute scalar product a*grad(N)
            AGradN = zeros(self.nnodes)
            for i in range(0,self.nnodes):
                AGradN[i]       = advVel[0] * DN_DX[i,0] + advVel[1] * DN_DX[i,1]
                
            for i in range(0,self.nnodes):
                rv[i*2]         += tau1 * weight * density * f[0] * density * AGradN[i]
                rv[i*2+1]       += tau1 * weight * density * f[1] * density * AGradN[i]
                
               
    def computeExternalForce(self, rv, Ns, weights):
        
        density = self.prop[DENSITY]
        f       = zeros(2)
        f[0]    = self.prop[BODY_FORCE_X]
        f[1]    = self.prop[BODY_FORCE_Y]
        
        for gauss in range(0,len(Ns)):
            weight      = weights[gauss]
            N           = Ns[gauss]
            
            for i in range(0,self.nnodes):
                rv[i*2]         += weight * density * N[i] * f[0]
                rv[i*2+1]       += weight * density * N[i] * f[1]
        
    def computePressureStabilizationRHS(self, rp, derivatives, weights, number_of_gauss, tauOne):
        
        density = self.prop[DENSITY]
        f       = zeros(2)
        f[0]    = self.prop[BODY_FORCE_X]
        f[1]    = self.prop[BODY_FORCE_Y]
        
        for gauss in range(0,number_of_gauss):
            DN_DX       = derivatives[gauss]
            weight      = weights[gauss]
            tau1        = tauOne[gauss]
            
            for i in range(0,self.nnodes):
                rp[i]   += weight * tau1 * ( DN_DX[i, 0] * density * f[0] + DN_DX[i, 1] * f[1] )
        
        
    def computeMassStabilizationTermPressure(self, Mp , Ns, derivatives, weights, tau1, number_of_gauss, adVels):
            
        density     = self.prop[DENSITY]
        dim         = 2
        
        for gauss in range(0, number_of_gauss):
            N       = Ns[gauss]
            DN_DX   = derivatives[gauss]
            weight  = weights[gauss]
            tauOne  = tau1[gauss]
            AdvVel   = adVels[gauss]
            
            coeff   = weight * tauOne * density
            
            #compute scalar product a*grad(N)
            AGradN = zeros(self.nnodes)
            for i in range(0,self.nnodes):
                AGradN[i] = AdvVel[0] * DN_DX[i,0] + AdvVel[1] * DN_DX[i,1]
                
            for i in range(0,self.nnodes):
                for j in range(0,self.nnodes):
                    Mp[i,2*j]   += coeff * DN_DX[i,0] * N[j]
                    Mp[i,2*j+1] += coeff * DN_DX[i,1] * N[j] 
                    
     
    def computeRHSInertia(self, rv, rp, lumpedM, Mp, alpha, gamma, dt):
        
        coeff1  = ( 1.0 - alpha ) / ( gamma * dt )
        coeff2  = ( 1.0 - gamma - alpha ) / ( gamma )
        
        u       = zeros(2*self.nnodes)
        uDot    = zeros(2*self.nnodes)
        for k in range(0,self.nnodes):
            
            #explicite values for RHS inertia term
            u[k*2]        = self.geometry[k].GetSolutionStepValue(VELOCITY_X, 1)
            u[k*2+1]      = self.geometry[k].GetSolutionStepValue(VELOCITY_Y, 1)
            
            uDot[k*2]     = self.geometry[k].GetSolutionStepValue(ACCELERATION_X, 1)
            uDot[k*2+1]   = self.geometry[k].GetSolutionStepValue(ACCELERATION_Y, 1)

        rv  += dot( lumpedM, ( coeff1 * u + coeff2 * uDot ) )
        rp  += dot( Mp, ( coeff1 * u + coeff2 * uDot )  )
        
        
    def computeTauOneAndTwo(self, number_of_gauss, advVels, dt):

        nu              = self.prop[VISCOSITY]
        density         = self.prop[DENSITY]
        
        x0 = self.geometry[0].coordinates[0]
        x1 = self.geometry[1].coordinates[0] 
        x2 = self.geometry[2].coordinates[0] 
        y0 = self.geometry[0].coordinates[1]
        y1 = self.geometry[1].coordinates[1] 
        y2 = self.geometry[2].coordinates[1] 
    
        detJ = (x1-x0) * (y2-y0) - (y1-y0) * (x2-x0)
        Area = 0.5*detJ
        
        #create empty container
        tau1 = []
        tau2 = []
        
        #ElemSize = sqrt(2.0*Area)
        ElemSize = 1.128379167 * sqrt(Area)

        for gauss in range(0,number_of_gauss):

            #compute norm of advection velocity for gauss point
            AdvVelNorm = 0.0
            for k in range(0,2):
                AdvVelNorm += advVels[gauss][k] * advVels[gauss][k]
                
            AdvVelNorm = sqrt(AdvVelNorm)

            tauOne = 1.0 / (density * (self.tauDyn / dt + 2.0 * AdvVelNorm / ElemSize + 4.0 * nu / (ElemSize**2)))
            tauTwo = density*(nu + 0.5*ElemSize*AdvVelNorm)

            tau1.append(tauOne)
            tau2.append(tauTwo)
            
        return [tau1, tau2]
    

    def computeAdvectionVelocity(self, Ns, derivatives, weights, stepRelative):

        #create empty container
        AdvVel = []

        number_of_gauss = len(Ns)
        for gauss in range(0, number_of_gauss):
            weight = weights[gauss]
            N      = Ns[gauss]
            DN_DX  = derivatives[gauss]

            AdvVelGauss    = zeros(2)
            AdvVelGauss[0] = N[0] * self.geometry[0].GetSolutionStepValue(VELOCITY_X,stepRelative) \
                           + N[1] * self.geometry[1].GetSolutionStepValue(VELOCITY_X,stepRelative) \
                           + N[2] * self.geometry[2].GetSolutionStepValue(VELOCITY_X,stepRelative)
            AdvVelGauss[1] = N[0] * self.geometry[0].GetSolutionStepValue(VELOCITY_Y,stepRelative) \
                           + N[1] * self.geometry[1].GetSolutionStepValue(VELOCITY_Y,stepRelative) \
                           + N[2] * self.geometry[2].GetSolutionStepValue(VELOCITY_Y,stepRelative)

            AdvVel.append(AdvVelGauss)
        
        return AdvVel
    
    def computeVelocityNorm(self, adVels, ngauss):
        
        velNorm = []
        
        for gauss in range(0,ngauss):
            norm = 0.0
            adVel   = adVels[gauss]
            #consider 2 dimensions here
            for i in range(0,2):
                norm += adVel[i]**2
            norm = sqrt(norm)
            velNorm.append(norm)
            
        return velNorm
    
    
    def computeLumpedMassMatrix(self, Ns, weights):
    
        M       = zeros((2*self.nnodes,2*self.nnodes))
        density = self.prop[DENSITY]
        
        for gauss in range(0,len(Ns)):
            N           = Ns[gauss]
            weight      = weights[gauss]
            for i in range(0,self.nnodes):
                M[i*2, i*2]         += density * weight * N[i]
                M[i*2+1, i*2+1]     += density * weight * N[i]

        return M
    
    
    #direct reimplementation from kratos vms element
    def computeConsistentMassMatrixContribution(self, Ns, weights):
        
        density     = self.prop[DENSITY]
        dimension   = 2
        M           = zeros((2*self.nnodes,2*self.nnodes))
        
        for gauss in range(0,len(Ns)):
            
            N       = Ns[gauss]
            weight  = weights[gauss]

            coef = density * weight;

            # Note: Dof order is (vx,vy,[vz,]p) for each node
            for i in range(0,self.nnodes):
                for j in range(0,self.nnodes):
                    temp = coef * N[i] * N[j]
    
                    M[i*2,j*2]      += temp
                    M[i*2+1,j*2+1]  += temp
                    
        return M
    
    
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

    def GetValues(self, step):
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
