from __future__ import print_function, absolute_import, division 
import math
from pyKratos import *
from numpy import *
from scipy import linalg


def Create(Id, prop, list_of_nodes):
    geom = triangle.Triangle(list_of_nodes)
    return ConvDiffElement2D(Id, prop, geom)


class ConvDiffElement2D:
    integration_order = 2  # this is like a c "static variable" one for all of the objects of this type

    def __init__(self, Id, prop, geometry):
        self.Id = Id
        self.prop = prop
        self.geometry = geometry

    def GetDofsPerNode(self):
        return 1


    def GetGaussPointVelocity(self, N,step=0):
        value = zeros(2)
        for i in range(0, 3):
            value[0] += N[i] * self.geometry[i].GetSolutionStepValue(VELOCITY_X, step)
            value[1] += N[i] * self.geometry[i].GetSolutionStepValue(VELOCITY_Y, step)
        return value

    def GetGaussPointValue(self, N, variable,step=0):
        value = 0.0
        for i in range(0, 3):
            value += N[i] * self.geometry[i].GetSolutionStepValue(variable, step)
        return value

    def CalculateLocalSystem(self,ProcessInfo):
        #print("line41")
        order = self.integration_order
        [Ns, derivatives, weights] = self.geometry.ShapeFunctions(order)
        number_of_gauss = len(Ns)

        RHS = zeros(3)  # no external forces so far
        LHS = zeros((3, 3))

        conductivity = self.prop[CONDUCTIVITY]
        
        dt = ProcessInfo[DELTA_TIME]
        if(dt == 0):
            raise Exception("Dt can not be zero!!")


        # estimate element lenght
        Area = 0.0
        for gauss in range(0, number_of_gauss):
            Area += weights[gauss]
        h = sqrt(2.0 * Area)

        for gauss in range(0, number_of_gauss):
            N = Ns[gauss]
            DN_DX = derivatives[gauss]
            

            vgauss = self.GetGaussPointVelocity(N)

            # compute contribution of convection matrix
            conv_tmp = zeros((3, 3))
            for i in range(0, 3):
                for j in range(0, 3):
                    vgauss_dot_DNj = 0.0
                    for k in range(0, 2):  # sum onto the components
                        vgauss_dot_DNj += vgauss[k] * DN_DX[j, k]
                    conv_tmp[i, j] += N[i] * vgauss_dot_DNj
            
            # compute diffusion contribution matrix
            diff_tmp = conductivity * dot(DN_DX, DN_DX.transpose())

            # compute convection stabilization term
            tau_denom = linalg.norm(vgauss)/h + conductivity/(h*h)
            tau_denom = max(tau_denom, 0.001/dt) #this is a trick... 0.001 is arbitrary            
            Tau_SUPG = 1.0/tau_denom
            
            v_outer_v = outer(vgauss, vgauss)
            tmp1 = dot(DN_DX, v_outer_v)
            stab_tmp = Tau_SUPG *  dot(tmp1, DN_DX.transpose())

            # add up contributions to form the LHS
            LHS += weights[gauss] * (conv_tmp + diff_tmp + stab_tmp)

            # source term - rhs contribution
            # note that the part of the stabilization that is in the LHS will
            # be added at the end
            source_gauss = self.GetGaussPointValue(N, VOLUME_HEAT_SOURCE)
            source_tmp = source_gauss * N

            # compute stabilization contribution related to the source term
            vDN = zeros(3)
            for i in range(0, 3):
                for k in range(0, 2):
                    vDN[i] += DN_DX[i, k] * vgauss[k]
            stab_rhs_tmp = (Tau_SUPG * source_gauss) * vDN

            # add source and related contribution to the RHS
            RHS += weights[gauss] * (source_tmp + stab_rhs_tmp)
            
            #########################################################
            ## add dynamic contribution - use BDF2
            c0 = 1.5/dt #3/(2dt)
            c1 = -2.0/dt
            c2 = 0.5/dt
                            
            #part of "acc" to the RHS
            ##0vec = self.GetValues(0) #current step
            T1gauss = self.GetGaussPointValue(N,TEMPERATURE,1) #old step
            T2gauss = self.GetGaussPointValue(N,TEMPERATURE,2) #two steps ago
            arhs = c1*T1gauss + c2*T2gauss
            
            
            RHS -= (weights[gauss]*arhs ) *N
            
            
            #part of "acc" to the LHS
            M = zeros((3,3))
            for i in range(0, 3):
                for j in range(0, 3):
                    M[i, j] += N[i] * N[j] 
            LHS += (weights[gauss]*c0)*M
            
            #part of stabilization in the LHS
            RHS -= (weights[gauss]*Tau_SUPG) * arhs * vDN #part of stabilization in the uRHS
            LHS += (weights[gauss]*Tau_SUPG*c0)*outer(vDN,N)
            

                    
        # compute RESIDUAL for the RHS
        # since the problem is linear this can be done as
        # RHS = fext - LHS*values
        # note that this is done out of the integration loop!
        values = self.GetValues()
        RHS -= dot(LHS, values)
        #print("line150")
        return [LHS, RHS]

    # this function returns a list with the node and unkowns to be solved for
    def GetDofList(self):
        unknowns = []
        unknowns.append(Dof(self.geometry[0], TEMPERATURE))
        unknowns.append(Dof(self.geometry[1], TEMPERATURE))
        unknowns.append(Dof(self.geometry[2], TEMPERATURE))
        return unknowns

    def EquationId(self):
        equation_ids = []
        equation_ids.append(self.geometry[0].EquationId(TEMPERATURE))
        equation_ids.append(self.geometry[1].EquationId(TEMPERATURE))
        equation_ids.append(self.geometry[2].EquationId(TEMPERATURE))
        return equation_ids

    def GetValues(self, step=0):
        values = zeros(3)
        values[0] = self.geometry[0].GetSolutionStepValue(TEMPERATURE, step)
        values[1] = self.geometry[1].GetSolutionStepValue(TEMPERATURE, step)
        values[2] = self.geometry[2].GetSolutionStepValue(TEMPERATURE, step)
        return values

    def GetFirstTimeDerivatives(self, step=0):
        raise Exception("not implemented")

    def GetSecondTimeDerivatives(self, step=0):
        raise Exception("not implemented")
