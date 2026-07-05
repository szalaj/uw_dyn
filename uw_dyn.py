# -*- coding: utf-8 -*-
# Program do obliczen dynamiki 3D 
# metoda ukladow wieloczlonowych

# stworzone w ramach pracy magisterskiej, 2016

#autor: Marcin Szalajski



#uzyte moduly (biblioteki)
from scipy.integrate import odeint,ode
import scipy.linalg
import scipy.optimize
import numpy as np
from math import cos,sin,pi
import matplotlib.pyplot as plt
import csv
from numpy import unravel_index
import functools
import copy
import sys



#funkcja tworząca wektor r czlonu i-tego
def r_i(i,q):
    r = q[3*(i-1):3*(i-1)+3]
    rr = wektor(r[0],r[1],r[2])
    return rr

#funkcja tworząca wektor dr czlonu i-tego
def dr_i(i,dq):
    dr = dq[3*(i-1):3*(i-1)+3]
    drr = wektor(dr[0],dr[1],dr[2])
    return drr    


#funkcja tworząca wektor p czlonu i-tego
def p_i(i,q,N):
    p = q[3*N+4*(i-1):3*N+4*(i-1)+4]
    #pp = np.array([ [p[0]], [p[1]], [p[2]], [p[3]] ])
    pp = wektor_p(p[0],p[1],p[2],p[3])
    return pp


#funkcja tworząca wektor dp czlonu i-tego
def dp_i(i,dq,N):
    dp = dq[3*N+4*(i-1):3*N+4*(i-1)+4]
    dpp = wektor_p(dp[0],dp[1],dp[2],dp[3])
    return dpp

#funkcja obliczajaca parametry eulera
#u - wektor tworzaca os obrotu
#chi - kat obrotu
def u2p(u,chi):
    e0=cos(chi/2)
    e=u*sin(chi/2)
    p=np.array([e0,e[0],e[1],e[2]])
    return p
    
#funkcja obliczaja parametry eulera
#z katow eulera    
def EA_to_EP(fi1,fi2,fi3):
    e0 = cos(0.5*(fi1+fi3))*cos(0.5*fi2)
    e1 = cos(0.5*(fi1-fi3))*sin(0.5*fi2)
    e2 = sin(0.5*(fi1-fi3))*sin(0.5*fi2)
    e3 = sin(0.5*(fi1+fi3))*cos(0.5*fi2)
    return e0,e1,e2,e3

# macierz obrotowa
def R(p):
    e0=p[0][0]
    e1=p[1][0]
    e2=p[2][0]
    e3=p[3][0]
    Rot =2*np.array([[e0*e0+e1*e1-0.5, e1*e2-e0*e3, e1*e3+e0*e2],
                     [e1*e2+e0*e3, e0*e0+e2*e2-0.5, e2*e3-e0*e1],
                     [e1*e3-e0*e2, e2*e3+e0*e1, e0*e0+e3*e3-0.5]])
    return Rot

# macierz pomocnicza
def G(p):
    #print(p)
    e0=p[0][0]
    e1=p[1][0]
    e2=p[2][0]
    e3=p[3][0]
    #print("ee",e0)
    Gie =np.array([[-e1, e0, e3, -e2],
                   [-e2, -e3, e0, e1],
                   [-e3, e2, -e1, e0]])
    #print(Gie)
    return Gie
    
# pochodna macierzy pomocniczej
def dG(dp):
    de0=dp[0][0]
    de1=dp[1][0]
    de2=dp[2][0]
    de3=dp[3][0]
    dGie =np.array([[-de1, de0, de3, -de2],
                    [-de2, -de3, de0, de1],
                    [-de3, de2, -de1, de0]])
    return dGie
    
# tworzenie z wektoru macierzy
# skosno-symetrycznej
def skew(a):
    ax=a[0]
    ay=a[1]
    az=a[2]
    A = np.array([[0, -az, ay],
    [az, 0, -ax],
    [-ay, ax, 0]])
    return A
  
# tworzenie wektoru pionowego r    
def wektor(ax,ay,az):
    return np.array([[ax],[ay],[az]])
    
# tworzenie wektoru pionowego p    
def wektor_p(e0,e1,e2,e3):
    return np.array([[e0],[e1],[e2],[e3]])
    
        
class Czlon:
    
    def __init__(self, numer, masa, tensor_bez):
        self.i = numer 
        self.m = masa
        self.J = tensor_bez
    
    # macierz masowa    
    def M(self):
        return self.m*np.eye(3)
    
    
    
# klasy par kinematycznych i kierujacych

    
class Polaczenie:
    
    def __init__(self):
        self.m = 0
        

        
class Para_Prostopadla(Polaczenie):
    
    def __init__(self, i, j, ai, aj):
        self.i = i
        self.j = j
        self.ai = ai
        self.aj = aj
        self.m = 1
        
    def wiezyK(self,q,N):
        pj = p_i(self.j,q,N)
        Rj = R(pj)
        F=0
        if self.i == 0:

            F = self.ai.transpose().dot(Rj).dot(self.aj)
            
        else:
            pi = p_i(self.i, q, N)     
            Ri = R(pi)       
            F = self.ai.transpose().dot(Ri.transpose()).dot(Rj).dot(self.aj)
            
        return F
        

        
    def jakobianK(self,q,N):
        Fqi=np.zeros([1,7])
        Fqj=np.zeros([1,7])
        pj = p_i(self.j,q,N)
        Rj = R(pj)
        Gj = G(pj)
        
        if self.i == 0:
            Fqj[0,3:7]= -2*self.ai.transpose().dot(Rj).dot(skew(self.aj)).dot(Gj)
        else: 
            pi=p_i(self.i,q,N)
            Ri=R(pi)
            Gi=G(pi)
            Fqi[0,3:7] = -2*self.aj.transpose().dot(Rj.transpose()).dot(Ri).dot(skew(self.ai)).dot(Gi)      
            Fqj[0,3:7] = -2*self.ai.transpose().dot(Ri.transpose()).dot(Rj).dot(skew(self.aj)).dot(Gj) 
            
        return Fqi, Fqj
        
    def gammaK(self,q,dq,N):
        
        pj = p_i(self.j,q,N)
        dpj = dp_i(self.j,dq,N)
        
        Rj = R(pj)
        Gj = G(pj)
        dGj = dG(dpj)
        
        gamK=0
        
        Om_j_skew = 2*Gj.dot(dGj.transpose())
        Om_j = 2*Gj.dot(pj)
          
        if self.i == 0:        
            pom = Om_j_skew.dot(Om_j_skew).dot(Rj.transpose())
            gamK = -self.aj.transpose().dot(pom).dot(self.ai)
        else:
            pi = p_i(self.i,q,N)
            dpi = dp_i(self.i,dq,N)
            Gi = G(pi)
            dGi = dG(dpi)
            
            Om_i_skew = 2*Gi.dot(dGi.transpose())
            Om_i = 2*Gi.dot(pi)
            
            Ri = R(pi)
            Rj = R(pj)

            aj_skew = skew(self.aj)    
            ai_skew = skew(self.ai)
            
            pom1 = Rj.transpose().dot(Ri).dot(Om_i_skew).dot(Om_i_skew)
            pom2 = Om_j_skew.dot(Om_j_skew).dot(Rj.transpose()).dot(Ri)
            pom3 = 2*Om_j.transpose().dot(aj_skew).dot(Rj.transpose()).dot(Ri).dot(ai_skew).dot(Om_i)

            gamK= -self.aj.transpose().dot(pom1+pom2).dot(self.ai) + pom3
            
        return gamK        
        

class Para_Prostopadla_D(Polaczenie):
    
    def __init__(self, i, j, sA_i, sB_j, ai):
        self.i = i
        self.j = j
        self.sA_i = sA_i
        self.sB_j = sB_j
        self.ai = ai
        self.m = 1
        
    def wiezyK(self,q,N):
        rj = r_i(self.j,q)
        pj = p_i(self.j,q,N)
        Rj = R(pj)
        F=0
        if self.i == 0:
            
            F = self.ai.transpose().dot(rj+Rj.dot(self.sB_j))-self.ai.transpose().dot(self.sA_i)
            
        else:
            ri = r_i(self.i,q)
            pi = p_i(self.i,q,N)     
            Ri = R(pi)       
            F = self.ai.transpose().dot(Ri.transpose()).dot(rj+Rj.dot(self.sB_j)-ri)-self.ai.transpose().dot(self.sA_i)
            
        return F

    def jakobianK(self,q,N):
        Fqi=np.zeros([1,7])
        Fqj=np.zeros([1,7])
        pj = p_i(self.j,q,N)
        Rj = R(pj)
        Gj = G(pj)
        
        if self.i == 0:
            Fqj[0,0:3]= self.ai.transpose() 
            Fqj[0,3:7]= -2*self.ai.transpose().dot(Rj).dot(skew(self.sB_j)).dot(Gj)
        else: 
            pi=p_i(self.i,q,N)
            Ri=R(pi)
            Gi=G(pi)
            
            Fqi[0,0:3]= -self.ai.transpose().dot(Ri.transpose()) 
            dij = rj+Rj.dot(self.sB_j)-ri-Ri.dot(self.sA_i)
            pom = self.ai.transpose().dot(skew(self.sA_i))-dij.transpose().dot(Ri).dot(skew(self.ai))
            Fqi[0,3:7] = 2*pom.dot(Gi)     
            
            Fqj[0,0:3]= self.ai.transpose().dot(Ri.transpose()) 
            Fqj[0,3:7] = -2*self.ai.transpose().dot(Ri.transpose()).dot(Rj).dot(skew(self.sB_j)).dot(Gj) 
            
        return Fqi, Fqj
        
        
    def gammaK(self,q,dq,N):
        
        pj = p_i(self.j,q,N)
        dpj = dp_i(self.j,dq,N)
        Rj = R(pj)
        Gj = G(pj)
        dGj = dG(dpj)
        
        gamK=0
        
        Om_j_skew = 2*Gj.dot(dGj.transpose())
        Om_j = 2*Gj.dot(pj)
        
        
        if self.i == 0:        
            
            pom = Om_j_skew.dot(Om_j_skew).dot(Rj.transpose())
            

            gamK = -self.sB_j.transpose().dot(pom).dot(self.ai)

        else:
            ri = r_i(self.i,q)
            rj = r_i(self.j,q)
            dri = dr_i(self.i,dq)
            drj = dr_i(self.j,dq)
            pi = p_i(self.i,q,N)
            dpi = dp_i(self.i,dq,N)
            
            Gi = G(pi)
            dGi = dG(dpi)
            
            Om_i_skew = 2*Gi.dot(dGi.transpose())
            Om_i = 2*Gi.dot(pi)
            
            Ri = R(pi)
            Rj = R(pj)
            
            uj_skew = skew(uj)    
            vi_skew = skew(vi)
            wi_skew = skew(wi)
            vj_skew = skew(vj)
            
       
            dij = rj + Rj.dot(self.sB_j) - ri - Ri.dot(self.sA_i)
            
            pom1 = 2*Om_i.transpose().dot(skew(self.ai)).dot(Ri.transpose()).dot(dri-drj)
            pom2 = 2*self.sB_j.transpose().dot(Om_j_skew).dot(Rj.transpose()).dot(Ri).dot(Om_i_skew).dot(self.ai)
            pom3 = self.sA_i.transpose().dot(Om_i_skew).dot(Om_i_skew).dot(self.ai)
            pom4 = self.sB_j.transpose().dot(Om_j_skew).dot(Om_j_skew).dot(Rj.transpose()).dot(Ri).dot(self.ai)
            pom5 = dij.transpose().dot(Ri).dot(Om_i_skew).dot(Om_i_skew).dot(self.ai)
            
            gamK = pom1 + pom2 - pom3 - pom4 - pom5
            
        return gamK
        

class Para_Sferyczna(Polaczenie):
    
    def __init__(self, i, j, sA_i, sB_j):
        self.i = i
        self.j = j
        self.sA_i = sA_i
        self.sB_j = sB_j
        self.m = 3    
        
    def wiezyK(self,q,N):
        rj = r_i(self.j,q)
        pj = p_i(self.j,q,N)
        Rj = R(pj)
        F=np.zeros([3,1])
        if self.i == 0:
            
            F[0:3,:] = rj + Rj.dot(self.sB_j) - self.sA_i
            
        else:
            ri = r_i(self.i,q)
            pi = p_i(self.i,q,N)     
            Ri = R(pi)       
            F[0:3,:] = rj + Rj.dot(self.sB_j) - ri - Ri.dot(self.sA_i)
            
        return F
        
    def jakobianK(self,q,N):
        pj = p_i(self.j,q,N)
        Rj = R(pj)
        Gj = G(pj)
        Fqi=np.zeros([3,7])
        Fqj=np.zeros([3,7])
        
        if self.i == 0:
            
            Fqj[0:3,0:3]=np.eye(3)
            Fqj[0:3,3:7]=-2*Rj.dot(skew(self.sB_j)).dot(Gj)
            
        else: 
            pi=p_i(self.i,q,N)
            Ri = R(pi)
            Gi = G(pi)
            Fqi[0:3,0:3]=-np.eye(3)
            Fqi[0:3,3:7]=2*Ri.dot(skew(self.sA_i) ).dot(Gi)
            
            Fqj[0:3,0:3]=np.eye(3)
            Fqj[0:3,3:7]=-2*Rj.dot(skew(self.sB_j)).dot(Gj)

        return Fqi, Fqj
    
    def gammaK(self,q,dq,N):
        
        pj = p_i(self.j,q,N)
        dpj = dp_i(self.j,dq,N)
        Rj = R(pj)
        Gj = G(pj)
        dGj = dG(dpj)
        
        gamK=np.zeros([3,1])
        
        Om_j_skew = 2*Gj.dot(dGj.transpose())
        Om_j = 2*Gj.dot(pj)
    
        if self.i == 0:        
            
            gamK[0:3] = -Rj.dot(Om_j_skew).dot(Om_j_skew).dot(self.sB_j)

        else:
            
            pi = p_i(self.i,q,N)
            dpi = dp_i(self.i,dq,N)
            
            Gi = G(pi)
            dGi = dG(dpi)
            
            Om_i_skew = 2*Gi.dot(dGi.transpose())
            Om_i = 2*Gi.dot(pi)
            
            Ri = R(pi)
            Rj = R(pj)
            

            pom1 = Ri.dot(Om_i_skew).dot(Om_i_skew).dot(self.sA_i)
            pom2 = Rj.dot(Om_j_skew).dot(Om_j_skew).dot(self.sB_j)
            
            gamK[0:3] = pom1 - pom2

        
        return gamK  
        
class Polaczenie_Obr(Para_Sferyczna):
    
    def __init__(self, i, j, sA_i, sB_j, vi, wi, uj):
        super().__init__(i,j,sA_i,sB_j)
        self.pp1 = Para_Prostopadla(i,j,vi,uj)
        self.pp2 = Para_Prostopadla(i,j,wi,uj)
        self.m = 5
        
    def wiezyK(self,q,N):

        F= super().wiezyK(q,N)
        
        F = np.append(F,self.pp1.wiezyK(q,N),0)
        F = np.append(F,self.pp2.wiezyK(q,N),0)
        
        return F
        
        
    def jakobianK(self,q,N):

        
        Fqi,Fqj = super().jakobianK(q,N)
        
        Fqi = np.append(Fqi, self.pp1.jakobianK(q,N)[0],0)
        Fqi = np.append(Fqi, self.pp2.jakobianK(q,N)[0],0)
        
        Fqj = np.append(Fqj, self.pp1.jakobianK(q,N)[1],0)
        Fqj = np.append(Fqj, self.pp2.jakobianK(q,N)[1],0)
        
        
        return Fqi, Fqj
        
    def gammaK(self,q,dq,N):
        
        gamK=super().gammaK(q,dq,N)
        
        gamK = np.append(gamK,self.pp1.gammaK(q,dq,N),0)
        gamK = np.append(gamK,self.pp2.gammaK(q,dq,N),0)

        return gamK         
     
  
class Polaczenie_Cyl(Polaczenie):
    
    def __init__(self, i, j, sA_i, sB_j, vi, wi, uj):
        self.i = i
        self.j = j
        self.p_prost1 = Para_Prostopadla(i,j,vi,uj)
        self.p_prost2 = Para_Prostopadla(i,j,wi,uj)
        self.p_prost_dij_1 = Para_Prostopadla_D(i,j,sA_i,sB_j,vi)
        self.p_prost_dij_2 = Para_Prostopadla_D(i,j,sA_i,sB_j,wi)
        self.m = 4
        
    def wiezyK(self,q,N):

        F= np.zeros([4,1])
                 
        F[0,:] = self.p_prost1.wiezyK(q,N)
        F[1,:] = self.p_prost2.wiezyK(q,N)
        F[2,:] = self.p_prost_dij_1.wiezyK(q,N)
        F[3,:] = self.p_prost_dij_2.wiezyK(q,N)

        return F
        
    def jakobianK(self,q,N):
        

        Fqi=np.zeros([4,7])
        Fqj=np.zeros([4,7])

        Fqi[0,:]= self.p_prost1.jakobianK(q,N)[0]
        Fqi[1,:]= self.p_prost2.jakobianK(q,N)[0]
        Fqi[2,:]= self.p_prost_dij_1.jakobianK(q,N)[0]
        Fqi[3,:]= self.p_prost_dij_2.jakobianK(q,N)[0]
        
        Fqj[0,:]= self.p_prost1.jakobianK(q,N)[1]
        Fqj[1,:]= self.p_prost2.jakobianK(q,N)[1]
        Fqj[2,:]= self.p_prost_dij_1.jakobianK(q,N)[1]
        Fqj[3,:]= self.p_prost_dij_2.jakobianK(q,N)[1]     
              
        return Fqi, Fqj
        
    def gammaK(self,q,dq,N):
    
        gamK=np.zeros([4,1])
        
        gamK[0]= self.p_prost1.gammaK(q,dq,N)
        gamK[1]=  self.p_prost2.gammaK(q,dq,N)
        gamK[2] = self.p_prost_dij_1.gammaK(q,dq,N)
        gamK[3] = self.p_prost_dij_2.gammaK(q,dq,N)
        
        return gamK
        
class Polaczenie_Przes(Polaczenie_Cyl):
    
    def __init__(self, i, j, sA_i, sB_j, vi, wi, vj, uj):
        
        super().__init__(i,j,sA_i,sB_j,vi,wi,uj)

        self.p_prost3 = Para_Prostopadla(i,j,vi,vj)
        self.m = 5
        
    def wiezyK(self,q,N):
        
        F=super().wiezyK(q,N)
        F = np.append(F, self.p_prost3.wiezyK(q,N),0)

        return F
        
    def jakobianK(self,q,N):


        Fqi,Fqj=super().jakobianK(q,N)
        Fqi = np.append(Fqi, self.p_prost3.jakobianK(q,N)[0], 0)
        Fqj = np.append(Fqj, self.p_prost3.jakobianK(q,N)[1], 0)   
              
        return Fqi, Fqj
        
    def gammaK(self,q,dq,N):
    
        gamK=super().gammaK(q,dq,N)
        gamK = np.append(gamK,  self.p_prost3.gammaK(q,dq,N), 0)
        
        return gamK
        
        
class Odleglosc(Para_Prostopadla_D):

    def __init__(self, i, j, sA_i, sB_j, ai, C):
        
        super().__init__(i, j, sA_i, sB_j, ai)

        self.C = C
        self.m = 1
        
    def wiezyD(self,q,N):
        F = super().wiezyK(q,N) - self.C
        
        return F

    def jakobianD(self,q,N):
        Fqi,Fqj = super().jakobianK(q,N) 
            
        return Fqi, Fqj

class Kat(Polaczenie):
    
    def __init__(self, i, j, vi, wi, uj,fi):
        self.i = i
        self.j = j
        self.vi = vi
        self.wi = wi
        self.uj = uj
        self.fi = fi
        self.p_prost1 = Para_Prostopadla(i,j,vi,uj)
        self.p_prost2 = Para_Prostopadla(i,j,wi,uj)
        self.m = 1     
        
    # wiezy kierujace pary
    def wiezyD(self,q,N):
 
        F=0
        if self.i == 0:
 
            if np.fabs(np.cos(self.fi)) <= np.sqrt(2)*0.5:
                F = self.p_prost1.wiezyK(q,N) - np.cos(self.fi)
            else:
                F = self.p_prost2.wiezyK(q,N) - np.cos(self.fi + np.pi/2)

            
        else:       
            if np.fabs(np.cos(self.fi)) <= np.sqrt(2)*0.5:
                F = self.p_prost1.wiezyK(q,N) - np.cos(self.fi)
            else:
                F = self.p_prost2.wiezyK(q,N) - np.cos(self.fi + np.pi/2)
        return F
    
    # jakobian wiezow kierujacych pary    
    def jakobianD(self,q,N):

        if self.i == 0:
            if np.fabs(np.cos(self.fi)) <= np.sqrt(2)*0.5:
                Fqj= self.p_prost1.jakobianK(q,N)[1]
            else:
                Fqj= self.p_prost2.jakobianK(q,N)[1]
        else: 

            if np.fabs(np.cos(self.fi)) <= np.sqrt(2)*0.5:
                Fqi = self.p_prost1.jakobianK(q,N)[0]    
                Fqj = self.p_prost1.jakobianK(q,N)[1]
            else:
                Fqi = self.p_prost2.jakobianK(q,N)[0]     
                Fqj = self.p_prost2.jakobianK(q,N)[1] 
        return Fqi, Fqj        
        

# klasa sily wewnetrznej w ukladzie
class SilaWewnProst:
    def __init__(self, i, j, sA_i, sB_j, k,l0, c, F):
            self.i = i
            self.j = j
            self.sA_i = sA_i
            self.sB_j = sB_j
            self.k = k
            self.l0 = l0
            self.c = c
            self.F = F
            
    def sila(self,q,dq,N):
        rj = r_i(self.j,q)
        drj = dr_i(self.j,dq)
        
        pj = p_i(self.j,q,N)
        Rj = R(pj)
        Gj = G(pj)

        om_j = 2*Gj.dot(pj)

        Qr_i = np.zeros((3,1))
        Qp_i = np.zeros((4,1))         
        Qr_j = np.zeros((3,1))
        Qp_j = np.zeros((4,1))  
              
        if self.i==0:
            
            dij = rj + Rj.dot(self.sB_j) - self.sA_i
            l = np.sqrt( dij.transpose().dot(dij) )

            if l>0.01:
                pom = drj - Rj.dot(skew(self.sB_j)).dot(om_j)

                dl = (dij/l).transpose().dot(pom)
                f = self.k*(l-self.l0) + self.c*dl + self.F

                Qr_j = -(f/l)*dij
                Qp_j = -(f/l)*2*Gj.transpose().dot(skew(self.sB_j)).dot(Rj.transpose()).dot(dij)
             
        else:
            ri = r_i(self.i,q)
            dri = dr_i(self.i,dq)
            
            pi = p_i(self.i,q,N)
            Ri = R(pi)
            Gi = G(pi)
    

            om_i = 2*Gi.dot(pi)            
            
            dij = rj + Rj.dot(self.sB_j) - ri - Ri.dot(self.sA_i)
            l = np.sqrt( dij.transpose().dot(dij) )

            if l>0.01:
                pom = drj - Rj.dot(skew(self.sB_j)).dot(om_j) - dri + Ri.dot(skew(self.sA_i)).dot(om_i)

                dl = (dij/l).transpose().dot(pom)
                f = self.k*(l-self.l0) + self.c*dl + self.F
                #print('i: ',self.i,' j: ', self.j, ' f: ',f)

                Qr_i = (f/l)*dij
                Qp_i = (f/l)*2*Gi.transpose().dot(skew(self.sA_i)).dot(Ri.transpose()).dot(dij)
                
                Qr_j = -(f/l)*dij
                Qp_j = -(f/l)*2*Gj.transpose().dot(skew(self.sB_j)).dot(Rj.transpose()).dot(dij)
            
        return Qr_i, Qp_i, Qr_j, Qp_j
        
        
class SilaZewn:
    def __init__(self,czlon, rodzaj, wielkosc):
        self.czlon = czlon
        self.rodzaj = rodzaj
        self.wielkosc = wielkosc


def jakobian_p_kolumny(N):
    for i in range(0,N):
        yield i*7+3
        yield i*7+4        
        yield i*7+5
        yield i*7+6
        
def jakobian_r_kolumny(N):
    for i in range(0,N):
        yield i*7
        yield i*7+1       
        yield i*7+2


class Uklad:

    def __init__(self):
        self.czlony = []
        self.wiezy_k = [] # wiezy kinematyczne
        self.wiezy_d = [] # wiezy kierujace
        self.silyWewn = []
        self.silyZewn = []
        self.N = 0 #ilosc czlonow
        self.M = 0 #ilosc wiezow kinematycznych
        self.Mi = 0 #ilosc wiezow kierujacych (war.poczatkowe)
        self.Y = [] #wyniki symulacji
        self.grawitacja = True
        
        
    def dodajCzlon(self, czlon):
        #sprawdzenie czy 'czlon' jest typu Czlon
        if isinstance(czlon, Czlon):
            self.czlony.append(czlon)
            self.N += 1
        else:
            raise Exception("Obiekt nie jest typu Czlon")
            
    def dodajWiez(self, wiez):
        if isinstance(wiez, Polaczenie):
            self.wiezy_k.append(wiez)
            self.M += wiez.m
        else:
            raise Exception('Obiekt nie jest typu Polaczenie')
            
    def dodajWiezD(self, wiez):
        if isinstance(wiez, Polaczenie):
            self.wiezy_d.append(wiez)
            self.Mi += wiez.m
        else:
            raise Exception('Obiekt nie jest typu Polaczenie')        
            
    def dodajSileWewn(self, sila):
        self.silyWewn.append(sila)
        
    def dodajSileZewn(self, sila):
        self.silyZewn.append(sila)
            

    # wiezy kinematyczne zbiorczo
    def wiezyK(self, q):
    
        #
        F= np.zeros([self.M,1])
        k=0
        for w in self.wiezy_k:
            F[k:k+w.m,:] = w.wiezyK(q,self.N)
            k=k+w.m    

        return F
    
    # wiezy kierujace zbiorczo    
    def wiezyD(self, q):
    
        #
        F= np.zeros([self.Mi,1])
        k=0
        for w in self.wiezy_d:
            F[k:k+w.m,:] = w.wiezyD(q,self.N)
            k=k+w.m    

        return F
        
    
    # wiezy parametrow eulera zbiorczo
    def wiezyP(self, q):
        
        #liczba czlonow w ukladzie
        N = self.N
        
        Fp = np.zeros([N,1])
        
        for k in range(0,N):
            pi = p_i(k+1,q,N)
            Fp[k,:] = pi.transpose().dot(pi) - 1
    
        return Fp
        
    # wiezy kinematyczne i parametrow Eulera zbiorczo   
    def wiezyKP(self, q):
        
        wiezy_K = self.wiezyK(q)
        wiezy_P = self.wiezyP(q)

        F = np.concatenate((wiezy_K, wiezy_P), axis=0)
        return F
    
    # wiezy kinematyczne, kierujace i parametrow Eulera zbiorczo       
    def wiezyKPD(self, q):
        
        wiezy_K = self.wiezyK(q)
        wiezy_P = self.wiezyP(q)
        wiezy_D = self.wiezyD(q)
        
        F = np.concatenate((wiezy_K, wiezy_P, wiezy_D), axis=0)
        return F
        
    # jakobian wiezow kinematycznych zbiorczo
    def jakobianK(self, q):
        
        M=self.M
        N=self.N
        
        Fq=np.zeros([M,7*N])
        
        k=0
        for w in self.wiezy_k:
            if w.i == 0:
                Fq[k:k+w.m, 7*(w.j-1):7*(w.j-1)+7]=w.jakobianK(q,N)[1]
            else:
                Fq[k:k+w.m, 7*(w.i-1):7*(w.i-1)+7]=w.jakobianK(q,N)[0]
                Fq[k:k+w.m, 7*(w.j-1):7*(w.j-1)+7]=w.jakobianK(q,N)[1]
            k+=w.m
        
    
        return Fq
     
    # jakobian wiezow kierujacych zbiorczo
    def jakobianD(self, q):
        
        Mi=self.Mi
        N=self.N
        
        Fq=np.zeros([Mi,7*N])
        
        k=0
        for w in self.wiezy_d:
            if w.i == 0:
                Fq[k:k+w.m, 7*(w.j-1):7*(w.j-1)+7]=w.jakobianD(q,N)[1]
            else:
                Fq[k:k+w.m, 7*(w.i-1):7*(w.i-1)+7]=w.jakobianD(q,N)[0]
                Fq[k:k+w.m, 7*(w.j-1):7*(w.j-1)+7]=w.jakobianD(q,N)[1]
            k+=w.m
        
    
        return Fq
        
    # jakobian wiezow parametrow eulera zbiorczo
    def jakobianP(self, q):
        N=self.N
        Fq_p=np.zeros([N,7*N])
        for k in range(0,N):
            pi=q[3*N+4*k:3*N+4*k+4]
            #pi=p_i(k,q)
            Fq_p[k,k*7+3:k*7+7]=2*pi
            
        return Fq_p
    
    # jakobian wiezow kinematycznych i parametrow eulera zbiorczo
    def jakobianKP(self, q):
        
        jak_K = self.jakobianK(q)
        jak_P = self.jakobianP(q)
        
        jak = np.concatenate((jak_K, jak_P), axis=0)
    
        return jak
        
    # jakobian wiezow kin., kier. i parametrow eulera zbiorczo    
    def jakobianKPD(self, q):
        
        jak_K = self.jakobianK(q)
        jak_P = self.jakobianP(q)
        jak_D = self.jakobianD(q)
        
        jak = np.concatenate((jak_K, jak_P, jak_D), axis=0)
    
        return jak
       
    # wektor gamma wiezow kin. zbiorczo
    def gammaK(self, q,dq):
        
        M = self.M
    
        gamK= np.zeros([M,1])
        k=0
        for w in self.wiezy_k:
            gamK[k:k+w.m,:] = w.gammaK(q,dq,self.N)
            k=k+w.m
            
        return gamK
        
    # wektor gamma wiezow par. eulera zbiorczo
    def gammaP(self, dq):
        
        N=self.N
        
        gammP = np.zeros([N,1])
        
        for k in range(0,N):
            dpi=dq[3*N+4*k:3*N+4*k+4]
            gammP[k,:]=-2*dpi.transpose().dot(dpi)
    
        return gammP
    
    # macierz masowa zbiorczo   
    def zbiorczeM(self):
        
        czM =[]
        for cz in self.czlony:
            czM.append(cz.M())
            
        zbM = functools.reduce(scipy.linalg.block_diag, czM)
        
        return zbM
        
    # tensor bezwladnosci zbiorczo  
    def zbiorczeJ(self):
        
        czJ =[]
        for cz in self.czlony:
            czJ.append(cz.J)
            
        zbJ = functools.reduce(scipy.linalg.block_diag, czJ)
    
        return zbJ
    
    # macierz G zbiorczo 
    def zbiorczeG(self,q):
        
        N = self.N
        
        czG = []
        for k in range(0,N):
            pi=q[3*N+4*k:3*N+4*k+4]
            pi = wektor_p(pi[0],pi[1],pi[2],pi[3])
            #pi = p_i(k,q)
            
            czG.append(G(pi))
            
        zbG = functools.reduce(scipy.linalg.block_diag, czG)
    
        return zbG
        
    # macierz dG zbiorczo     
    def zbiorcze_dG(self,dq):
        N = self.N
        czG = []
        for k in range(0,N):
            dpi=dq[3*N+4*k:3*N+4*k+4]
            dpi = wektor_p(dpi[0],dpi[1],dpi[2],dpi[3])
            #dpi=dp_i(k,dq)
            czG.append(dG(dpi))
            
        zbdG = functools.reduce(scipy.linalg.block_diag, czG)
    
        return zbdG
        

    # jakobianK polozen zbiorczo
    def zbiorczeF_r(self,q):
        
        kol = list(jakobian_p_kolumny(self.N))
            
        F_q = self.jakobianK(q)
        #usuniecie kolumn odpowiadajacych z p
        F_r = scipy.delete(F_q, kol, 1)
        #F_r = F_q[:,0:3]
        return F_r
    
    # jakobianK par. eulera zbiorczo
    def zbiorczeF_p(self,q):
        

        kol = list(jakobian_r_kolumny(self.N))
            
        F_q = self.jakobianK(q)
        #usuniecie kolumn odpowiadajacych z r
        F_p = scipy.delete(F_q, kol, 1)
        return F_p
        
    # jakobianP par. eulera zbiorczo
    def zbiorczeFp_p(self,q):
        
        kol = list(jakobian_r_kolumny(self.N))
            
        Fp = self.jakobianP(q)
        #usuniecie kolumn odpowiadajacych z r
        Fp_p = scipy.delete(Fp, kol, 1)
        return Fp_p
        
    # parametry eulera zbiorczo    
    def zbiorcze_p(self,q):
        N=self.N
        p=[]
        for k in range(0,N):
            pi=q[3*N+4*k:3*N+4*k+4]
            #pi = p_i(k,dq)
            p.append(pi[0])
            p.append(pi[1])
            p.append(pi[2])
            p.append(pi[3])
            
        p = np.array(p)
        #zamiana z array 1D na 2D 
        return p.reshape((-1,1))
    
    # lewa strona rownania dynamiki
    def Lstrona(self,q,dq):
        N = self.N
        M = self.M
        #M - masa , może sie mylic z liczba wiezow kinematycznych M
        zbM = self.zbiorczeM()
        J = self.zbiorczeJ()
        G = self.zbiorczeG(q)
        F_r = self.zbiorczeF_r(q)
        F_p = self.zbiorczeF_p(q)
        Fp_p = self.zbiorczeFp_p(q)
        
        GT = G.transpose()
        Iloczyn = 4*GT.dot(J).dot(G)
        
        # wiersze duzej macierzy lewej strony 
        Lstr1 = np.concatenate((zbM, np.zeros([3*N,4*N]), F_r.transpose(), np.zeros([3*N,N])), axis=1)
        Lstr2 = np.concatenate((np.zeros([4*N,3*N]), Iloczyn, F_p.transpose(), Fp_p.transpose() ), axis=1)
        Lstr3 = np.concatenate((F_r, F_p, np.zeros([M,M]), np.zeros([M,N])), axis=1)
        Lstr4 = np.concatenate((np.zeros([N,3*N]), Fp_p, np.zeros([N,M]), np.zeros([N,N]) ), axis=1)
        
        # zlozenie duzej macierzy lewje strony
        Lstr = np.concatenate((Lstr1, Lstr2, Lstr3, Lstr4), axis=0)
    
        return Lstr
        
        
    # prawa strona rownania dynamiki
    def Pstrona(self,q,dq):
        N = self.N
        
        G = self.zbiorczeG(q)
        dG = self.zbiorcze_dG(dq)
        dGT = dG.transpose()
        
        J = self.zbiorczeJ()
        
        
        FA = np.zeros([3*N,1])
        nA = np.zeros([3*N,1]) 

        # dodanie sil zewnetrznych jesli istnieja
        if len(self.silyZewn):
            for s in self.silyZewn:
                i = s.czlon-1
                if s.rodzaj == 'Fx':
                    FA[3*i,:] = s.wielkosc
                elif s.rodzaj == 'Fy':
                    FA[3*i+1,:] = s.wielkosc
                elif s.rodzaj == 'Fz':
                    FA[3*i+2,:] = s.wielkosc
                elif s.rodzaj == 'nx':
                    nA[3*i,:] = s.wielkosc
                elif s.rodzaj == 'ny':
                    nA[3*i+1,:] = s.wielkosc
                elif s.rodzaj == 'nz':
                    nA[3*i+2,:] = s.wielkosc
                else:
                    raise Exception('zla sila')
                    
        #dodanie sily grawitacyjnej            
        if self.grawitacja:
            for cz in self.czlony:
                FA[3*(cz.i-1)+2,:] += -9.80665*cz.m
        
        p = self.zbiorcze_p(q)
                                        
        Ilo = 2*G.transpose().dot(nA) + 8*dGT.dot(J).dot(dG).dot(p)
        
        #jesli sa dodane jakies sily wewn
        if len(self.silyWewn):
            for s in self.silyWewn:
                if s.i == 0:
                    Qr_i, Qp_i, Qr_j, Qp_j = s.sila(q,dq,self.N)
                    #print(Qr_j)
                    FA[3*(s.j-1):3*(s.j-1)+3] += Qr_j
                    Ilo[4*(s.j-1):4*(s.j-1)+4] += Qp_j
                else:
                   Qr_i, Qp_i, Qr_j, Qp_j = s.sila(q,dq,self.N)
                   FA[3*(s.i-1):3*(s.i-1)+3] += Qr_i
                   Ilo[4*(s.i-1):4*(s.i-1)+4] += Qp_i
                   FA[3*(s.j-1):3*(s.j-1)+3] += Qr_j
                   Ilo[4*(s.j-1):4*(s.j-1)+4] += Qp_j
                   

        
        # zamiana kolejnosci wspolrzednych dopasowanych do jakobianu
        # q-> r1 r2 rn p1 p2 pn. bedzie r1 p1 r2 p2 rn pn
        
        dq_jak = []
        
        for k in range(0,N):
            dri=dq[3*k:3*k+3]
            
            dq_jak.append(dri[0])
            dq_jak.append(dri[1])
            dq_jak.append(dri[2])
            
            dpi=dq[3*N+4*k:3*N+4*k+4]
        
            dq_jak.append(dpi[0])
            dq_jak.append(dpi[1])
            dq_jak.append(dpi[2])
            dq_jak.append(dpi[3])
        
        dq_jak = np.array(dq_jak).reshape([-1,1])


        Ps = np.concatenate((self.gammaK(q,dq), self.gammaP(dq)), axis=0)
        
        #obliczenie wspolczynnikow metody Baugmarte'a
        Baug = -2*self.alfa*self.jakobianKP(q).dot(dq_jak) - self.beta*self.beta*self.wiezyKP(q)

        Ps += Baug
    
        Pstr = np.concatenate((FA, Ilo, Ps), axis=0)
    
        return Pstr
            
    # integracja metoda Newtona        
    def sym2(self, y0,t0,tK,dt, alfa, beta):
        N = self.N
        self.Y = []
        self.alfa = alfa
        self.beta = beta
        y=y0         
        self.Y.append(copy.copy(y))
        
        print('--trwa obliczanie symulacji--')
        
        for t in np.arange(t0,tK,dt):
            #print("t",t)
            q = y[0:7*N]
            dq = y[7*N:14*N]
            
            
            
            LS = self.Lstrona(q,dq)
            PS = self.Pstrona(q,dq)
            #print(np.linalg.matrix_rank(LS))
            ddq = np.linalg.solve(LS, PS).ravel()
            

            
            for i in range(0,7*N):
                dq[i] += ddq[i]*dt
                q[i] += dq[i]*dt

            for j in range(0,7*N):
                y[j] = q[j]
                y[7*N+j] = dq[j]

            self.Y.append(copy.copy(y))
    
        self.Y = np.array(self.Y)
        

    # integracja funkcja 'ode' metoda 'dopri5'
    def sym(self,y0,t0,tK,dt,alfa,beta):
        N = self.N
        self.Y = []
        self.alfa = alfa
        self.beta = beta
        
        def funode(t,y):
            N = self.N
            q = y[0:7*N]
            dq = y[7*N:14*N]
            
            LS = self.Lstrona(q,dq)
            PS = self.Pstrona(q,dq)
 
            ddq = np.linalg.solve(LS, PS).ravel()  
            ddq = ddq[0:7*N]

            ret = list(dq)+list(ddq)
            return ret

        solver = ode(funode).set_integrator('dopri5')
        solver.set_initial_value(y0, t0)

        while solver.successful() and solver.t < tK:
            print(solver.t)
            solver.integrate(solver.t+dt)
            self.Y.append(solver.y)

        self.Y = np.array(self.Y)

    # rozwiazanie nieliniowego ukladu rownan
    # metoda Newtona-Raphsona
    # dla obliczenia warunkow pocz.
    def newraph(self,q0):
        q=q0
        F=self.wiezyKPD(q)
        it=1
        norma =[]
        norma.append(np.linalg.norm(F))
        while (np.linalg.norm(F)>0.001 and it<25):
            F=self.wiezyKPD(q)
            J=self.jakobianKPD(q)
            
            dq = np.linalg.solve(J,F).ravel()
            
            r = []
            p = []
            for i in range(0, self.N):
            
                ri = dq[7*i:7*i+3]
                pi = dq[7*i+3:7*i+7]
                r.append(ri[0])
                r.append(ri[1])
                r.append(ri[2])
                
                p.append(pi[0])
                p.append(pi[1])
                p.append(pi[2])
                p.append(pi[3])
        
            dq2 = np.concatenate((r,p))
            
            q=q-dq2
            norma.append(np.linalg.norm(F))
            it=it+1  
            if(it==24):
                print('cos nie tak')  
                
        print('warunki początkowe obliczone')
        return q        


    def zapiszWyniki(self, nazwaPliku):
        if len(self.Y) == 0:
            raise Exception('brak wynikow symulacji, najpierw zasymuluj uklad')
        else:
            np.savetxt(nazwaPliku, self.Y, delimiter=";")     
    

        





        