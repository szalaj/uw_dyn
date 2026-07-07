# -*- coding: utf-8 -*-
# uw_dyn: dynamika 3D ukladow wieloczlonowych
# autor: Marcin Szalajski (2016; pakiet od 2026)

"""Wiezy: pary kinematyczne i wiezy kierujace."""

import numpy as np

from uw_dyn.algebra import r_i, dr_i, p_i, dp_i, R, G, dG, skew, wektor, wektor_p


class Polaczenie:
    
    def __init__(self):
        self.m = 0
        

        
class Para_Prostopadla(Polaczenie):
    """Wiez prostopadlosci wektorow ai (czlon i) i aj (czlon j): ai . R aj = 0."""
    
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
        Om_j = 2*Gj.dot(dpj)
          
        if self.i == 0:        
            pom = Om_j_skew.dot(Om_j_skew).dot(Rj.transpose())
            gamK = -self.aj.transpose().dot(pom).dot(self.ai)
        else:
            pi = p_i(self.i,q,N)
            dpi = dp_i(self.i,dq,N)
            Gi = G(pi)
            dGi = dG(dpi)
            
            Om_i_skew = 2*Gi.dot(dGi.transpose())
            Om_i = 2*Gi.dot(dpi)
            
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
    """Wiez prostopadlosci wektora ai do wektora laczacego punkty A (czlon i) i B (czlon j)."""
    
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
            ri = r_i(self.i,q)
            rj = r_i(self.j,q)
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
        Om_j = 2*Gj.dot(dpj)
        
        
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
            Om_i = 2*Gi.dot(dpi)
            
            Ri = R(pi)
            Rj = R(pj)
            
            dij = rj + Rj.dot(self.sB_j) - ri - Ri.dot(self.sA_i)
            
            pom1 = 2*Om_i.transpose().dot(skew(self.ai)).dot(Ri.transpose()).dot(dri-drj)
            pom2 = 2*self.sB_j.transpose().dot(Om_j_skew).dot(Rj.transpose()).dot(Ri).dot(Om_i_skew).dot(self.ai)
            pom3 = self.sA_i.transpose().dot(Om_i_skew).dot(Om_i_skew).dot(self.ai)
            pom4 = self.sB_j.transpose().dot(Om_j_skew).dot(Om_j_skew).dot(Rj.transpose()).dot(Ri).dot(self.ai)
            pom5 = dij.transpose().dot(Ri).dot(Om_i_skew).dot(Om_i_skew).dot(self.ai)
            
            gamK = pom1 + pom2 - pom3 - pom4 - pom5
            
        return gamK
        

class Para_Sferyczna(Polaczenie):
    """Przegub kulisty: punkt A czlonu i pokrywa sie z punktem B czlonu j (3 wiezy)."""
    
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
        Om_j = 2*Gj.dot(dpj)
    
        if self.i == 0:        
            
            gamK[0:3] = -Rj.dot(Om_j_skew).dot(Om_j_skew).dot(self.sB_j)

        else:
            
            pi = p_i(self.i,q,N)
            dpi = dp_i(self.i,dq,N)
            
            Gi = G(pi)
            dGi = dG(dpi)
            
            Om_i_skew = 2*Gi.dot(dGi.transpose())
            Om_i = 2*Gi.dot(dpi)
            
            Ri = R(pi)
            Rj = R(pj)
            

            pom1 = Ri.dot(Om_i_skew).dot(Om_i_skew).dot(self.sA_i)
            pom2 = Rj.dot(Om_j_skew).dot(Om_j_skew).dot(self.sB_j)
            
            gamK[0:3] = pom1 - pom2

        
        return gamK  
        
class Polaczenie_Obr(Para_Sferyczna):
    """Przegub obrotowy: przegub kulisty + os obrotu (uj czlonu j prostopadle do vi i wi czlonu i); 5 wiezow."""
    
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

        J1i, J1j = self.pp1.jakobianK(q,N)
        J2i, J2j = self.pp2.jakobianK(q,N)

        Fqi = np.concatenate((Fqi, J1i, J2i), axis=0)
        Fqj = np.concatenate((Fqj, J1j, J2j), axis=0)

        return Fqi, Fqj
        
    def gammaK(self,q,dq,N):
        
        gamK=super().gammaK(q,dq,N)
        
        gamK = np.append(gamK,self.pp1.gammaK(q,dq,N),0)
        gamK = np.append(gamK,self.pp2.gammaK(q,dq,N),0)

        return gamK         
     
  
class Polaczenie_Cyl(Polaczenie):
    """Para cylindryczna: translacja + obrot wzdluz wspolnej osi; 4 wiezy."""
    
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

        J1i, J1j = self.p_prost1.jakobianK(q,N)
        J2i, J2j = self.p_prost2.jakobianK(q,N)
        J3i, J3j = self.p_prost_dij_1.jakobianK(q,N)
        J4i, J4j = self.p_prost_dij_2.jakobianK(q,N)

        Fqi[0,:]= J1i
        Fqi[1,:]= J2i
        Fqi[2,:]= J3i
        Fqi[3,:]= J4i

        Fqj[0,:]= J1j
        Fqj[1,:]= J2j
        Fqj[2,:]= J3j
        Fqj[3,:]= J4j

        return Fqi, Fqj
        
    def gammaK(self,q,dq,N):
    
        gamK=np.zeros([4,1])
        
        gamK[0]= self.p_prost1.gammaK(q,dq,N)
        gamK[1]=  self.p_prost2.gammaK(q,dq,N)
        gamK[2] = self.p_prost_dij_1.gammaK(q,dq,N)
        gamK[3] = self.p_prost_dij_2.gammaK(q,dq,N)
        
        return gamK
        
class Polaczenie_Przes(Polaczenie_Cyl):
    """Para przesuwna (pryzmatyczna): jak cylindryczna, ale bez obrotu; 5 wiezow."""
    
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
        J5i, J5j = self.p_prost3.jakobianK(q,N)
        Fqi = np.append(Fqi, J5i, 0)
        Fqj = np.append(Fqj, J5j, 0)

        return Fqi, Fqj
        
    def gammaK(self,q,dq,N):
    
        gamK=super().gammaK(q,dq,N)
        gamK = np.append(gamK,  self.p_prost3.gammaK(q,dq,N), 0)
        
        return gamK
        
        
class Odleglosc(Para_Prostopadla_D):
    """Wiez kierujacy: zadana odleglosc C punktu B czlonu j od punktu A czlonu i wzdluz kierunku ai (uzywany przez newraph)."""

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
    """Wiez kierujacy: zadany kat fi miedzy wektorem vi czlonu i a wektorem uj czlonu j (uzywany przez newraph)."""
    
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
        Fqi = np.zeros([1,7])
        Fqj = np.zeros([1,7])
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
