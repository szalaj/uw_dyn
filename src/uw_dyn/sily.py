# -*- coding: utf-8 -*-
# uw_dyn: dynamika 3D ukladow wieloczlonowych
# autor: Marcin Szalajski (praca magisterska, 2016; pakiet od 2026)

"""Sily: elementy sprezysto-tlumiace i sily zewnetrzne."""

import numpy as np

from uw_dyn.algebra import r_i, dr_i, p_i, dp_i, R, G, skew


class SilaWewnProst:
    """Element sprezysto-tlumiacy z sila stala (sprezyna k o dlugosci swobodnej l0, tlumik c, sila F) miedzy punktem A czlonu i oraz punktem B czlonu j."""
    def __init__(self, i, j, sA_i, sB_j, k,l0, c, F):
            self.i = i
            self.j = j
            self.sA_i = sA_i
            self.sB_j = sB_j
            self.k = k
            self.l0 = l0
            self.c = c
            self.F = F
            
    def dlugosc(self,q,N):
        """Aktualna dlugosc elementu (odleglosc punktow zaczepienia)."""
        rj = r_i(self.j,q)
        pj = p_i(self.j,q,N)
        Rj = R(pj)
        if self.i == 0:
            dij = rj + Rj.dot(self.sB_j) - self.sA_i
        else:
            ri = r_i(self.i,q)
            pi = p_i(self.i,q,N)
            Ri = R(pi)
            dij = rj + Rj.dot(self.sB_j) - ri - Ri.dot(self.sA_i)
        return float(np.linalg.norm(dij))

    def sila(self,q,dq,N):
        rj = r_i(self.j,q)
        drj = dr_i(self.j,dq)
        
        pj = p_i(self.j,q,N)
        dpj = dp_i(self.j,dq,N)
        Rj = R(pj)
        Gj = G(pj)

        om_j = 2*Gj.dot(dpj)

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
            dpi = dp_i(self.i,dq,N)
            Ri = R(pi)
            Gi = G(pi)


            om_i = 2*Gi.dot(dpi)
            
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
    """Sila (Fx/Fy/Fz) lub moment (nx/ny/nz) o stalej wielkosci dzialajacy na czlon."""
    def __init__(self,czlon, rodzaj, wielkosc):
        self.czlon = czlon
        self.rodzaj = rodzaj
        self.wielkosc = wielkosc


