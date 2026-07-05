# -*- coding: utf-8 -*-
# uw_dyn: dynamika 3D ukladow wieloczlonowych
# autor: Marcin Szalajski (praca magisterska, 2016; pakiet od 2026)

"""Sily: elementy sprezysto-tlumiace i sily zewnetrzne."""

import numpy as np

from uw_dyn.algebra import r_i, dr_i, p_i, dp_i, R, G, skew


class SilaWewnProst:
    """Element sprezysto-tlumiacy z sila stala (sprezyna k o dlugosci
    swobodnej l0, tlumik c, sila F) miedzy punktem A czlonu i oraz punktem B
    czlonu j. Z flaga tylko_rozciaganie=True dziala jak lina: sila znika,
    gdy element jest krotszy niz dlugosc swobodna."""
    def __init__(self, i, j, sA_i, sB_j, k,l0, c, F, tylko_rozciaganie=False):
            self.i = i
            self.j = j
            self.sA_i = sA_i
            self.sB_j = sB_j
            self.k = k
            self.l0 = l0
            self.c = c
            self.F = F
            self.tylko_rozciaganie = tylko_rozciaganie

    def energia_potencjalna(self,q,N):
        """Energia sprezysta elementu (czlon tlumika i sily stalej pominiety)."""
        if self.k == 0:
            return 0.0
        l = self.dlugosc(q,N)
        if self.tylko_rozciaganie and l < self.l0:
            return 0.0
        return 0.5*self.k*(l - self.l0)**2
            
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

            if self.tylko_rozciaganie and l < self.l0:
                return Qr_i, Qp_i, Qr_j, Qp_j

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

            if self.tylko_rozciaganie and l < self.l0:
                return Qr_i, Qp_i, Qr_j, Qp_j

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
        
        
class SilaWPunkcie:
    """Sila zaczepiona w punkcie ciala, zadana w ukladzie ciala
    (follower force), np. ciag wirnika drona.

    czlon: numer czlonu; s_punkt: punkt zaczepienia w ukladzie ciala;
    f_lokalna: wektor sily w ukladzie ciala (podazy za orientacja).
    Wektor sily mozna podmieniac miedzy segmentami symulacji
    (dyskretny regulator)."""

    def __init__(self, czlon, s_punkt, f_lokalna):
        self.i = 0
        self.j = czlon
        self.s_punkt = s_punkt
        self.f_lokalna = f_lokalna

    def energia_potencjalna(self,q,N):
        return 0.0

    def sila(self,q,dq,N):
        pj = p_i(self.j,q,N)
        Rj = R(pj)
        Gj = G(pj)

        Qr_i = np.zeros((3,1))
        Qp_i = np.zeros((4,1))
        # sila globalna: R f'; moment w ukladzie ciala: s' x f'
        Qr_j = Rj.dot(self.f_lokalna)
        n_lok = skew(self.s_punkt).dot(self.f_lokalna)
        Qp_j = 2*Gj.transpose().dot(n_lok)
        return Qr_i, Qp_i, Qr_j, Qp_j


class SilaKontaktu:
    """Jednostronny kontakt punktu ciala z podlozem z = 0 (model penalty).

    Gdy punkt zaczepienia jest pod powierzchnia: sila normalna
    N = max(0, k*wnikanie - c*predkosc_pionowa) oraz tarcie styczne
    (regularyzowany Coulomb): T = -mu*N * v_t/(|v_t| + eps).

    czlon: numer czlonu; s_punkt: punkt stopy w ukladzie ciala;
    k, c: sztywnosc i tlumienie podloza; mu: wspolczynnik tarcia."""

    def __init__(self, czlon, s_punkt, k=2.0e4, c=200.0, mu=0.8, eps=0.01):
        self.i = 0
        self.j = czlon
        self.s_punkt = s_punkt
        self.k = k
        self.c = c
        self.mu = mu
        self.eps = eps

    def _punkt(self,q,N):
        rj = r_i(self.j,q)
        pj = p_i(self.j,q,N)
        Rj = R(pj)
        return rj + Rj.dot(self.s_punkt), Rj

    def energia_potencjalna(self,q,N):
        pw, _ = self._punkt(q,N)
        wnikanie = -float(pw[2,0])
        if wnikanie <= 0:
            return 0.0
        return 0.5*self.k*wnikanie**2

    def sila(self,q,dq,N):
        Qr_i = np.zeros((3,1))
        Qp_i = np.zeros((4,1))
        Qr_j = np.zeros((3,1))
        Qp_j = np.zeros((4,1))

        pw, Rj = self._punkt(q,N)
        wnikanie = -float(pw[2,0])
        if wnikanie <= 0:
            return Qr_i, Qp_i, Qr_j, Qp_j

        # predkosc punktu stopy: dr + R (om' x s')
        drj = dr_i(self.j,dq)
        pj = p_i(self.j,q,N)
        dpj = dp_i(self.j,dq,N)
        Gj = G(pj)
        om_j = 2*Gj.dot(dpj)
        v = drj + Rj.dot(skew(om_j)).dot(self.s_punkt)

        Fn = max(0.0, self.k*wnikanie - self.c*float(v[2,0]))
        vt = np.array([float(v[0,0]), float(v[1,0])])
        Ft = -self.mu*Fn*vt/(np.linalg.norm(vt) + self.eps)

        F_glob = np.array([[Ft[0]], [Ft[1]], [Fn]])
        Qr_j = F_glob
        # moment w ukladzie ciala: s' x (R^T F)
        n_lok = skew(self.s_punkt).dot(Rj.transpose().dot(F_glob))
        Qp_j = 2*Gj.transpose().dot(n_lok)
        return Qr_i, Qp_i, Qr_j, Qp_j


class SilaZewn:
    """Sila (Fx/Fy/Fz) lub moment (nx/ny/nz) o stalej wielkosci dzialajacy na czlon."""
    def __init__(self,czlon, rodzaj, wielkosc):
        self.czlon = czlon
        self.rodzaj = rodzaj
        self.wielkosc = wielkosc


