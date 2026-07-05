# -*- coding: utf-8 -*-
# uw_dyn: dynamika 3D ukladow wieloczlonowych
# autor: Marcin Szalajski (praca magisterska, 2016; pakiet od 2026)

"""Algebra: wektory, kwaterniony (parametry Eulera), macierze obrotu."""

import numpy as np
from math import cos,sin,pi


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


# iloczyn Hamiltona kwaternionow (parametrow Eulera) [e0,e1,e2,e3]
def mnoz_kwaterniony(pa, pb):
    a = np.asarray(pa, dtype=float).ravel()
    b = np.asarray(pb, dtype=float).ravel()
    return np.array([
        a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3],
        a[0]*b[1] + a[1]*b[0] + a[2]*b[3] - a[3]*b[2],
        a[0]*b[2] - a[1]*b[3] + a[2]*b[0] + a[3]*b[1],
        a[0]*b[3] + a[1]*b[2] - a[2]*b[1] + a[3]*b[0],
    ])


# sprzezenie kwaternionu (obrot odwrotny dla kwaternionu jednostkowego)
def sprzezenie_kwaternionu(p):
    a = np.asarray(p, dtype=float).ravel()
    return np.array([a[0], -a[1], -a[2], -a[3]])


# wektor obrotu (os * kat) z kwaternionu; wybiera obrot krotszy (kat w [0,pi])
def kwaternion_na_wektor_obrotu(p):
    a = np.asarray(p, dtype=float).ravel()
    if a[0] < 0:
        a = -a
    v = a[1:4]
    nv = np.linalg.norm(v)
    if nv < 1e-12:
        return np.zeros(3)
    kat = 2*np.arctan2(nv, a[0])
    return kat*v/nv

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
    b = np.asarray(a).ravel()
    ax = b[0]
    ay = b[1]
    az = b[2]
    A = np.array([[0.0, -az, ay],
    [az, 0.0, -ax],
    [-ay, ax, 0.0]])
    return A
  
# tworzenie wektoru pionowego r    
def wektor(ax,ay,az):
    return np.array([[ax],[ay],[az]])
    
# tworzenie wektoru pionowego p    
def wektor_p(e0,e1,e2,e3):
    return np.array([[e0],[e1],[e2],[e3]])
    
        

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


