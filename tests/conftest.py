# -*- coding: utf-8 -*-
# wspolne funkcje pomocnicze testow

import numpy as np
import pytest

from uw_dyn import Uklad, Czlon, Polaczenie_Obr, wektor, u2p, G, p_i, dp_i, dr_i

GRAWITACJA = 9.80665


def zbuduj_wahadlo(masa=1.0, J_yy=10.0, dlugosc=2.0, kat0=0.0):
    """Wahadlo fizyczne: jeden czlon zawieszony przegubem obrotowym w poczatku
    ukladu, os obrotu y, srodek masy w odleglosci `dlugosc` od przegubu.
    Zwraca (uklad, y0) dla kata poczatkowego `kat0` (obrot wokol osi y)."""
    J = np.diag([J_yy, J_yy, J_yy])
    ukl = Uklad()
    ukl.dodajCzlon(Czlon(1, masa, J))
    ukl.dodajWiez(Polaczenie_Obr(0, 1,
                                 wektor(0, 0, 0), wektor(0, 0, dlugosc),
                                 wektor(1, 0, 0), wektor(0, 0, 1), wektor(0, 1, 0)))
    ukl.grawitacja = True

    # warunki poczatkowe spelniajace wiezy: r = -R(p) * sB
    p0 = u2p(np.array([0.0, 1.0, 0.0]), kat0)
    q0 = np.zeros(7)
    q0[0] = -dlugosc * np.sin(kat0)
    q0[2] = -dlugosc * np.cos(kat0)
    q0[3:7] = p0
    dq0 = np.zeros(7)
    y0 = np.concatenate((q0, dq0))
    return ukl, y0


def energia_calkowita(ukl, y):
    """Energia mechaniczna ukladu dla wektora stanu y = [q, dq]."""
    N = ukl.N
    q = y[0:7 * N]
    dq = y[7 * N:14 * N]
    E = 0.0
    for cz in ukl.czlony:
        i = cz.i
        dr = dr_i(i, dq)
        p = p_i(i, q, N)
        dp = dp_i(i, dq, N)
        om = 2 * G(p).dot(dp)  # predkosc katowa w ukladzie ciala
        E += 0.5 * cz.m * dr.T.dot(dr).item()
        E += 0.5 * om.T.dot(cz.J).dot(om).item()
        E += cz.m * GRAWITACJA * float(q[3 * (i - 1) + 2])
    return E


def normy_kwaternionow(ukl, y):
    """Normy kwaternionow wszystkich czlonow dla wektora stanu y."""
    N = ukl.N
    q = y[0:7 * N]
    return [float(np.linalg.norm(q[3 * N + 4 * k:3 * N + 4 * k + 4])) for k in range(N)]
