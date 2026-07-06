# -*- coding: utf-8 -*-
# testy czlonu calkujacego aktuatorow (regulator PID)
# PID znosi blad ustalony pod stalym obciazeniem (grawitacja), ktory PD zostawia

import numpy as np
import pytest

from uw_dyn import (Uklad, Czlon, Polaczenie_Obr, Para_Sferyczna,
                    MomentWzgledny, MomentSferyczny, wektor, u2p)


def _wahadlo(ki, theta_cel=0.8, T=12.0):
    """Wahadlo pod grawitacja, slaba sprezyna (widoczny sag); PID gdy ki>0."""
    ukl = Uklad()
    ukl.dodajCzlon(Czlon(1, 1.0, np.diag([1.0, 1.0, 1.0])))
    ukl.dodajWiez(Polaczenie_Obr(0, 1, wektor(0, 0, 0), wektor(0, 0, 2),
                                 wektor(1, 0, 0), wektor(0, 0, 1), wektor(0, 1, 0)))
    ukl.dodajSileWewn(MomentWzgledny(0, 1, wektor(0, 1, 0), wektor(1, 0, 0),
                                     40.0, theta_cel, 8.0, ki=ki, calka_max=200.0))
    ukl.grawitacja = True
    q0 = np.zeros(7)
    q0[2] = -2
    q0[3] = 1
    ukl.sym2(np.concatenate((q0, np.zeros(7))), 0.0, T, 0.002)
    p = ukl.Y[-1][3:7]
    return 2*np.arctan2(p[2], p[0])


def test_pd_ma_blad_ustalony():
    """PD (ki=0) pod grawitacja nie dochodzi do celu (sag > 0.2 rad)."""
    assert 0.8 - _wahadlo(0.0) > 0.2


def test_pid_znosi_blad_ustalony():
    """PID (ki>0) dochodzi do celu mimo grawitacji (blad < 0.06 rad)."""
    assert abs(0.8 - _wahadlo(60.0)) < 0.06


def test_pid_zgodnosc_wsteczna():
    """Bez ki (domyslnie 0) aktuator to czysty PD - zachowanie bez zmian."""
    ukl = Uklad()
    ukl.dodajCzlon(Czlon(1, 1.0, np.diag([1.0, 1.0, 1.0])))
    ukl.dodajWiez(Polaczenie_Obr(0, 1, wektor(0, 0, 0), wektor(0, 0, 2),
                                 wektor(1, 0, 0), wektor(0, 0, 1), wektor(0, 1, 0)))
    a = MomentWzgledny(0, 1, wektor(0, 1, 0), wektor(1, 0, 0), 40.0, 0.5, 8.0)
    assert a.ki == 0.0
    ukl.dodajSileWewn(a)
    ukl.grawitacja = False
    ukl.sym2(np.concatenate((np.array([0, 0, -2, 1, 0, 0, 0.0]), np.zeros(7))),
             0.0, 8.0, 0.005)
    p = ukl.Y[-1][3:7]
    assert 2*np.arctan2(p[2], p[0]) == pytest.approx(0.5, abs=0.02)


def _ramie_sferyczne(ki, T=8.0):
    """Ramie na stawie kulistym, cel poziomy, pod grawitacja."""
    ukl = Uklad()
    ukl.dodajCzlon(Czlon(1, 2.0, np.diag([0.02, 0.05, 0.05])))
    ukl.dodajWiez(Para_Sferyczna(0, 1, wektor(0, 0, 1.4), wektor(0, 0, 0.15)))
    p_cel = u2p(np.array([0.0, 1.0, 0.0]), -np.pi/2)
    akt = MomentSferyczny(0, 1, 120.0, 6.0, p_cel=p_cel, ki=ki, calka_max=60.0)
    ukl.dodajSileWewn(akt)
    ukl.grawitacja = True
    q0 = np.zeros(7)
    q0[2] = 1.4 - 0.15
    q0[3] = 1
    ukl.sym2(np.concatenate((q0, np.zeros(7))), 0.0, T, 0.002)
    return akt.kat(ukl.Y[-1][0:7], 1)


def test_pid_sferyczny_zmniejsza_blad():
    """PID w stawie kulistym zmniejsza blad ustalony orientacji vs samo PD."""
    blad_pd = _ramie_sferyczne(0.0)
    blad_pid = _ramie_sferyczne(40.0)
    assert blad_pid < 0.5*blad_pd
