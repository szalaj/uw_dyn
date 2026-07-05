# -*- coding: utf-8 -*-
# testy napedzanego stawu kulistego 3 DOF (MomentSferyczny)
# Etap A: fundament pod biomechaniczny model czlowieka.

import numpy as np
import pytest

from uw_dyn import (Uklad, Czlon, Para_Sferyczna, MomentSferyczny,
                    wektor, u2p, p_i,
                    mnoz_kwaterniony, sprzezenie_kwaternionu,
                    kwaternion_na_wektor_obrotu)


def _blad_orientacji(p_osiagniete, p_cel):
    err = mnoz_kwaterniony(p_osiagniete, sprzezenie_kwaternionu(p_cel))
    return float(np.linalg.norm(kwaternion_na_wektor_obrotu(err)))


@pytest.mark.parametrize("os,kat", [
    ([0, 0, 1], 0.6), ([0, 1, 0], -0.8), ([1, 0, 0], 1.0),
    ([1, 1, 0], 0.9), ([0.3, 0.5, 0.8], 1.2),
])
def test_osiaga_dowolna_orientacje(os, kat):
    """Staw kulisty sprowadza czlon do zadanej orientacji 3D (dowolna oś)."""
    ukl = Uklad()
    ukl.dodajCzlon(Czlon(1, 1.0, np.diag([0.5, 0.5, 0.5])))
    ukl.dodajWiez(Para_Sferyczna(0, 1, wektor(0, 0, 0), wektor(0, 0, -1)))
    p_cel = u2p(np.asarray(os, float), kat)
    ukl.dodajSileWewn(MomentSferyczny(0, 1, k=40.0, c=12.0, p_cel=p_cel))
    ukl.grawitacja = False
    q0 = np.zeros(7)
    q0[2] = 1.0
    q0[3] = 1.0
    ukl.sym2(np.concatenate((q0, np.zeros(7))), 0.0, 8.0, 0.005)
    p_kon = p_i(1, ukl.Y[-1][0:7], 1).ravel()
    assert np.degrees(_blad_orientacji(p_kon, p_cel)) < 0.5


def test_orientacja_wzgledna_miedzy_czlonami():
    """Cel zadany wzgledem czlonu i (i != 0): sprawdzenie obrotu wzglednego."""
    ukl = Uklad()
    ukl.dodajCzlon(Czlon(1, 1.0, np.diag([0.4, 0.4, 0.4])))
    ukl.dodajCzlon(Czlon(2, 1.0, np.diag([0.4, 0.4, 0.4])))
    # czlon 1 przypiety i sztywno trzymany w obrocie (na skos), czlon 2 na nim
    ukl.dodajWiez(Para_Sferyczna(0, 1, wektor(0, 0, 0), wektor(0, 0, -1)))
    ukl.dodajWiez(Para_Sferyczna(1, 2, wektor(0, 0, 1), wektor(0, 0, -1)))
    p1_cel = u2p(np.array([1.0, 0.0, 0.0]), 0.5)
    ukl.dodajSileWewn(MomentSferyczny(0, 1, k=200.0, c=25.0, p_cel=p1_cel))
    p_rel_cel = u2p(np.array([0.0, 1.0, 0.0]), 0.7)  # czlon 2 wzgledem 1
    akt2 = MomentSferyczny(1, 2, k=60.0, c=15.0, p_cel=p_rel_cel)
    ukl.dodajSileWewn(akt2)
    ukl.grawitacja = False

    q0 = np.zeros(14)
    q0[2] = 1.0
    q0[5] = 2.0
    q0[6] = 1.0
    q0[10] = 1.0
    ukl.sym2(np.concatenate((q0, np.zeros(14))), 0.0, 10.0, 0.004)
    assert akt2.kat(ukl.Y[-1][0:14], 2) == pytest.approx(0.0, abs=2e-2)


def test_trzyma_ramie_pod_grawitacja():
    """Bark (staw kulisty) trzyma ramie w zadanej pozie mimo grawitacji."""
    ukl = Uklad()
    ukl.dodajCzlon(Czlon(1, 2.0, np.diag([0.02, 0.05, 0.05])))
    # bark w (0,0,1.4); srodek ramienia wisi ponizej (wzdluz lokalnej -z)
    ukl.dodajWiez(Para_Sferyczna(0, 1, wektor(0, 0, 1.4), wektor(0, 0, 0.15)))
    # cel: ramie uniesione poziomo (obrot o -90 st. wokol y -> lokalna -z na +x)
    p_cel = u2p(np.array([0.0, 1.0, 0.0]), -np.pi/2)
    akt = MomentSferyczny(0, 1, k=200.0, c=8.0, p_cel=p_cel)
    ukl.dodajSileWewn(akt)
    ukl.grawitacja = True

    q0 = np.zeros(7)
    q0[2] = 1.4 - 0.15
    q0[3] = 1.0
    ukl.sym2(np.concatenate((q0, np.zeros(7))), 0.0, 6.0, 0.002)
    # utrzymuje poze blisko celu mimo ciezaru (maly zwis statyczny)
    assert np.degrees(akt.kat(ukl.Y[-1][0:7], 1)) < 12.0
    assert not np.isnan(ukl.Y).any()


def test_energia_zero_na_celu():
    """Energia potencjalna aktuatora = 0, gdy czlon jest w orientacji docelowej."""
    ukl = Uklad()
    ukl.dodajCzlon(Czlon(1, 1.0, np.diag([0.5, 0.5, 0.5])))
    p_cel = u2p(np.array([0.0, 0.0, 1.0]), 0.4)
    akt = MomentSferyczny(0, 1, k=50.0, c=5.0, p_cel=p_cel)
    q = np.zeros(7)
    q[0:4] = [0, 0, 0, 0]
    q[3:7] = p_cel   # ustaw orientacje = cel
    assert akt.energia_potencjalna(q, 1) == pytest.approx(0.0, abs=1e-12)
    assert akt.kat(q, 1) == pytest.approx(0.0, abs=1e-9)
