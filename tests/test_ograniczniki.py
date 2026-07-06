# -*- coding: utf-8 -*-
# testy miekkich ogranicznikow zakresu ruchu w stawach
# Etap A: OgranicznikKata (zawias) i OgranicznikStozka (staw kulisty)

import numpy as np
import pytest

from uw_dyn import (Uklad, Czlon, Polaczenie_Obr, Para_Sferyczna,
                    MomentWzgledny, MomentSferyczny,
                    OgranicznikKata, OgranicznikStozka, wektor, u2p, p_i, R)


def _wahadlo_z_napedem(theta_cel, kat_min, kat_max, k_ogr=5000.0):
    """Czlon na przegubie obrotowym (os y), naped do theta_cel + ogranicznik."""
    ukl = Uklad()
    ukl.dodajCzlon(Czlon(1, 1.0, np.diag([1.0, 1.0, 1.0])))
    ukl.dodajWiez(Polaczenie_Obr(0, 1, wektor(0, 0, 0), wektor(0, 0, 2),
                                 wektor(1, 0, 0), wektor(0, 0, 1), wektor(0, 1, 0)))
    ukl.dodajSileWewn(MomentWzgledny(0, 1, wektor(0, 1, 0), wektor(1, 0, 0),
                                     50.0, theta_cel, 10.0))
    ukl.dodajSileWewn(OgranicznikKata(0, 1, wektor(0, 1, 0), wektor(1, 0, 0),
                                      kat_min, kat_max, k_ogr, 40.0))
    ukl.grawitacja = False
    q0 = np.zeros(7)
    q0[2] = -2
    q0[3] = 1
    ukl.sym2(np.concatenate((q0, np.zeros(7))), 0.0, 8.0, 0.004)
    return ukl


def _kat_y(q):
    p = q[3:7]
    return 2*np.arctan2(p[2], p[0])


def test_ogranicznik_kata_gorna_granica():
    """Naped ciagnie do 1.0, ogranicznik 0.5: przegub staje przy 0.5."""
    ukl = _wahadlo_z_napedem(1.0, -0.5, 0.5)
    katy = np.array([_kat_y(w[0:7]) for w in ukl.Y])
    assert katy.max() < 0.56           # nie przekracza istotnie granicy
    assert _kat_y(ukl.Y[-1][0:7]) == pytest.approx(0.5, abs=0.05)


def test_ogranicznik_kata_dolna_granica():
    """Naped ciagnie do -1.2, ogranicznik -0.4: przegub staje przy -0.4."""
    ukl = _wahadlo_z_napedem(-1.2, -0.4, 0.9)
    katy = np.array([_kat_y(w[0:7]) for w in ukl.Y])
    # miekki ogranicznik dopuszcza male przejsciowe wniknienie pod silnym napedem
    assert katy.min() > -0.5
    assert _kat_y(ukl.Y[-1][0:7]) == pytest.approx(-0.4, abs=0.05)


def test_ogranicznik_kata_w_zakresie_nieaktywny():
    """Cel 0.3 wewnatrz zakresu [-0.5,0.5]: przegub osiaga cel bez zmian."""
    ukl = _wahadlo_z_napedem(0.3, -0.5, 0.5)
    assert _kat_y(ukl.Y[-1][0:7]) == pytest.approx(0.3, abs=1e-3)


def test_ogranicznik_stozka_trzyma_w_stozku():
    """Naped kulisty przechyla os o 1.0 rad, stozek 0.5: os zostaje w stozku."""
    ukl = Uklad()
    ukl.dodajCzlon(Czlon(1, 1.0, np.diag([0.5, 0.5, 0.5])))
    ukl.dodajWiez(Para_Sferyczna(0, 1, wektor(0, 0, 0), wektor(0, 0, -1)))
    ukl.dodajSileWewn(MomentSferyczny(0, 1, 40.0, 10.0,
                                      p_cel=u2p(np.array([0.0, 1.0, 0.0]), 1.0)))
    ogr = OgranicznikStozka(0, 1, wektor(0, 0, 1), 0.5, 3000.0, 40.0)
    ukl.dodajSileWewn(ogr)
    ukl.grawitacja = False
    q0 = np.zeros(7)
    q0[2] = 1.0
    q0[3] = 1.0
    ukl.sym2(np.concatenate((q0, np.zeros(7))), 0.0, 8.0, 0.004)

    katy = np.array([ogr.kat(w[0:7], 1) for w in ukl.Y])
    assert katy.max() < 0.56
    assert ogr.kat(ukl.Y[-1][0:7], 1) == pytest.approx(0.5, abs=0.06)


def test_ogranicznik_stozka_w_zakresie_nieaktywny():
    """Przechylenie 0.3 rad < stozek 0.6: staw osiaga cel bez ograniczania."""
    ukl = Uklad()
    ukl.dodajCzlon(Czlon(1, 1.0, np.diag([0.5, 0.5, 0.5])))
    ukl.dodajWiez(Para_Sferyczna(0, 1, wektor(0, 0, 0), wektor(0, 0, -1)))
    ukl.dodajSileWewn(MomentSferyczny(0, 1, 40.0, 10.0,
                                      p_cel=u2p(np.array([0.0, 1.0, 0.0]), 0.3)))
    ogr = OgranicznikStozka(0, 1, wektor(0, 0, 1), 0.6, 3000.0, 40.0)
    ukl.dodajSileWewn(ogr)
    ukl.grawitacja = False
    q0 = np.zeros(7)
    q0[2] = 1.0
    q0[3] = 1.0
    ukl.sym2(np.concatenate((q0, np.zeros(7))), 0.0, 8.0, 0.004)
    assert ogr.kat(ukl.Y[-1][0:7], 1) == pytest.approx(0.3, abs=0.02)
