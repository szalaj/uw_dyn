# -*- coding: utf-8 -*-
# testy stabilizacji rzutowaniem i metod energii

import numpy as np
import pytest

from conftest import zbuduj_wahadlo, energia_calkowita, GRAWITACJA


def test_rzutowanie_wiezy_scisle():
    """Stabilizacja rzutowaniem: wiezy spelnione na poziomie tolerancji Newtona."""
    ukl, y0 = zbuduj_wahadlo(kat0=0.5)
    ukl.sym2(y0, 0.0, 5.0, 0.005)
    for k in range(0, len(ukl.Y), 100):
        q = ukl.Y[k][0:7]
        assert np.linalg.norm(ukl.wiezyKP(q)) < 1e-9


def test_rzutowanie_energia_i_okres():
    """Rzutowanie nie psuje fizyki: energia stala, okres jak w Baumgarte."""
    masa, J_yy, L = 1.0, 10.0, 2.0
    ukl, y0 = zbuduj_wahadlo(masa=masa, J_yy=J_yy, dlugosc=L, kat0=0.1)
    dt = 0.005
    ukl.sym2(y0, 0.0, 12.0, dt)

    E0 = ukl.energia(ukl.Y[0])
    energie = [ukl.energia(ukl.Y[k]) for k in range(0, len(ukl.Y), 50)]
    assert max(abs(E - E0) for E in energie) / (GRAWITACJA * L) < 0.05

    x = ukl.Y[:, 0]
    t = np.arange(len(x)) * dt
    przejscia = []
    for k in range(len(x) - 1):
        if x[k] * x[k + 1] < 0:
            przejscia.append(t[k] + dt * x[k] / (x[k] - x[k + 1]))
    T_zmierzony = 2 * np.mean(np.diff(przejscia))
    I_przegub = J_yy + masa * L ** 2
    T_analityczny = 2 * np.pi * np.sqrt(I_przegub / (masa * GRAWITACJA * L))
    assert T_zmierzony == pytest.approx(T_analityczny, rel=0.02)


def test_projekcja_polozen_naprawia_wiezy():
    """Zaburzone polozenie wraca na rozmaitosc wiezow."""
    ukl, y0 = zbuduj_wahadlo(kat0=0.3)
    q = y0[0:7].copy()
    q[0] += 0.05   # zaburzenie lamiace przegub
    q[3] += 0.02   # zaburzenie normy kwaternionu
    assert np.linalg.norm(ukl.wiezyKP(q)) > 1e-2
    q_popr = ukl.projekcja_polozen(q)
    assert np.linalg.norm(ukl.wiezyKP(q_popr)) < 1e-10


def test_projekcja_predkosci_uderzenie():
    """Rzutowanie predkosci: znosi skladowa lamiaca wiezy i nie dodaje energii."""
    ukl, y0 = zbuduj_wahadlo(kat0=0.3)
    q = y0[0:7]
    # predkosc radialna (wzdluz preta): niedopuszczalna dla przegubu
    dq = np.zeros(7)
    dq[0], dq[2] = 0.5, 0.5

    dq_po = ukl.projekcja_predkosci(q, dq)

    # wiezy predkosciowe spelnione: J * dq_jak = 0
    Jq = ukl._jakobianKP_q(q)
    assert np.linalg.norm(Jq.dot(dq_po)) < 1e-10

    Ek_przed = ukl.energia_kinetyczna(np.concatenate((q, dq)))
    Ek_po = ukl.energia_kinetyczna(np.concatenate((q, dq_po)))
    assert Ek_po <= Ek_przed + 1e-12


def test_energia_zgodna_z_definicja():
    """Metody energii Uklad zgodne z niezalezna implementacja w testach."""
    ukl, y0 = zbuduj_wahadlo(kat0=0.4)
    ukl.sym2(y0, 0.0, 1.0, 0.005)
    for k in (0, len(ukl.Y) // 2, -1):
        y = ukl.Y[k]
        assert ukl.energia(y) == pytest.approx(energia_calkowita(ukl, y), abs=1e-9)
