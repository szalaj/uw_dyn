# -*- coding: utf-8 -*-
# testy walidacyjne na wahadle fizycznym:
# rownowaga, okres drgan, energia, wiezy

import numpy as np
import pytest

from conftest import zbuduj_wahadlo, energia_calkowita, normy_kwaternionow, GRAWITACJA


def test_rownowaga_w_spoczynku():
    """Wahadlo wiszace pionowo w spoczynku nie powinno sie ruszyc."""
    ukl, y0 = zbuduj_wahadlo(kat0=0.0)
    ukl.sym2(y0.copy(), 0.0, 2.0, 0.01, 1, 1)
    koniec = ukl.Y[-1]
    assert np.allclose(koniec[0:3], [0.0, 0.0, -2.0], atol=1e-9)
    assert np.allclose(koniec[7:14], 0.0, atol=1e-9)


def test_okres_drgan_male_wychylenie():
    """Okres malych drgan zgodny ze wzorem T = 2*pi*sqrt(I_p/(m*g*L))."""
    masa, J_yy, L = 1.0, 10.0, 2.0
    ukl, y0 = zbuduj_wahadlo(masa=masa, J_yy=J_yy, dlugosc=L, kat0=0.1)
    dt = 0.005
    ukl.sym2(y0.copy(), 0.0, 12.0, dt, 1, 1)

    # okres z przejsc x przez zero (interpolacja liniowa)
    x = ukl.Y[:, 0]
    t = np.arange(len(x)) * dt
    przejscia = []
    for k in range(len(x) - 1):
        if x[k] * x[k + 1] < 0:
            przejscia.append(t[k] + dt * x[k] / (x[k] - x[k + 1]))
    assert len(przejscia) >= 4, "za malo przejsc przez zero"
    # kolejne przejscia przez zero sa co pol okresu
    T_zmierzony = 2 * np.mean(np.diff(przejscia))

    I_przegub = J_yy + masa * L ** 2
    T_analityczny = 2 * np.pi * np.sqrt(I_przegub / (masa * GRAWITACJA * L))
    assert T_zmierzony == pytest.approx(T_analityczny, rel=0.02)


def test_zachowanie_energii():
    """Bez tlumienia energia mechaniczna powinna byc w przyblizeniu stala."""
    ukl, y0 = zbuduj_wahadlo(kat0=0.5)
    ukl.sym2(y0.copy(), 0.0, 10.0, 0.005, 1, 1)
    E0 = energia_calkowita(ukl, ukl.Y[0])
    energie = [energia_calkowita(ukl, ukl.Y[k]) for k in range(0, len(ukl.Y), 50)]
    # skala odniesienia: m*g*L, bo E0 moze byc bliskie zeru
    skala = GRAWITACJA * 2.0
    assert max(abs(E - E0) for E in energie) / skala < 0.05


def test_sym2_nie_modyfikuje_y0():
    ukl, y0 = zbuduj_wahadlo(kat0=0.3)
    kopia = y0.copy()
    ukl.sym2(y0, 0.0, 0.5, 0.01, 1, 1)
    assert np.array_equal(y0, kopia)


def test_normy_kwaternionow_dokladne():
    """Po normalizacji w sym2 normy kwaternionow sa jednostkowe."""
    ukl, y0 = zbuduj_wahadlo(kat0=0.5)
    ukl.sym2(y0, 0.0, 10.0, 0.005, 1, 1)
    for n in normy_kwaternionow(ukl, ukl.Y[-1]):
        assert n == pytest.approx(1.0, abs=1e-12)


def test_sym_adaptacyjny_okres_i_energia():
    """Integrator adaptacyjny (solve_ivp): okres i energia dokladniejsze."""
    masa, J_yy, L = 1.0, 10.0, 2.0
    ukl, y0 = zbuduj_wahadlo(masa=masa, J_yy=J_yy, dlugosc=L, kat0=0.1)
    dt = 0.01
    ukl.sym(y0, 0.0, 12.0, dt, 1, 1)

    x = ukl.Y[:, 0]
    t = np.arange(len(x)) * dt
    przejscia = []
    for k in range(len(x) - 1):
        if x[k] * x[k + 1] < 0:
            przejscia.append(t[k] + dt * x[k] / (x[k] - x[k + 1]))
    T_zmierzony = 2 * np.mean(np.diff(przejscia))
    I_przegub = J_yy + masa * L ** 2
    T_analityczny = 2 * np.pi * np.sqrt(I_przegub / (masa * GRAWITACJA * L))
    assert T_zmierzony == pytest.approx(T_analityczny, rel=0.005)

    E0 = energia_calkowita(ukl, ukl.Y[0])
    energie = [energia_calkowita(ukl, ukl.Y[k]) for k in range(0, len(ukl.Y), 50)]
    assert max(abs(E - E0) for E in energie) / (GRAWITACJA * L) < 1e-4


def test_wiezy_spelnione():
    """Wiezy kinematyczne i norma kwaternionu spelnione w trakcie ruchu."""
    ukl, y0 = zbuduj_wahadlo(kat0=0.5)
    ukl.sym2(y0.copy(), 0.0, 10.0, 0.005, 1, 1)
    for k in range(0, len(ukl.Y), 100):
        q = ukl.Y[k][0:7]
        assert np.linalg.norm(ukl.wiezyK(q)) < 5e-3
        for n in normy_kwaternionow(ukl, ukl.Y[k]):
            assert n == pytest.approx(1.0, abs=1e-3)
