# -*- coding: utf-8 -*-
# testy pozostalych wiezow i sil:
# para przesuwna, sila zewnetrzna, sprezyna/tlumik, wiezy kierujace (newraph)

import numpy as np
import pytest

from uw_dyn import (Uklad, Czlon, Polaczenie_Obr, Polaczenie_Przes,
                    Odleglosc, Kat, SilaWewnProst, SilaZewn, wektor, u2p)
from conftest import zbuduj_wahadlo, energia_calkowita, GRAWITACJA


def zbuduj_suwak_pionowy(masa=1.0, z0=-1.5):
    """Czlon na parze przesuwnej wzdluz osi z (1 stopien swobody)."""
    ukl = Uklad()
    ukl.dodajCzlon(Czlon(1, masa, np.diag([1.0, 1.0, 1.0])))
    ukl.dodajWiez(Polaczenie_Przes(0, 1, wektor(0, 0, 0), wektor(0, 0, 0),
                                   wektor(1, 0, 0), wektor(0, 1, 0),
                                   wektor(0, 1, 0), wektor(0, 0, 1)))
    ukl.grawitacja = True
    q0 = np.zeros(7)
    q0[2] = z0
    q0[3] = 1.0
    y0 = np.concatenate((q0, np.zeros(7)))
    return ukl, y0


def test_przesuwna_spadek_swobodny():
    """Para przesuwna wzdluz z: swobodny spadek z(t) = z0 - g*t^2/2."""
    ukl, y0 = zbuduj_suwak_pionowy(z0=0.0)
    ukl.grawitacja = True
    ukl.sym(y0, 0.0, 1.0, 0.01, 1, 1)
    t = 1.0
    z_analityczne = -0.5 * GRAWITACJA * t ** 2
    koniec = ukl.Y[-1]
    assert koniec[2] == pytest.approx(z_analityczne, rel=1e-5)
    # pozostale wspolrzedne bez zmian
    assert np.allclose(koniec[0:2], 0.0, atol=1e-9)
    assert np.allclose(koniec[3:7], [1, 0, 0, 0], atol=1e-6)


def test_sila_zewnetrzna_rownowaga_momentu():
    """Moment ny rownowazacy grawitacje: wahadlo stoi w wychyleniu."""
    masa, L, kat = 1.0, 2.0, 0.2
    ukl, y0 = zbuduj_wahadlo(masa=masa, dlugosc=L, kat0=kat)
    ukl.dodajSileZewn(SilaZewn(1, 'ny', masa * GRAWITACJA * L * np.sin(kat)))
    ukl.sym2(y0.copy(), 0.0, 2.0, 0.005, 1, 1)
    assert np.allclose(ukl.Y[-1][0:3], y0[0:3], atol=1e-6)
    assert np.allclose(ukl.Y[-1][7:14], 0.0, atol=1e-6)


def test_sprezyna_rownowaga_statyczna():
    """Suwak z sprezyna i tlumikiem osiada w x_eq = l0 + m*g/k."""
    masa, k, c, l0 = 1.0, 100.0, 10.0, 1.0
    ukl, y0 = zbuduj_suwak_pionowy(masa=masa, z0=-1.5)
    ukl.dodajSileWewn(SilaWewnProst(0, 1, wektor(0, 0, 0), wektor(0, 0, 0),
                                    k, l0, c, 0))
    ukl.sym2(y0, 0.0, 8.0, 0.002, 1, 1)
    z_eq = -(l0 + masa * GRAWITACJA / k)
    assert ukl.Y[-1][2] == pytest.approx(z_eq, abs=1e-4)
    assert abs(ukl.Y[-1][9]) < 1e-4  # predkosc dz ~ 0


def test_sprezyna_czestosc_drgan():
    """Bez tlumienia okres drgan suwaka T = 2*pi*sqrt(m/k)."""
    masa, k, l0 = 1.0, 100.0, 1.0
    z_eq = -(l0 + masa * GRAWITACJA / k)
    ukl, y0 = zbuduj_suwak_pionowy(masa=masa, z0=z_eq - 0.1)
    ukl.dodajSileWewn(SilaWewnProst(0, 1, wektor(0, 0, 0), wektor(0, 0, 0),
                                    k, l0, 0, 0))
    dt = 0.001
    ukl.sym2(y0, 0.0, 3.0, dt, 1, 1)

    z = ukl.Y[:, 2] - z_eq
    t = np.arange(len(z)) * dt
    przejscia = []
    for i in range(len(z) - 1):
        if z[i] * z[i + 1] < 0:
            przejscia.append(t[i] + dt * z[i] / (z[i] - z[i + 1]))
    assert len(przejscia) >= 4
    T_zmierzony = 2 * np.mean(np.diff(przejscia))
    T_analityczny = 2 * np.pi * np.sqrt(masa / k)
    assert T_zmierzony == pytest.approx(T_analityczny, rel=0.01)


def test_kat_warunki_poczatkowe():
    """Wiez kierujacy Kat + newraph: wahadlo ustawione pod zadanym katem.

    Kat miedzy osia x podstawy a osia z ciala rowny pi/3 oznacza
    wychylenie wahadla o pi/6 od pionu."""
    ukl, _ = zbuduj_wahadlo()
    ukl.dodajWiezD(Kat(0, 1, wektor(1, 0, 0), wektor(0, 0, 1),
                       wektor(0, 0, 1), np.pi / 3))

    kat_zgadywany = 0.4
    q_guess = np.zeros(7)
    q_guess[0] = -2 * np.sin(kat_zgadywany)
    q_guess[2] = -2 * np.cos(kat_zgadywany)
    q_guess[3:7] = u2p(np.array([0.0, 1.0, 0.0]), kat_zgadywany)

    q = ukl.newraph(q_guess)
    assert np.linalg.norm(ukl.wiezyKPD(q)) < 1e-3
    # sin(pi/6) = 0.5 -> x = -2*0.5 = -1, z = -2*cos(pi/6) = -sqrt(3)
    assert q[0] == pytest.approx(-1.0, abs=1e-3)
    assert q[2] == pytest.approx(-np.sqrt(3), abs=1e-3)


def test_odleglosc_warunki_poczatkowe():
    """Wiez kierujacy Odleglosc + newraph: zadana wspolrzedna x wahadla."""
    ukl, _ = zbuduj_wahadlo()
    ukl.dodajWiezD(Odleglosc(0, 1, wektor(0, 0, 0), wektor(0, 0, 0),
                             wektor(1, 0, 0), -1.0))

    kat_zgadywany = 0.4
    q_guess = np.zeros(7)
    q_guess[0] = -2 * np.sin(kat_zgadywany)
    q_guess[2] = -2 * np.cos(kat_zgadywany)
    q_guess[3:7] = u2p(np.array([0.0, 1.0, 0.0]), kat_zgadywany)

    q = ukl.newraph(q_guess)
    assert np.linalg.norm(ukl.wiezyKPD(q)) < 1e-3
    assert q[0] == pytest.approx(-1.0, abs=1e-3)
    assert q[2] == pytest.approx(-np.sqrt(3), abs=1e-3)


def test_wahadlo_podwojne_energia():
    """Wahadlo podwojne (przegub czlon-czlon, i != 0): zachowanie energii.

    Cwiczy czlony sprzezenia predkosci katowych w gammaK par miedzy
    dwoma ruchomymi czlonami."""
    J = np.diag([10.0, 10.0, 10.0])
    ukl = Uklad()
    ukl.dodajCzlon(Czlon(1, 1.0, J))
    ukl.dodajCzlon(Czlon(2, 1.0, J))
    ukl.dodajWiez(Polaczenie_Obr(0, 1, wektor(0, 0, 0), wektor(0, 0, 2),
                                 wektor(1, 0, 0), wektor(0, 0, 1), wektor(0, 1, 0)))
    ukl.dodajWiez(Polaczenie_Obr(1, 2, wektor(0, 0, -2), wektor(0, 0, 2),
                                 wektor(1, 0, 0), wektor(0, 0, 1), wektor(0, 1, 0)))
    ukl.grawitacja = True

    # oba czlony wychylone sztywno o kat0 (lancuch prosty, obrocony)
    kat0 = 0.6
    kierunek = np.array([-np.sin(kat0), 0.0, -np.cos(kat0)])
    p0 = u2p(np.array([0.0, 1.0, 0.0]), kat0)
    q0 = np.zeros(14)
    q0[0:3] = 2 * kierunek
    q0[3:6] = 6 * kierunek
    q0[6:10] = p0
    q0[10:14] = p0
    y0 = np.concatenate((q0, np.zeros(14)))

    ukl.sym2(y0, 0.0, 4.0, 0.002, 1, 1)

    E0 = energia_calkowita(ukl, ukl.Y[0])
    energie = [energia_calkowita(ukl, ukl.Y[k]) for k in range(0, len(ukl.Y), 100)]
    skala = 2 * GRAWITACJA * 2.0  # 2 czlony * m*g*L
    assert max(abs(E - E0) for E in energie) / skala < 0.02
    assert not np.isnan(ukl.Y).any()
