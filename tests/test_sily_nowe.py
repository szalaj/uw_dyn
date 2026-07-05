# -*- coding: utf-8 -*-
# testy nowych typow sil: SilaWPunkcie, SilaKontaktu, lina (tylko_rozciaganie)

import numpy as np
import pytest

from uw_dyn import (Uklad, Czlon, SilaWPunkcie, SilaKontaktu, SilaWewnProst,
                    wektor, u2p)
from conftest import GRAWITACJA


def cialo_swobodne(masa=2.0, J=None, z0=1.0, kat_y=0.0):
    """Pojedynczy czlon bez wiezow kinematycznych."""
    ukl = Uklad()
    ukl.dodajCzlon(Czlon(1, masa, J if J is not None else np.diag([0.5, 0.5, 0.5])))
    ukl.grawitacja = True
    q0 = np.zeros(7)
    q0[2] = z0
    q0[3:7] = u2p(np.array([0.0, 1.0, 0.0]), kat_y)
    return ukl, np.concatenate((q0, np.zeros(7)))


def test_sila_w_punkcie_zawis():
    """Ciag rowny ciezarowi w srodku masy: cialo wisi nieruchomo."""
    masa = 2.0
    ukl, y0 = cialo_swobodne(masa=masa)
    ukl.dodajSileWewn(SilaWPunkcie(1, wektor(0, 0, 0),
                                   wektor(0, 0, masa * GRAWITACJA)))
    ukl.sym2(y0, 0.0, 1.0, 0.005)
    assert np.allclose(ukl.Y[-1][0:3], [0, 0, 1.0], atol=1e-9)
    assert np.allclose(ukl.Y[-1][7:14], 0.0, atol=1e-9)


def test_sila_w_punkcie_podaza_za_orientacja():
    """Cialo obrocone o 90 st. wokol y: ciag 'w osi z ciala' pcha w osi x."""
    masa = 2.0
    ukl, y0 = cialo_swobodne(masa=masa, kat_y=np.pi / 2)
    ukl.grawitacja = False
    T = 4.0
    ukl.dodajSileWewn(SilaWPunkcie(1, wektor(0, 0, 0), wektor(0, 0, T)))
    t = 0.5
    ukl.sym2(y0, 0.0, t, 0.001)
    # os z ciala po obrocie o +90 st. wokol y pokazuje +x globalne
    assert ukl.Y[-1][0] == pytest.approx(0.5 * T / masa * t**2, rel=0.01)
    assert abs(ukl.Y[-1][2] - 1.0) < 1e-6


def test_sila_w_punkcie_moment():
    """Ciag zaczepiony z ramieniem d daje przyspieszenie katowe tau/J."""
    masa, Jyy, d, T = 2.0, 0.5, 0.2, 3.0
    ukl, y0 = cialo_swobodne(masa=masa, J=np.diag([0.5, Jyy, 0.5]))
    ukl.grawitacja = False
    ukl.dodajSileWewn(SilaWPunkcie(1, wektor(d, 0, 0), wektor(0, 0, T)))
    t = 0.05  # krotko: kat maly, orientacja ~stala
    ukl.sym2(y0, 0.0, t, 0.0005)
    # moment s' x f' = (d,0,0)x(0,0,T) = (0, -d*T, 0)
    kat = 2 * np.arctan2(ukl.Y[-1][5], ukl.Y[-1][3])
    kat_teoria = 0.5 * (-d * T / Jyy) * t**2
    assert kat == pytest.approx(kat_teoria, rel=0.02)


def test_kontakt_osiadanie():
    """Cialo upuszczone na podloze osiada na wnikaniu mg/k."""
    masa, k, c = 2.0, 2.0e4, 400.0
    ukl, y0 = cialo_swobodne(masa=masa, z0=0.05)
    ukl.dodajSileWewn(SilaKontaktu(1, wektor(0, 0, 0), k=k, c=c, mu=0.5))
    ukl.sym2(y0, 0.0, 2.0, 0.0005)
    z_rownowagi = -masa * GRAWITACJA / k
    assert ukl.Y[-1][2] == pytest.approx(z_rownowagi, abs=1e-4)
    assert abs(ukl.Y[-1][9]) < 1e-3  # predkosc pionowa wygaszona


def test_kontakt_tarcie_hamuje():
    """Cialo slizgajace sie po podlozu traci predkosc pozioma."""
    masa = 2.0
    ukl, y0 = cialo_swobodne(masa=masa, z0=0.0)
    ukl.dodajSileWewn(SilaKontaktu(1, wektor(0, 0, 0), k=2.0e4, c=400.0, mu=0.5))
    y0[7] = 2.0  # predkosc pozioma
    ukl.sym2(y0, 0.0, 1.0, 0.0005)
    assert ukl.Y[-1][7] < 0.2  # niemal zatrzymane
    # kierunek ruchu bez zmian (tarcie nie odwraca predkosci)
    assert ukl.Y[-1][7] > -1e-6


def test_lina_luzna_nie_pcha():
    """Lina (tylko_rozciaganie): luzna nie dziala, napieta trzyma."""
    masa, k, l0 = 2.0, 5.0e3, 1.0
    # zaczep liny w (0,0,0), cialo tuz pod nim: lina luzna -> swobodny spadek
    ukl, y0 = cialo_swobodne(masa=masa, z0=-0.2)
    ukl.dodajSileWewn(SilaWewnProst(0, 1, wektor(0, 0, 0), wektor(0, 0, 0),
                                    k, l0, 50.0, 0, tylko_rozciaganie=True))
    t = 0.1
    ukl.sym2(y0, 0.0, t, 0.001)
    assert ukl.Y[-1][2] == pytest.approx(-0.2 - 0.5 * GRAWITACJA * t**2, rel=0.02)

    # po opadnieciu ponizej l0 lina sie napina i zatrzymuje cialo
    ukl2, y02 = cialo_swobodne(masa=masa, z0=-0.2)
    ukl2.dodajSileWewn(SilaWewnProst(0, 1, wektor(0, 0, 0), wektor(0, 0, 0),
                                     k, l0, 200.0, 0, tylko_rozciaganie=True))
    ukl2.sym2(y02, 0.0, 3.0, 0.001)
    z_rownowagi = -(l0 + masa * GRAWITACJA / k)
    assert ukl2.Y[-1][2] == pytest.approx(z_rownowagi, abs=1e-3)
