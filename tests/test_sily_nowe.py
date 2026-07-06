# -*- coding: utf-8 -*-
# testy nowych typow sil: SilaWPunkcie, SilaKontaktu, lina (tylko_rozciaganie)

import numpy as np
import pytest

from uw_dyn import (Uklad, Czlon, Polaczenie_Obr, Para_Sferyczna, SilaWPunkcie,
                    SilaKontaktu, SilaUderzenia, SilaWewnProst, MomentWzgledny,
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


@pytest.mark.parametrize("theta_cel", [0.3, -0.5, 0.8])
def test_moment_wzgledny_osiaga_kat(theta_cel):
    """Aktuator obrotowy w przegubie ustawia człon na zadanym kącie."""
    ukl = Uklad()
    ukl.dodajCzlon(Czlon(1, 1.0, np.diag([1., 1., 1.])))
    ukl.dodajWiez(Polaczenie_Obr(0, 1, wektor(0, 0, 0), wektor(0, 0, 2),
                                 wektor(1, 0, 0), wektor(0, 0, 1), wektor(0, 1, 0)))
    ukl.dodajSileWewn(MomentWzgledny(0, 1, wektor(0, 1, 0), wektor(1, 0, 0),
                                     k=50.0, theta_cel=theta_cel, c=15.0))
    ukl.grawitacja = False
    q0 = np.zeros(7)
    q0[2] = -2
    q0[3] = 1
    ukl.sym2(np.concatenate((q0, np.zeros(7))), 0.0, 8.0, 0.005)
    p = ukl.Y[-1][3:7]
    assert 2 * np.arctan2(p[2], p[0]) == pytest.approx(theta_cel, abs=1e-3)


def test_moment_wzgledny_para_wewnetrzna():
    """Moment jest parą wewnętrzną: nie zmienia pędu ani nie działa dla i=0
    poza swoim członem (zwięzły test bilansu na dwóch wolnych członach)."""
    ukl = Uklad()
    ukl.dodajCzlon(Czlon(1, 1.0, np.diag([1., 1., 1.])))
    ukl.dodajCzlon(Czlon(2, 1.0, np.diag([1., 1., 1.])))
    ukl.dodajWiez(Polaczenie_Obr(1, 2, wektor(0, 0, 1), wektor(0, 0, -1),
                                 wektor(1, 0, 0), wektor(0, 0, 1), wektor(0, 1, 0)))
    akt = MomentWzgledny(1, 2, wektor(0, 1, 0), wektor(1, 0, 0),
                         k=10.0, theta_cel=0.6, c=2.0)
    ukl.dodajSileWewn(akt)
    ukl.grawitacja = False
    q0 = np.zeros(14)
    q0[2] = 1
    q0[5] = -1
    q0[6] = 1
    q0[10] = 1
    ukl.sym2(np.concatenate((q0, np.zeros(14))), 0.0, 6.0, 0.005)
    # kąt względny człon2-człon1 dąży do celu
    assert akt.kat(ukl.Y[-1][0:14], 2) == pytest.approx(0.6, abs=2e-2)


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


def test_uderzenie_rownowaga_sil_i_potencjal():
    """SilaUderzenia (kontakt brya-brya): sila na oba ciaa jest rowna i
    przeciwna, a jest gradientem energii penalty (Q = -dV/dq stycznie)."""
    N, EPS = 2, 1e-6
    rng = np.random.default_rng(1)
    s = SilaUderzenia(1, wektor(0.1, 0, 0), 2, promien=0.3, polowa_wys=0.2,
                      k=1.0e4, c=0.0)

    def losq():
        q = np.zeros(7 * N)
        q[0:3 * N] = rng.normal(size=3 * N) * 0.1
        for k in range(N):
            p = rng.normal(size=4)
            q[3 * N + 4 * k:3 * N + 4 * k + 4] = p / np.linalg.norm(p)
        q[0:3] = q[3:6] + rng.normal(size=3) * 0.05     # wymus penetracje
        return q

    maxF, maxV, traf = 0.0, 0.0, 0
    for _ in range(40):
        q = losq()
        Qr_i, Qp_i, Qr_j, Qp_j = s.sila(q, np.zeros(7 * N), N)
        maxF = max(maxF, np.linalg.norm(np.ravel(Qr_i) + np.ravel(Qr_j)))
        if s.energia_potencjalna(q, N) > 0:
            traf += 1
            Q = np.zeros(7 * N)
            Q[0:3] = np.ravel(Qr_i)
            Q[3 * N:3 * N + 4] = np.ravel(Qp_i)
            Q[3:6] = np.ravel(Qr_j)
            Q[3 * N + 4:3 * N + 8] = np.ravel(Qp_j)
            for _ in range(4):
                d = np.zeros(7 * N)
                d[0:3 * N] = rng.normal(size=3 * N)
                for k in range(N):
                    p = q[3 * N + 4 * k:3 * N + 4 * k + 4]
                    dp = rng.normal(size=4)
                    dp -= p * np.dot(p, dp)
                    d[3 * N + 4 * k:3 * N + 4 * k + 4] = dp
                d /= np.linalg.norm(d)
                dV = (s.energia_potencjalna(q + EPS * d, N)
                      - s.energia_potencjalna(q - EPS * d, N)) / (2 * EPS)
                maxV = max(maxV, abs(dV + Q @ d))
    assert traf > 20                 # wiekszosc prob to trafienia
    assert maxF < 1e-9               # rownowaga sil (F + (-F) = 0)
    assert maxV < 1e-4               # Q = -dV/dq


def test_uderzenie_wprawia_worek_w_ruch():
    """Worek jako wahadlo (staw kulisty) uderzony poruszajacym sie cialem
    wychyla sie; w spoczynku (bez uderzenia) wisi nieruchomo."""
    def scena(v_uderzenia):
        ukl = Uklad()
        ukl.dodajCzlon(Czlon(1, 1.0, np.diag([0.05, 0.05, 0.05])))   # piesc
        ukl.dodajCzlon(Czlon(2, 5.0, np.diag([0.1, 0.1, 0.02])))     # worek
        # worek: wahadlo na stawie kulistym, zaczep nad workiem
        ukl.dodajWiez(Para_Sferyczna(0, 2, wektor(0.3, 0, 1.4),
                                     wektor(0, 0, 0.4)))
        ukl.dodajSileWewn(SilaUderzenia(1, wektor(0, 0, 0), 2,
                                        promien=0.2, polowa_wys=0.4,
                                        k=2.0e4, c=50.0))
        ukl.grawitacja = True
        q0 = np.zeros(14)
        q0[0:3] = [0.05, 0, 1.0]      # piesc tuz przed workiem
        q0[3:6] = [0.3, 0, 1.0]       # worek (COM pod zaczepem)
        q0[6] = 1.0
        q0[10] = 1.0
        dq0 = np.zeros(14)
        dq0[0] = v_uderzenia          # piesc leci w worek (+x)
        ukl.sym2(np.concatenate((q0, dq0)), 0.0, 0.6, 5e-4)
        worek = ukl.Y[:, 3:6]
        return np.hypot(worek[:, 0] - 0.3, worek[:, 1]).max()

    assert scena(6.0) > 0.03          # uderzony worek sie wychyla
    assert scena(0.0) < 5e-3          # bez uderzenia wisi (prawie) nieruchomo
