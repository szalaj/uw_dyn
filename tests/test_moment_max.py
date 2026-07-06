# -*- coding: utf-8 -*-
# testy ograniczenia maksymalnego momentu aktuatorow (saturacja)
# Etap A: sily napedow stawow fizycznie sensowne (nie nieskonczona sprezyna)

import numpy as np
import pytest

from uw_dyn import (Uklad, Czlon, Polaczenie_Obr, Para_Sferyczna,
                    MomentWzgledny, MomentSferyczny, wektor, u2p)

GRAW = 9.80665


def _wahadlo(moment_max, theta_cel=1.5, masa=1.0, L=2.0, grawitacja=True):
    """Wahadlo (przegub obrotowy os y) z napedem o ograniczonym momencie."""
    ukl = Uklad()
    ukl.dodajCzlon(Czlon(1, masa, np.diag([1.0, 1.0, 1.0])))
    ukl.dodajWiez(Polaczenie_Obr(0, 1, wektor(0, 0, 0), wektor(0, 0, L),
                                 wektor(1, 0, 0), wektor(0, 0, 1), wektor(0, 1, 0)))
    ukl.dodajSileWewn(MomentWzgledny(0, 1, wektor(0, 1, 0), wektor(1, 0, 0),
                                     200.0, theta_cel, 20.0, moment_max=moment_max))
    ukl.grawitacja = grawitacja
    q0 = np.zeros(7)
    q0[2] = -L
    q0[3] = 1
    ukl.sym2(np.concatenate((q0, np.zeros(7))), 0.0, 12.0, 0.003)
    p = ukl.Y[-1][3:7]
    return 2*np.arctan2(p[2], p[0])


def _wahadlo_maks_kat(moment_max, theta_cel=1.5, masa=1.0, L=2.0):
    """Jak _wahadlo, ale zwraca maksymalny kat osiagniety w trakcie ruchu."""
    ukl = Uklad()
    ukl.dodajCzlon(Czlon(1, masa, np.diag([1.0, 1.0, 1.0])))
    ukl.dodajWiez(Polaczenie_Obr(0, 1, wektor(0, 0, 0), wektor(0, 0, L),
                                 wektor(1, 0, 0), wektor(0, 0, 1), wektor(0, 1, 0)))
    ukl.dodajSileWewn(MomentWzgledny(0, 1, wektor(0, 1, 0), wektor(1, 0, 0),
                                     200.0, theta_cel, 20.0, moment_max=moment_max))
    ukl.grawitacja = True
    q0 = np.zeros(7)
    q0[2] = -L
    q0[3] = 1
    ukl.sym2(np.concatenate((q0, np.zeros(7))), 0.0, 12.0, 0.003)
    katy = np.array([2*np.arctan2(w[5], w[3]) for w in ukl.Y])
    return katy.max()


def test_moment_max_ogranicza_zasieg():
    """Slaby naped nie jest w stanie dzwignac ciezaru do celu (limit momentu
    ponizej momentu grawitacji przy celu); silny naped cel osiaga."""
    masa, L, cel = 1.0, 2.0, 1.5
    slaby = 0.5*masa*GRAW*L         # < moment grawitacji przy celu (m*g*L*sin1.5)
    maks_slaby = _wahadlo_maks_kat(slaby, theta_cel=cel, masa=masa, L=L)
    maks_silny = _wahadlo_maks_kat(5000.0, theta_cel=cel, masa=masa, L=L)
    assert maks_slaby < cel - 0.3      # slaby nie dochodzi do celu (saturacja)
    assert maks_silny > cel - 0.05     # silny cel osiaga


def test_moment_max_duzy_osiaga_cel():
    """Duzy limit momentu: naped bez trudu osiaga cel (saturacja nieaktywna)."""
    theta = _wahadlo(moment_max=500.0, theta_cel=1.0, grawitacja=False)
    assert theta == pytest.approx(1.0, abs=0.02)


def test_moment_max_brak_limitu_zgodnosc():
    """Bez moment_max (None) zachowanie jak dawniej: osiaga cel."""
    theta = _wahadlo(moment_max=None, theta_cel=1.0, grawitacja=False)
    assert theta == pytest.approx(1.0, abs=0.02)


def test_moment_max_sferyczny_slaby_nie_utrzyma():
    """Slaby staw kulisty nie utrzyma ramienia w celu (grawitacja wygrywa),
    silny utrzyma. Porownanie potwierdza dzialanie limitu normy momentu."""
    def tilt_koncowy(moment_max):
        ukl = Uklad()
        ukl.dodajCzlon(Czlon(1, 1.0, np.diag([0.05, 0.05, 0.05])))
        # COM ponizej stawu (stabilny zwis): staw w origin, punkt ciala (0,0,1)
        ukl.dodajWiez(Para_Sferyczna(0, 1, wektor(0, 0, 0), wektor(0, 0, 1)))
        # cel: przechyl w poziom (obrot o pi/2 wokol y)
        ukl.dodajSileWewn(MomentSferyczny(0, 1, 300.0, 20.0,
                                          p_cel=u2p(np.array([0.0, 1.0, 0.0]), np.pi/2),
                                          moment_max=moment_max))
        ukl.grawitacja = True
        q0 = np.zeros(7)
        q0[2] = -1.0
        q0[3] = 1
        ukl.sym2(np.concatenate((q0, np.zeros(7))), 0.0, 10.0, 0.002)
        p = ukl.Y[-1][3:7]
        # kat przechylenia lokalnej osi z od pionu
        from uw_dyn import R, wektor_p
        zloc = R(wektor_p(*p))[:, 2]
        return float(np.arccos(np.clip(zloc[2], -1, 1)))

    slaby = tilt_koncowy(3.0)
    silny = tilt_koncowy(300.0)
    assert silny > 1.3          # prawie poziom (pi/2 ~ 1.57)
    assert slaby < silny - 0.3  # slaby wyraznie nie dochodzi (saturacja)
