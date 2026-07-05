# -*- coding: utf-8 -*-
# test regresyjny: lancuch czterech czlonow (konfiguracja z przykladu lancuch02)

import numpy as np
import pytest

from uw_dyn import Uklad, Czlon, Polaczenie_Obr, SilaWewnProst, wektor
from conftest import normy_kwaternionow


@pytest.fixture
def lancuch():
    J1 = np.diag([10.0, 10.0, 10.0])
    ukl = Uklad()
    for i in range(1, 5):
        ukl.dodajCzlon(Czlon(i, 1, J1))

    ukl.dodajWiez(Polaczenie_Obr(0, 1, wektor(0, 0, 0), wektor(0, 0, 2),
                                 wektor(1, 0, 0), wektor(0, 0, 1), wektor(0, 1, 0)))
    for i in range(1, 4):
        ukl.dodajWiez(Polaczenie_Obr(i, i + 1, wektor(0, 0, -2), wektor(0, 0, 2),
                                     wektor(1, 0, 0), wektor(0, 0, 1), wektor(0, 1, 0)))

    # elementy sprezyste i wymuszenie jak w przykladzie lancuch02
    k, c, f1, t = 0, 0, 40, 1
    ukl.dodajSileWewn(SilaWewnProst(0, 1, wektor(-t, 0, 0), wektor(-t, 0, -2), k, 4, c, f1))
    for i in range(1, 4):
        ukl.dodajSileWewn(SilaWewnProst(i, i + 1, wektor(-t, 0, 2), wektor(-t, 0, -2), k, 4, c, 0))

    ukl.grawitacja = True

    q0 = np.zeros(7 * ukl.N)
    q0[2], q0[5], q0[8], q0[11] = -2, -6, -10, -14
    q0[12] = q0[16] = q0[20] = q0[24] = 1
    y0 = np.concatenate((q0, np.zeros(7 * ukl.N)))
    return ukl, y0


def test_lancuch_stabilny(lancuch):
    """Krotka symulacja lancucha: brak NaN, wiezy i normy kwaternionow OK."""
    ukl, y0 = lancuch
    ukl.sym2(y0.copy(), 0.0, 2.0, 0.01, 1, 1)

    assert not np.isnan(ukl.Y).any()
    assert ukl.Y.shape == (201, 14 * ukl.N)

    q_koniec = ukl.Y[-1][0:7 * ukl.N]
    assert np.linalg.norm(ukl.wiezyK(q_koniec)) < 1e-2
    for n in normy_kwaternionow(ukl, ukl.Y[-1]):
        assert n == pytest.approx(1.0, abs=5e-3)


def test_lancuch_reaguje_na_sile(lancuch):
    """Wymuszenie f1 dziala w kierunku x: lancuch musi sie wychylic."""
    ukl, y0 = lancuch
    ukl.sym2(y0.copy(), 0.0, 2.0, 0.01, 1, 1)
    x4 = ukl.Y[-1][9]  # wspolrzedna x czlonu 4
    assert abs(x4) > 0.01
