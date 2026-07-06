# -*- coding: utf-8 -*-
# testy antropometrii: parametry segmentow (Winter) i builder pelnej sylwetki
# Etap A: pelna sylwetka skalowana wzrostem i masa

import numpy as np
import pytest

from uw_dyn.antropometria import (segmenty, segment, masa_calkowita,
                                  tensor_segmentu, zbuduj_postac, TABELA_WINTER)


def test_suma_mas_rowna_masie_ciala():
    """Masy segmentow sumuja sie do calkowitej masy ciala (tablica Wintera)."""
    for M in (55.0, 75.0, 95.0):
        assert masa_calkowita(M, 1.80) == pytest.approx(M, rel=1e-9)


def test_skalowanie_wzrostem_i_masa():
    """Dlugosci skaluja sie wzrostem, masy masa (liniowo)."""
    s1 = segment('udo', 70.0, 1.70)
    s2 = segment('udo', 70.0, 1.90)
    assert s2.dlugosc/s1.dlugosc == pytest.approx(1.90/1.70)
    s3 = segment('udo', 90.0, 1.70)
    assert s3.masa/s1.masa == pytest.approx(90.0/70.0)


def test_tensor_dodatnio_okreslony():
    """Tensory bezwladnosci wszystkich segmentow sa dodatnio okreslone."""
    for s in segmenty(75.0, 1.80).values():
        w = np.linalg.eigvalsh(s.tensor)
        assert (w > 0).all()
        # nierownosc trojkata bezwladnosci (osie poprzeczne >= podluzna)
        assert w[0] > 0


def test_liczba_segmentow():
    """12 segmentow: glowa, tulow + 5 parzystych x 2 strony."""
    s = segmenty(75.0, 1.80)
    assert len(s) == 12
    assert 'ramie_L' in s and 'ramie_P' in s and 'stopa_L' in s


def test_tensor_segmentu_wzor():
    """Tensor: osie poprzeczne m*(rg*L)^2, os podluzna 0.5*m*prom^2."""
    m, L, rg, prom = 2.0, 0.3, 0.3, 0.03
    T = tensor_segmentu(m, L, rg, prom)
    assert T[0, 0] == pytest.approx(m*(rg*L)**2)
    assert T[2, 2] == pytest.approx(0.5*m*prom**2)


def test_builder_sklada_poprawnie():
    """Builder tworzy 12 czlonow, 48 wiezow kinematycznych, 12 aktuatorow."""
    ukl, nry, q0, akt = zbuduj_postac(75.0, 1.80)
    assert ukl.N == 12
    assert ukl.M == 48          # bilans stopni swobody (patrz PLAN)
    assert len(akt) == 12
    assert len(q0) == 14*ukl.N // 2  # 7*N wspolrzednych


def test_postac_trzyma_poze_bez_grawitacji():
    """Bez grawitacji aktuatory trzymaja poze neutralna idealnie (spojnosc)."""
    ukl, nry, q0, akt = zbuduj_postac(75.0, 1.80)
    ukl.grawitacja = False
    N = ukl.N
    ukl.sym2(np.concatenate((q0, np.zeros(7*N))), 0.0, 0.1, 0.0005)
    dryf = max(np.linalg.norm(ukl.Y[-1][3*(nry[n]-1):3*(nry[n]-1)+3]
                              - q0[3*(nry[n]-1):3*(nry[n]-1)+3]) for n in nry)
    assert dryf < 1e-3


def test_postac_na_stopach_sklada_sie():
    """Tryb 'stopy': swobodna podstawa, brak pinu miednicy, kontakt stop."""
    ukl, nry, q0, akt = zbuduj_postac(75.0, 1.80, podparcie='stopy')
    assert ukl.N == 12
    assert ukl.M == 45          # 48 - 3 (bez stawu kulistego w miednicy)
    kontakty = [s for s in ukl.silyWewn if type(s).__name__ == 'SilaKontaktu']
    assert len(kontakty) == 8   # 4 punkty pod kazda stopa


def test_stopy_podpieraja_ciezar():
    """Wcisniete stopy daja pionowa sile kontaktu podpierajaca ciezar ciala."""
    ukl, nry, q0, akt = zbuduj_postac(75.0, 1.80, podparcie='stopy')
    N = ukl.N
    kontakty = [s for s in ukl.silyWewn if type(s).__name__ == 'SilaKontaktu']
    dq = np.zeros(7*N)

    # bez wciiescia (sola na z=0): brak sily
    Fz0 = sum(s.sila(q0, dq, N)[2][2, 0] for s in kontakty)
    assert Fz0 == pytest.approx(0.0, abs=1e-6)

    # obnizenie calej postaci o 5 mm: stopy wnikaja, sila w gore ~ 8*k*delta
    delta = 0.005
    q = q0.copy()
    for i in range(N):
        q[3*i+2] -= delta
    Fz = sum(s.sila(q, dq, N)[2][2, 0] for s in kontakty)
    assert Fz > 0                                  # sila skierowana w gore
    assert Fz == pytest.approx(8*4.0e4*delta, rel=0.05)   # zgodnie z modelem


def test_postac_stoi_pod_grawitacja():
    """Pod grawitacja postac stoi w pozie neutralnej (maly krok dt=1e-4,
    uklad sztywny: ciezki tulow + wiele stawow)."""
    ukl, nry, q0, akt = zbuduj_postac(75.0, 1.80)
    N = ukl.N
    ukl.sym2(np.concatenate((q0, np.zeros(7*N))), 0.0, 0.04, 0.0001)
    assert not np.isnan(ukl.Y).any()
    dryf = max(np.linalg.norm(ukl.Y[-1][3*(nry[n]-1):3*(nry[n]-1)+3]
                              - q0[3*(nry[n]-1):3*(nry[n]-1)+3]) for n in nry)
    assert dryf < 0.02          # stoi (dryf < 2 cm)
