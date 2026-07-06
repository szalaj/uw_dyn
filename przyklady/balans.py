# -*- coding: utf-8 -*-
# Przyklad (Etap B): postac stoi na stopach z regulatorem balansu PID.
#
# Pelna sylwetka (antropometria) na swobodnej podstawie, stopy w kontakcie
# z podlozem. Stawy trzymane regulatorami PID (czlon calkujacy znosi sag pod
# grawitacja -> mniejsza sztywnosc wystarcza, wiec wiekszy krok dt). Na to
# nalozony jest REGULATOR BALANSU PID: sprzezenie od poziomego polozenia
# srodka masy (CoM) wzgledem srodka stop do:
#   - kostek (przod-tyl, os y) - strategia kostki,
#   - bioder (bok, przechyl wokol x) - strategia biodra.
# Postac dostaje pchniecie (poczatkowa predkosc pozioma) i regulator utrzymuje
# ja nad wielobokiem podparcia (bez balansu przewraca sie znacznie bardziej).
#
# Wynik: web/dane_balans.js do wizualizacji Three.js (web/balans.html).

import json
import os

import numpy as np

from uw_dyn.antropometria import zbuduj_postac, segmenty
from uw_dyn import u2p, mnoz_kwaterniony

OSX = np.array([1.0, 0.0, 0.0])

# ----- parametry -----
MASA, WZROST = 75.0, 1.80
PCHNIECIE = 0.35        # pozioma predkosc poczatkowa (pchniecie) [m/s]
SEGMENT = 0.01          # takt regulatora balansu [s]
DT = 2.0e-4             # krok calkowania [s]
CZAS = 1.0

# PID trzymania pozy (stawy): nizsza sztywnosc + czlon calkujacy (bez sagu)
K_HOLD, C_HOLD, KI_HOLD, CMAX_HOLD = 140.0, 14.0, 180.0, 120.0

# PID balansu (CoM -> kostki/biodra)
KP_BAL, KI_BAL, KD_BAL = 10.0, 12.0, 2.2
WZM_BIODRO = 1.6        # wzmocnienie kanalu bocznego (biodra)


def _com(ukl, fr):
    m = 0.0
    c = np.zeros(3)
    for cz in ukl.czlony:
        i = cz.i
        c += cz.m*fr[3*(i-1):3*(i-1)+3]
        m += cz.m
    return c/m


def symuluj(balans=True, czas=CZAS):
    ukl, nry, q0, akt = zbuduj_postac(MASA, WZROST, podparcie='stopy')
    N = ukl.N
    for a in akt.values():          # stawy jako PID trzymajacy poze
        a.k, a.c, a.ki, a.calka_max = K_HOLD, C_HOLD, KI_HOLD, CMAX_HOLD

    baza_kostka = {b: akt['kostka_'+b].theta_cel for b in 'LP'}
    baza_biodro = {b: akt['biodro_'+b].p_cel.copy() for b in 'LP'}
    srodek = np.mean([q0[3*(nry['stopa_'+b]-1):3*(nry['stopa_'+b]-1)+2]
                      for b in 'LP'], axis=0)

    q = q0.copy()
    dq = np.zeros(7*N)
    for i in range(N):              # pchniecie: pozioma predkosc wszystkich czlonow
        dq[3*i] = PCHNIECIE

    klatki = []
    c_prev = _com(ukl, q)
    calka = np.zeros(2)
    maks_dryf = 0.0
    t = 0.0
    while t < czas:
        c = _com(ukl, q)
        v = (c - c_prev)/SEGMENT if t > 0 else np.zeros(3)
        c_prev = c
        e = c[0:2] - srodek
        maks_dryf = max(maks_dryf, float(np.hypot(*e)))

        if balans:
            calka += e*SEGMENT
            u = KP_BAL*e + KI_BAL*calka + KD_BAL*v[0:2]
            for b in 'LP':
                akt['kostka_'+b].theta_cel = baza_kostka[b] - u[0]   # przod-tyl
                akt['biodro_'+b].p_cel = mnoz_kwaterniony(
                    u2p(OSX, WZM_BIODRO*u[1]), baza_biodro[b])        # bok

        ukl.sym2(np.concatenate((q, dq)), 0.0, SEGMENT, DT)
        Y = ukl.Y
        for w in Y[:-1]:
            klatki.append(w[0:7*N].copy())
        q = Y[-1][0:7*N].copy()
        dq = Y[-1][7*N:14*N].copy()
        t += SEGMENT

    return ukl, nry, klatki, maks_dryf


def eksportuj(ukl, nry, klatki, co_ile=5, plik='web/dane_balans.js'):
    N = ukl.N
    inv = {v: k for k, v in nry.items()}
    S = segmenty(MASA, WZROST)
    wymiary = {inv[i+1]: [round(S[inv[i+1]].dlugosc, 4), round(S[inv[i+1]].promien, 4)]
               for i in range(N)}
    dane_klatki = []
    for q in klatki[::co_ile]:
        czlony = []
        for k in range(N):
            r = q[3*k:3*k+3]
            p = q[3*N + 4*k:3*N + 4*k + 4]
            czlony.append({'r': [round(float(v), 4) for v in r],
                           'p': [round(float(v), 5) for v in p]})
        dane_klatki.append(czlony)
    dane = {'dt': DT*co_ile,
            'kolejnosc': [inv[i+1] for i in range(N)],
            'wymiary': wymiary,
            'klatki': dane_klatki}
    katalog = os.path.dirname(os.path.abspath(plik))
    os.makedirs(katalog, exist_ok=True)
    with open(plik, 'w') as f:
        f.write('const DANE = ' + json.dumps(dane) + ';\n')
    print(f'zapisano {len(dane_klatki)} klatek do {plik}')


if __name__ == '__main__':
    _, _, _, dryf_bez = symuluj(balans=False, czas=CZAS)
    ukl, nry, klatki, dryf_pid = symuluj(balans=True, czas=CZAS)
    print(f'pchniecie {PCHNIECIE} m/s, maks. poziomy dryf CoM w {CZAS} s:')
    print(f'  bez balansu: {dryf_bez:.3f} m')
    print(f'  PID balans:  {dryf_pid:.3f} m')
    sciezka = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           '..', 'web', 'dane_balans.js')
    eksportuj(ukl, nry, klatki, plik=os.path.normpath(sciezka))
