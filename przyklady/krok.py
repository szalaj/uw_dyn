# -*- coding: utf-8 -*-
# Przyklad (Etap B): praca nog - krok w bok.
#
# Rozwiniecie balansu (balans.py) o zmiane podparcia: przeniesienie ciezaru
# na noge podporowa (lewa), uniesienie i odwiedzenie nogi wykrocznej (prawa),
# postawienie jej dalej w bok, wyrownanie. Kontakt penalty sam obsluguje
# oderwanie i postawienie stopy (sila znika, gdy stopa opuszcza podloze).
#
# Sterowanie:
#   - przeniesienie ciezaru: pochylenie bioder w bok (roll) -> CoM nad lewa noga,
#   - noga wykroczna: zgiecie kolana (unos) + odwiedzenie biodra (w bok),
#   - balans przod-tyl: regulator PID na kostkach (jak w balans.py).
#
# STAN (uczciwie): dziala mechanika kroku (przeniesienie ciezaru, uniesienie
# nogi przez zwolnienie kontaktu, wymach, postawienie) ORAZ balans - lateralny
# PID (bok, biodra) + PID przod-tyl (kostki) utrzymuja figure na nogach przez
# faze jednonozna, czego bez balansu brakowalo (przewracala sie).
# Nie osiagnieto jednak czystego, stabilnego kroku KIERUNKOWEGO: silne
# odwiedzenie stawia stope wyraznie w bok, ale figura osiada (kolana ugielaja
# sie pod dynamika); umiarkowane odwiedzenie trzyma balans, lecz pochyl
# rownowazacy sciaga noge z powrotem do srodka. To sprzezenie balans<->krok
# wymaga regulatora ustawienia stopy (capture point) i drobniejszego strojenia.
# Calkowanie: sym3 polniejawny (RATTLE + niejawne tlumienie/sprezyny), dt=1e-3
# (5x wiekszy krok niz sym2). UWAGA: ta sekwencja jest WRAZLIWA na dt - wynik
# koncowy stopy zmienia sie jakosciowo miedzy dt (bo to sterowanie w petli
# otwartej bez regulatora ustawienia stopy), np. dt=2e-3 daje absurdalny
# wyskok. Balans trzyma figure, ale sam krok wymaga jeszcze capture pointu -
# az do tego trzymamy zachowawczy dt=1e-3.
#
# Wynik: web/dane_krok.js do wizualizacji Three.js (web/krok.html).

import json
import os

import numpy as np

from uw_dyn.antropometria import zbuduj_postac, segmenty
from uw_dyn import u2p, mnoz_kwaterniony

OSX = np.array([1.0, 0.0, 0.0])

MASA, WZROST = 75.0, 1.80
SEGMENT = 0.01
DT = 1.0e-3          # sym3 polniejawny; regulator kroku jest wrazliwy na dt
                     # (patrz naglowek), wiec bezpieczny krok - nie 2e-3+
CZAS = 1.8

K_HOLD, C_HOLD, KI_HOLD, CMAX_HOLD = 140.0, 14.0, 180.0, 120.0
KP_ANK, KI_ANK, KD_ANK = 8.0, 10.0, 1.6      # balans przod-tyl (kostki)
KP_LAT, KI_LAT, KD_LAT = 10.0, 12.0, 2.2     # balans bok (biodra, roll)
WZM_BIODRO = 1.6

# fazy kroku w bok (w prawo): przeniesienie ciezaru, unos, odwiedzenie, postawienie
T_PRZENIES = 0.5
T_UNOS = 0.35
T_ODWODZ = 0.35
T_POSTAW = 0.35
CEL_LEWA = 0.10      # docelowe przesuniecie CoM w bok nad lewa noge [m]
CEL_KONC = -0.05     # docelowy CoM w bok po postawieniu (nowa, szersza baza)
KOLANO_UNOS = 0.9    # zgiecie kolana nogi wykrocznej (unos stopy)
ABDUKCJA = 0.6       # odwiedzenie biodra nogi wykrocznej (w bok)


def _com(ukl, fr):
    m = 0.0
    c = np.zeros(3)
    for cz in ukl.czlony:
        i = cz.i
        c += cz.m*fr[3*(i-1):3*(i-1)+3]
        m += cz.m
    return c/m


def _prof(t, t0, dt_f):
    """Gladkie 0->1 na [t0, t0+dt_f]."""
    if t <= t0:
        return 0.0
    if t >= t0 + dt_f:
        return 1.0
    return (1 - np.cos(np.pi*(t - t0)/dt_f))/2


def symuluj():
    ukl, nry, q0, akt = zbuduj_postac(MASA, WZROST, podparcie='stopy')
    N = ukl.N
    for a in akt.values():
        a.k, a.c, a.ki, a.calka_max = K_HOLD, C_HOLD, KI_HOLD, CMAX_HOLD

    baza_kostka = {b: akt['kostka_'+b].theta_cel for b in 'LP'}
    baza_biodro = {b: akt['biodro_'+b].p_cel.copy() for b in 'LP'}
    baza_kolano_P = akt['kolano_P'].theta_cel
    srodek = np.mean([q0[3*(nry['stopa_'+b]-1):3*(nry['stopa_'+b]-1)+2]
                      for b in 'LP'], axis=0)

    # harmonogram faz
    t1 = T_PRZENIES
    t2 = t1 + T_UNOS
    t3 = t2 + T_ODWODZ
    t4 = t3 + T_POSTAW

    q = q0.copy()
    dq = np.zeros(7*N)
    klatki = []
    c_prev = _com(ukl, q)
    calka_x = 0.0
    calka_y = 0.0
    t = 0.0
    while t < CZAS:
        c = _com(ukl, q)
        v = (c - c_prev)/SEGMENT if t > 0 else np.zeros(3)
        c_prev = c

        # balans przod-tyl: PID na kostkach
        ex = c[0] - srodek[0]
        calka_x += ex*SEGMENT
        u_ank = KP_ANK*ex + KI_ANK*calka_x + KD_ANK*v[0]
        for b in 'LP':
            akt['kostka_'+b].theta_cel = baza_kostka[b] - u_ank

        # boczny cel CoM: 0 -> nad lewa noga (single support) -> nowa baza
        cel_y = (srodek[1] + CEL_LEWA*(_prof(t, 0, T_PRZENIES) - _prof(t, t3, T_POSTAW))
                 + (CEL_KONC - srodek[1])*_prof(t, t4, 0.4))
        # balans bok: PID na roll bioder (utrzymuje CoM na celu)
        ey = c[1] - cel_y
        calka_y += ey*SEGMENT
        roll = WZM_BIODRO*(KP_LAT*ey + KI_LAT*calka_y + KD_LAT*v[1])
        for b in 'LP':
            akt['biodro_'+b].p_cel = mnoz_kwaterniony(u2p(OSX, roll), baza_biodro[b])

        # noga wykroczna (prawa): unos kolana + odwiedzenie (na roll nakladamy abd)
        unos = _prof(t, t1, T_UNOS) - _prof(t, t3, T_POSTAW)
        abd = _prof(t, t2, T_ODWODZ) - _prof(t, t3, T_POSTAW)
        akt['kolano_P'].theta_cel = baza_kolano_P + KOLANO_UNOS*unos
        akt['biodro_P'].p_cel = mnoz_kwaterniony(
            u2p(OSX, roll), mnoz_kwaterniony(baza_biodro['P'],
                                             u2p(OSX, ABDUKCJA*abd)))

        ukl.sym3(np.concatenate((q, dq)), 0.0, SEGMENT, DT, polniejawne=True)
        Y = ukl.Y
        for w in Y[:-1]:
            klatki.append(w[0:7*N].copy())
        q = Y[-1][0:7*N].copy()
        dq = Y[-1][7*N:14*N].copy()
        t += SEGMENT

    return ukl, nry, klatki, srodek


def eksportuj(ukl, nry, klatki, co_ile=5, plik='web/dane_krok.js'):
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
    ukl, nry, klatki, srodek = symuluj()
    N = ukl.N
    st0 = klatki[0][3*(nry['stopa_P']-1):3*(nry['stopa_P']-1)+3]
    stK = klatki[-1][3*(nry['stopa_P']-1):3*(nry['stopa_P']-1)+3]
    zt = klatki[-1][3*(nry['tulow']-1)+2]
    print(f'prawa stopa: y {st0[1]:.3f} -> {stK[1]:.3f}')
    print(f'wysokosc tulowia koncowa: {zt:.3f} m (postawa ~1.19); stoi: {zt > 0.85}')
    sciezka = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           '..', 'web', 'dane_krok.js')
    eksportuj(ukl, nry, klatki, plik=os.path.normpath(sciezka))
