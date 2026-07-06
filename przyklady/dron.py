# -*- coding: utf-8 -*-
# Przyklad: kwadrokopter z ladunkiem podwieszonym na linie.
#
# Dron to cialo sztywne z czterema ciagami wirnikow (SilaWPunkcie: sila
# w ukladzie ciala, zaczepiona na koncu ramienia). Ladunek wisi na linie
# (SilaWewnProst z tylko_rozciaganie=True: luzna lina nie pcha).
#
# Sterowanie jak w prawdziwym kontrolerze lotu: dyskretny regulator
# kaskadowy PD (pozycja -> zadane pochylenie i ciag; orientacja -> momenty
# -> rozdzial na 4 wirniki) taktowany co SEGMENT; miedzy segmentami
# podmieniane sa wektory ciagu wirnikow, a stan (q, dq) leci dalej.
#
# Misja: start z ziemi, wznoszenie, przelot z ladunkiem do punktu,
# zejscie do zawisu. Ladunek buja sie jak wahadlo i szarpie dronem.
#
# Wynik: web/dane_dron.js do wizualizacji Three.js (web/dron.html).

import json
import os

import numpy as np

from uw_dyn import (Uklad, Czlon, SilaWPunkcie, SilaWewnProst, SilaKontaktu,
                    wektor, G, wektor_p)
from uw_dyn.uklad import GRAWITACJA

# ----- dron -----
M_DRON = 1.5          # [kg]
J_DRON = np.diag([0.02, 0.02, 0.035])
RAMIE = 0.18          # odleglosc wirnika od srodka [m]
T_MAKS = 12.0         # maksymalny ciag jednego wirnika [N]

# ----- ladunek i lina -----
M_LADUNKU = 0.4
J_LADUNKU = np.diag([2e-4, 2e-4, 2e-4])
DL_LINY = 1.2
K_LINY = 2000.0
C_LINY = 20.0

# ----- regulator (dyskretny, 100 Hz) -----
SEGMENT = 0.01
DT = 0.001
KP_POZ, KD_POZ = 2.0, 2.4       # pozycja -> zadane przyspieszenie
KP_OR, KD_OR = 2.0, 0.30        # orientacja -> momenty
MAKS_POCHYLENIE = 0.45          # [rad]

# wirniki w ukladzie ciala: +x, -x, +y, -y
WIRNIKI = [wektor(RAMIE, 0, 0), wektor(-RAMIE, 0, 0),
           wektor(0, RAMIE, 0), wektor(0, -RAMIE, 0)]

# misja: punkty (x, y, z) i promien zaliczenia; przelot tam i z powrotem
# z ladunkiem na linie (wahadlo szarpie dronem w obu kierunkach)
MISJA = [(0.0, 0.0, 2.0), (4.0, 0.0, 2.2), (4.0, 0.0, 1.4),
         (0.0, 0.0, 2.0), (0.0, 0.0, 1.4)]
PROMIEN = 0.35
CZAS_MAKS = 24.0


def zbuduj():
    ukl = Uklad()
    ukl.dodajCzlon(Czlon(1, M_DRON, J_DRON))
    ukl.dodajCzlon(Czlon(2, M_LADUNKU, J_LADUNKU))

    ciagi = []
    for s_w in WIRNIKI:
        sila = SilaWPunkcie(1, s_w, wektor(0, 0, 0))
        ciagi.append(sila)
        ukl.dodajSileWewn(sila)

    ukl.dodajSileWewn(SilaWewnProst(1, 2, wektor(0, 0, -0.05), wektor(0, 0, 0),
                                    K_LINY, DL_LINY, C_LINY, 0,
                                    tylko_rozciaganie=True))
    # ladunek lezy na ziemi, dopoki lina go nie podniesie
    ukl.dodajSileWewn(SilaKontaktu(2, wektor(0, 0, 0), k=5000.0, c=100.0, mu=0.6))
    ukl.grawitacja = True
    return ukl, ciagi


def stan_poczatkowy():
    """Dron nisko nad ziemia, ladunek lezy pod nim (lina luzna)."""
    q = np.zeros(14)
    q[2] = 0.3                    # dron
    q[3:6] = [0.0, 0.0, 0.05]     # ladunek na ziemi pod dronem
    q[6] = 1.0                    # kwaternion drona
    q[10] = 1.0                   # kwaternion ladunku
    return np.concatenate((q, np.zeros(14)))


def katy_i_predkosc_katowa(q, dq):
    """Przyblizone katy przechylenia/pochylenia i predkosc katowa drona."""
    p = q[6:10]
    roll = 2 * np.arctan2(p[1], p[0])
    pitch = 2 * np.arctan2(p[2], p[0])
    dp = np.array(dq[6:10]).reshape(4, 1)
    om = (2 * G(wektor_p(*p)).dot(dp)).ravel()
    return roll, pitch, om


def regulator(q, dq, cel):
    """Kaskadowy PD: pozycja -> ciag i zadane katy -> momenty -> 4 wirniki."""
    poz = q[0:3]
    v = dq[0:3]
    m_calk = M_DRON + M_LADUNKU

    # petla pozycji: zadane przyspieszenie + kompensacja grawitacji
    a_zad = KP_POZ * (np.array(cel) - poz) - KD_POZ * v
    F_zad = m_calk * (a_zad + np.array([0, 0, GRAWITACJA]))
    F_zad[2] = max(F_zad[2], 0.3 * m_calk * GRAWITACJA)  # nie wylaczaj silnikow

    # zadane pochylenie z kierunku sily
    pitch_zad = np.clip(np.arctan2(F_zad[0], F_zad[2]),
                        -MAKS_POCHYLENIE, MAKS_POCHYLENIE)
    roll_zad = np.clip(np.arctan2(-F_zad[1], F_zad[2]),
                       -MAKS_POCHYLENIE, MAKS_POCHYLENIE)

    roll, pitch, om = katy_i_predkosc_katowa(q, dq)
    tau_x = KP_OR * (roll_zad - roll) - KD_OR * om[0]
    tau_y = KP_OR * (pitch_zad - pitch) - KD_OR * om[1]

    T = float(np.linalg.norm(F_zad))
    # rozdzial: wirniki +x/-x steruja tau_y, +y/-y steruja tau_x
    T1 = T / 4 - tau_y / (2 * RAMIE)
    T2 = T / 4 + tau_y / (2 * RAMIE)
    T3 = T / 4 + tau_x / (2 * RAMIE)
    T4 = T / 4 - tau_x / (2 * RAMIE)
    return [float(np.clip(Ti, 0.0, T_MAKS)) for Ti in (T1, T2, T3, T4)]


def symuluj():
    ukl, ciagi = zbuduj()
    y = stan_poczatkowy()
    q, dq = y[0:14].copy(), y[14:28].copy()

    klatki = []
    nr_celu = 0
    t = 0.0
    while t < CZAS_MAKS and nr_celu < len(MISJA):
        cel = MISJA[nr_celu]
        T = regulator(q, dq, cel)
        for sila, Ti in zip(ciagi, T):
            sila.f_lokalna = wektor(0, 0, Ti)

        ukl.sym2(np.concatenate((q, dq)), 0.0, SEGMENT, DT)
        Y = ukl.Y
        for w in Y[:-1]:
            klatki.append({'dron_r': [round(v, 4) for v in w[0:3]],
                           'dron_p': [round(v, 5) for v in w[6:10]],
                           'lad_r': [round(v, 4) for v in w[3:6]],
                           'cel': nr_celu})
        q = Y[-1][0:14].copy()
        dq = Y[-1][14:28].copy()
        t += SEGMENT

        blad = np.linalg.norm(q[0:3] - np.array(cel))
        if blad < PROMIEN and np.linalg.norm(dq[0:3]) < 0.6:
            nr_celu += 1

    return klatki, nr_celu, t, q, dq


def eksportuj(klatki, co_ile=20, plik='web/dane_dron.js'):
    dane = {
        'dt': DT * co_ile,
        'ramie': RAMIE,
        'dl_liny': DL_LINY,
        'misja': MISJA,
        'klatki': klatki[::co_ile],
    }
    katalog = os.path.dirname(os.path.abspath(plik))
    os.makedirs(katalog, exist_ok=True)
    with open(plik, 'w') as f:
        f.write('const DANE = ' + json.dumps(dane) + ';\n')
    print(f"zapisano {len(dane['klatki'])} klatek do {plik}")


if __name__ == '__main__':
    klatki, nr_celu, t, q, dq = symuluj()
    print(f'zaliczone punkty misji: {nr_celu}/{len(MISJA)} w {t:.1f} s')
    print(f'pozycja koncowa drona: {q[0:3].round(3)}, ladunku: {q[3:6].round(3)}')

    # wahanie ladunku: kat liny od pionu (sprzezenie wieloczlonowe dron-lina-ladunek)
    katy = []
    for kl in klatki:
        d = np.array(kl['lad_r']) - np.array(kl['dron_r'])
        n = np.linalg.norm(d)
        if n > 0.5*DL_LINY:     # lina napieta
            katy.append(np.degrees(np.arccos(np.clip(-d[2]/n, -1, 1))))
    print(f'maks. wychylenie liny od pionu: {max(katy):.1f} st. '
          f'(koncowe: {katy[-1]:.1f} st.)')
    if nr_celu < len(MISJA):
        raise SystemExit('misja niezaliczona: dostroic regulator')

    sciezka = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           '..', 'web', 'dane_dron.js')
    eksportuj(klatki, plik=os.path.normpath(sciezka))
