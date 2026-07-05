# -*- coding: utf-8 -*-
# Przyklad: mini-piesek (czworonog) robiacy przysiady.
#
# Rozwiniecie przykladu przysiadu na cztery nogi: tulow swobodny
# (floating base) podparty na czterech dwuczlonowych nogach; stopy stoja
# na podlozu (SilaKontaktu), a stawy (biodro i kolano kazdej nogi)
# napedzane sa aktuatorami MomentWzgledny. Piesek symetrycznie ugina
# wszystkie nogi, tulow opada, po czym prostuje je i wstaje.
#
# Kinematyka przysiadu (nogi o rownych czlonach L): przy kacie uda phi
# i kacie kolana -2*phi stopa jest dokladnie pod biodrem, a tulow na
# wysokosci 2*L*cos(phi). Zmieniajac phi w czasie, tulow opada i wstaje,
# a stopy pozostaja w miejscu (brak poslizgu).
#
# Wynik: web/dane_piesek.js do wizualizacji Three.js (web/piesek.html).

import json
import os

import numpy as np

from uw_dyn import (Uklad, Czlon, Polaczenie_Obr, SilaKontaktu, MomentWzgledny,
                    wektor, u2p, R, wektor_p)

# ----- wymiary (metry, kilogramy) -----
DL_TULOW, SZER_TULOW = 0.40, 0.18
M_TULOW = 4.0
L_NOGI = 0.15         # dlugosc czlonu nogi (udo = podudzie) [m]
M_NOGI = 0.4
PROMIEN = 0.03

# rozmieszczenie bioder w ukladzie tulowia (przod +x, lewo +y)
BIODRA = {
    'PL': wektor(DL_TULOW/2, SZER_TULOW/2, 0),
    'PP': wektor(DL_TULOW/2, -SZER_TULOW/2, 0),
    'TL': wektor(-DL_TULOW/2, SZER_TULOW/2, 0),
    'TP': wektor(-DL_TULOW/2, -SZER_TULOW/2, 0),
}
NOGI = list(BIODRA.keys())

# ----- konfiguracja przysiadu -----
PHI_STAC = 0.42       # kat uda na stojaco [rad] -> tulow wyzej
PHI_KUCA = 0.95       # kat uda w glebokim przysiadzie -> tulow nizej
OS_Y = np.array([0.0, 1.0, 0.0])
OSIE = (wektor(1, 0, 0), wektor(0, 0, 1), wektor(0, 1, 0))  # przegub obrotowy wokol y


def wys_tulowia(phi):
    return 2*L_NOGI*np.cos(phi)


H_STAC = wys_tulowia(PHI_STAC)

# ----- sterownik przysiadu -----
SEGMENT = 0.005
DT = 0.0005
OKRES = 2.2           # okres jednego przysiadu [s]
LICZBA_PRZYSIADOW = 3
K_STAWU, C_STAWU = 80.0, 3.0
CZAS_STAC = 1.0       # stabilizacja na stojaco przed przysiadami [s]
CZAS = CZAS_STAC + LICZBA_PRZYSIADOW*OKRES + 0.6


def tensor_preta(m, L):
    Jp = m*(3*PROMIEN**2 + L**2)/12
    Jo = m*PROMIEN**2/2
    return np.diag([Jp, Jp, Jo])


def tensor_tulowia():
    Jx = M_TULOW*(SZER_TULOW**2 + 0.1**2)/12
    Jy = M_TULOW*(DL_TULOW**2 + 0.1**2)/12
    Jz = M_TULOW*(DL_TULOW**2 + SZER_TULOW**2)/12
    return np.diag([Jx, Jy, Jz])


def numery():
    nr = {'tulow': 1}
    k = 2
    for noga in NOGI:
        nr[noga+'_udo'] = k
        nr[noga+'_podudzie'] = k+1
        k += 2
    return nr


NR = numery()
N_CZLONOW = 1 + 2*len(NOGI)


def katy_nogi(phi):
    """Kat biodra i kolana dla zadanego kata uda phi (stopa pod biodrem)."""
    return phi, -2*phi


def zbuduj():
    ukl = Uklad()
    ukl.dodajCzlon(Czlon(NR['tulow'], M_TULOW, tensor_tulowia()))
    for noga in NOGI:
        ukl.dodajCzlon(Czlon(NR[noga+'_udo'], M_NOGI, tensor_preta(M_NOGI, L_NOGI)))
        ukl.dodajCzlon(Czlon(NR[noga+'_podudzie'], M_NOGI, tensor_preta(M_NOGI, L_NOGI)))

    kat_b0, kat_k0 = katy_nogi(PHI_STAC)
    aktuatory = {}
    for noga in NOGI:
        u, p = NR[noga+'_udo'], NR[noga+'_podudzie']
        ukl.dodajWiez(Polaczenie_Obr(NR['tulow'], u,
                                     BIODRA[noga], wektor(0, 0, L_NOGI/2), *OSIE))
        ukl.dodajWiez(Polaczenie_Obr(u, p,
                                     wektor(0, 0, -L_NOGI/2), wektor(0, 0, L_NOGI/2), *OSIE))
        ukl.dodajSileWewn(SilaKontaktu(p, wektor(0, 0, -L_NOGI/2),
                                       k=8000.0, c=80.0, mu=1.2, eps=0.002))
        akt_b = MomentWzgledny(NR['tulow'], u, wektor(0, 1, 0), wektor(0, 0, 1),
                               K_STAWU, kat_b0, C_STAWU)
        akt_k = MomentWzgledny(u, p, wektor(0, 1, 0), wektor(0, 0, 1),
                               K_STAWU, kat_k0, C_STAWU)
        ukl.dodajSileWewn(akt_b)
        ukl.dodajSileWewn(akt_k)
        aktuatory[noga] = (akt_b, akt_k)

    ukl.grawitacja = True
    return ukl, aktuatory


def q_startowe(phi):
    """Konfiguracja stojaca dla kata uda phi (stopy pod biodrami)."""
    kat_b, kat_k = katy_nogi(phi)
    h = wys_tulowia(phi)
    q = np.zeros(7*N_CZLONOW)
    q[0:3] = [0, 0, h]
    q[3*N_CZLONOW:3*N_CZLONOW + 4] = [1, 0, 0, 0]  # kwaternion tulowia

    for noga in NOGI:
        u, p = NR[noga+'_udo'], NR[noga+'_podudzie']
        biodro = wektor(0, 0, h) + BIODRA[noga]
        Ru = R(wektor_p(*u2p(OS_Y, kat_b)))
        srodek_uda = biodro + Ru.dot(wektor(0, 0, -L_NOGI/2))
        kolano = biodro + Ru.dot(wektor(0, 0, -L_NOGI))
        Rp = R(wektor_p(*u2p(OS_Y, kat_b + kat_k)))
        srodek_pod = kolano + Rp.dot(wektor(0, 0, -L_NOGI/2))

        q[3*(u-1):3*(u-1)+3] = srodek_uda.ravel()
        q[3*(p-1):3*(p-1)+3] = srodek_pod.ravel()
        q[3*N_CZLONOW + 4*(u-1):3*N_CZLONOW + 4*(u-1)+4] = u2p(OS_Y, kat_b)
        q[3*N_CZLONOW + 4*(p-1):3*N_CZLONOW + 4*(p-1)+4] = u2p(OS_Y, kat_b + kat_k)
    return q


def phi_celu(t):
    """Kat uda w chwili t: stanie, potem plynne przysiady (kosinusoida)."""
    if t < CZAS_STAC:
        return PHI_STAC
    faza = (t - CZAS_STAC)/OKRES
    if faza >= LICZBA_PRZYSIADOW:
        return PHI_STAC
    # 0 na gorze, 1 w dole, plynnie
    zejscie = (1 - np.cos(2*np.pi*faza))/2
    return PHI_STAC + (PHI_KUCA - PHI_STAC)*zejscie


def symuluj():
    ukl, aktuatory = zbuduj()
    q = ukl.projekcja_polozen(q_startowe(PHI_STAC))
    dq = np.zeros(7*N_CZLONOW)

    klatki = []
    t = 0.0
    n_seg = int(CZAS/SEGMENT)
    for seg in range(n_seg):
        kat_b, kat_k = katy_nogi(phi_celu(t))
        for noga in NOGI:
            akt_b, akt_k = aktuatory[noga]
            akt_b.theta_cel = kat_b
            akt_k.theta_cel = kat_k

        ukl.sym2(np.concatenate((q, dq)), 0.0, SEGMENT, DT)
        Y = ukl.Y
        for w in Y[:-1]:
            klatki.append(w[0:7*N_CZLONOW].copy())
        q = Y[-1][0:7*N_CZLONOW].copy()
        dq = Y[-1][7*N_CZLONOW:14*N_CZLONOW].copy()
        t += SEGMENT

    return ukl, klatki


def eksportuj(klatki, co_ile=12, plik='web/dane_piesek.js'):
    dane_klatki = []
    for q in klatki[::co_ile]:
        czlony = []
        for k in range(N_CZLONOW):
            r = q[3*k:3*k+3]
            p = q[3*N_CZLONOW + 4*k:3*N_CZLONOW + 4*k + 4]
            czlony.append({'r': [round(float(v), 4) for v in r],
                           'p': [round(float(v), 5) for v in p]})
        dane_klatki.append(czlony)

    dane = {
        'dt': DT*co_ile,
        'wymiary': {'udo': L_NOGI, 'podudzie': L_NOGI,
                    'tulow': DL_TULOW, 'szer': SZER_TULOW},
        'nogi': NOGI,
        'indeksy': {'tulow': 0,
                    **{noga: [NR[noga+'_udo']-1, NR[noga+'_podudzie']-1] for noga in NOGI}},
        'klatki': dane_klatki,
    }
    katalog = os.path.dirname(os.path.abspath(plik))
    os.makedirs(katalog, exist_ok=True)
    with open(plik, 'w') as f:
        f.write('const DANE = ' + json.dumps(dane) + ';\n')
    print(f'zapisano {len(dane_klatki)} klatek do {plik}')


if __name__ == '__main__':
    ukl, klatki = symuluj()
    zt = [q[2] for q in klatki]
    print(f'wysokosc tulowia: stojac {H_STAC:.3f} m, '
          f'zakres w przysiadach {min(zt):.3f}..{max(zt):.3f} m')
    print(f'najglebszy przysiad: tulow o {H_STAC - min(zt):.3f} m nizej')

    sciezka = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           '..', 'web', 'dane_piesek.js')
    eksportuj(klatki, plik=os.path.normpath(sciezka))
