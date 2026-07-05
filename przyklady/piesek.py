# -*- coding: utf-8 -*-
# Przyklad: mini-piesek (czworonog) w trucie.
#
# Rozwiniecie najprostszego robota kroczacego: zamiast dwoch patykow
# pelny czworonog z podstawa swobodna (floating base):
#   - tulow (1 czlon),
#   - 4 nogi po 2 czlony (udo + podudzie), przeguby obrotowe (os y) w
#     biodrze i kolanie,
#   - stopy dotykaja podloza przez SilaKontaktu (model penalty z tarciem),
#   - stawy napedzane aktuatorami MomentWzgledny (sprezyna-tlumik z celem
#     katowym), a cele katow modulowane sa faza chodu (dyskretny sterownik
#     co SEGMENT).
#
# Chod: pelzajacy (crawl) - w danej chwili tylko jedna noga jest w wymachu,
# trzy pozostale podpieraja (chod statycznie stabilny). W fazie podporowej
# biodro przesuwa sie do tylu, odpychajac tulow do przodu; w fazie wymachu
# kolano zgina sie (unos stopy), a biodro wraca do przodu.
#
# Wynik: web/dane_piesek.js do wizualizacji Three.js (web/piesek.html).

import json
import os

import numpy as np

from uw_dyn import (Uklad, Czlon, Polaczenie_Obr, SilaKontaktu, MomentWzgledny,
                    wektor, u2p, R, wektor_p)
from uw_dyn.uklad import GRAWITACJA

# ----- wymiary (metry, kilogramy) -----
DL_TULOW, SZER_TULOW = 0.40, 0.18
M_TULOW = 4.0
L_UDO, L_PODUDZIE = 0.14, 0.14
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
# przesuniecia faz chodu pelzajacego (sekwencja: TL, PL, TP, PP),
# co cwierc okresu inna noga w wymachu -> zawsze 3 stopy na ziemi
PRZESUNIECIE = {'TL': 0.0, 'PL': 0.5, 'TP': 0.25, 'PP': 0.75}

# geometria nogi: udo i podudzie skierowane w dol, lekkie zgiecie kolana
KAT_BIODRA0 = 0.15    # wychylenie uda do przodu [rad]
KAT_KOLANA0 = -0.55   # zgiecie kolana (podudzie do tylu) [rad]
H_TULOW = L_UDO*np.cos(KAT_BIODRA0) + L_PODUDZIE*np.cos(KAT_BIODRA0+KAT_KOLANA0)

OS_Y = np.array([0.0, 1.0, 0.0])
OSIE = (wektor(1, 0, 0), wektor(0, 0, 1), wektor(0, 1, 0))  # przegub obrotowy wokol y

# ----- sterownik chodu -----
SEGMENT = 0.005
DT = 0.001
OKRES = 1.2           # okres pelnego cyklu chodu [s]
DUTY = 0.75           # ulamek cyklu w fazie podporowej (3 stopy na ziemi)
AMPL_BIODRA = 0.28    # zakres przesuwu biodra przod-tyl [rad]
UNOS_KOLANA = 0.7     # dodatkowe zgiecie kolana w wymachu (unos stopy) [rad]
K_BIODRA, C_BIODRA = 80.0, 2.5
K_KOLANA, C_KOLANA = 80.0, 2.5
CZAS = 6.0


def tensor_preta(m, L):
    Jp = m*(3*PROMIEN**2 + L**2)/12
    Jo = m*PROMIEN**2/2
    return np.diag([Jp, Jp, Jo])


def tensor_tulowia():
    Jx = M_TULOW*(SZER_TULOW**2 + (0.1)**2)/12
    Jy = M_TULOW*(DL_TULOW**2 + (0.1)**2)/12
    Jz = M_TULOW*(DL_TULOW**2 + SZER_TULOW**2)/12
    return np.diag([Jx, Jy, Jz])


def numery():
    """Mapa nazw czlonow na numery: tulow=1, potem udo/podudzie kazdej nogi."""
    nr = {'tulow': 1}
    k = 2
    for noga in NOGI:
        nr[noga+'_udo'] = k
        nr[noga+'_podudzie'] = k+1
        k += 2
    return nr


NR = numery()
N_CZLONOW = 1 + 2*len(NOGI)


def zbuduj():
    ukl = Uklad()
    ukl.dodajCzlon(Czlon(NR['tulow'], M_TULOW, tensor_tulowia()))
    for noga in NOGI:
        ukl.dodajCzlon(Czlon(NR[noga+'_udo'], M_NOGI, tensor_preta(M_NOGI, L_UDO)))
        ukl.dodajCzlon(Czlon(NR[noga+'_podudzie'], M_NOGI, tensor_preta(M_NOGI, L_PODUDZIE)))

    aktuatory = {}
    for noga in NOGI:
        u, p = NR[noga+'_udo'], NR[noga+'_podudzie']
        # biodro: tulow (punkt biodra) - gora uda
        ukl.dodajWiez(Polaczenie_Obr(NR['tulow'], u,
                                     BIODRA[noga], wektor(0, 0, L_UDO/2), *OSIE))
        # kolano: dol uda - gora podudzia
        ukl.dodajWiez(Polaczenie_Obr(u, p,
                                     wektor(0, 0, -L_UDO/2), wektor(0, 0, L_PODUDZIE/2), *OSIE))
        # kontakt stopy (dol podudzia); male eps -> sztywne tarcie, malo slizgu
        ukl.dodajSileWewn(SilaKontaktu(p, wektor(0, 0, -L_PODUDZIE/2),
                                       k=8000.0, c=80.0, mu=1.2, eps=0.002))
        # aktuatory: biodro (tulow-udo), kolano (udo-podudzie)
        akt_b = MomentWzgledny(NR['tulow'], u, wektor(0, 1, 0), wektor(0, 0, 1),
                               K_BIODRA, KAT_BIODRA0, C_BIODRA)
        akt_k = MomentWzgledny(u, p, wektor(0, 1, 0), wektor(0, 0, 1),
                               K_KOLANA, KAT_KOLANA0, C_KOLANA)
        ukl.dodajSileWewn(akt_b)
        ukl.dodajSileWewn(akt_k)
        aktuatory[noga] = (akt_b, akt_k)

    ukl.grawitacja = True
    return ukl, aktuatory


def q_startowe(h_tulow):
    """Konfiguracja stojaca: tulow na wysokosci h, nogi zgiete symetrycznie."""
    q = np.zeros(7*N_CZLONOW)
    q[0:3] = [0, 0, h_tulow]
    q[3*N_CZLONOW + 0:3*N_CZLONOW + 4] = [1, 0, 0, 0]  # kwaternion tulowia
    Rt = np.eye(3)

    for noga in NOGI:
        u, p = NR[noga+'_udo'], NR[noga+'_podudzie']
        biodro_glob = wektor(0, 0, h_tulow) + Rt.dot(BIODRA[noga])

        Ru = R(wektor_p(*u2p(OS_Y, KAT_BIODRA0)))
        srodek_uda = biodro_glob + Ru.dot(wektor(0, 0, -L_UDO/2))
        kolano_glob = biodro_glob + Ru.dot(wektor(0, 0, -L_UDO))

        Rp = R(wektor_p(*u2p(OS_Y, KAT_BIODRA0 + KAT_KOLANA0)))
        srodek_pod = kolano_glob + Rp.dot(wektor(0, 0, -L_PODUDZIE/2))

        q[3*(u-1):3*(u-1)+3] = srodek_uda.ravel()
        q[3*(p-1):3*(p-1)+3] = srodek_pod.ravel()
        q[3*N_CZLONOW + 4*(u-1):3*N_CZLONOW + 4*(u-1)+4] = u2p(OS_Y, KAT_BIODRA0)
        q[3*N_CZLONOW + 4*(p-1):3*N_CZLONOW + 4*(p-1)+4] = u2p(OS_Y, KAT_BIODRA0 + KAT_KOLANA0)
    return q


def cele_chodu(t):
    """Zadane katy biodra i kolana dla kazdej nogi w chwili t (chod pelzajacy).

    Faza [0, DUTY): podpora, biodro przesuwa sie do tylu (odpycha tulow).
    Faza [DUTY, 1): wymach, biodro wraca do przodu, kolano zgina sie (unos)."""
    cele = {}
    for noga in NOGI:
        faza = (t/OKRES + PRZESUNIECIE[noga]) % 1.0
        if faza < DUTY:
            # podpora: biodro od +ampl/2 (przod) do -ampl/2 (tyl) liniowo
            s = faza/DUTY
            kat_b = KAT_BIODRA0 + AMPL_BIODRA*(0.5 - s)
            kat_k = KAT_KOLANA0
        else:
            # wymach: biodro wraca do przodu, kolano sie zgina i prostuje
            s = (faza - DUTY)/(1.0 - DUTY)
            kat_b = KAT_BIODRA0 + AMPL_BIODRA*(-0.5 + s)
            kat_k = KAT_KOLANA0 - UNOS_KOLANA*np.sin(np.pi*s)
        cele[noga] = (kat_b, kat_k)
    return cele


def symuluj():
    ukl, aktuatory = zbuduj()
    q = q_startowe(H_TULOW)
    dq = np.zeros(7*N_CZLONOW)

    # osadzenie na wiezach i chwila stabilizacji stojac (bez nagrywania:
    # poczatkowy transient osiadania na kontakcie nie jest interesujacy)
    q = ukl.projekcja_polozen(q)
    klatki = []
    t = 0.0
    n_seg = int(CZAS/SEGMENT)
    start_chodu = 1.0

    for seg in range(n_seg):
        if t >= start_chodu:
            cele = cele_chodu(t - start_chodu)
            for noga in NOGI:
                akt_b, akt_k = aktuatory[noga]
                akt_b.theta_cel, akt_k.theta_cel = cele[noga]

        ukl.sym2(np.concatenate((q, dq)), 0.0, SEGMENT, DT)
        Y = ukl.Y
        if t >= start_chodu:
            for w in Y[:-1]:
                klatki.append(w[0:7*N_CZLONOW].copy())
        q = Y[-1][0:7*N_CZLONOW].copy()
        dq = Y[-1][7*N_CZLONOW:14*N_CZLONOW].copy()
        t += SEGMENT

    return ukl, klatki, q


def eksportuj(klatki, co_ile=15, plik='web/dane_piesek.js'):
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
        'wymiary': {'udo': L_UDO, 'podudzie': L_PODUDZIE,
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
    ukl, klatki, q = symuluj()
    x_koncowe = float(q[0])
    z_koncowe = float(q[2])
    print(f'tulow: x={x_koncowe:.3f} m (przebyta droga), wysokosc z={z_koncowe:.3f} m')
    if z_koncowe < 0.5*H_TULOW:
        print('UWAGA: piesek sie przewrocil (niska wysokosc tulowia)')

    sciezka = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           '..', 'web', 'dane_piesek.js')
    eksportuj(klatki, plik=os.path.normpath(sciezka))
