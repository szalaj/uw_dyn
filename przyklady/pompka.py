# -*- coding: utf-8 -*-
# Przyklad: zawodnik robi pompki (push-up).
#
# Model w plaszczyznie strzalkowej (x-z):
#   - TULOW: sztywna deska (tulow + nogi + glowa) obrotowo zamocowana w palcach
#     stop (przegub obrotowy do podstawy, os y) - cialo obraca sie wokol stop
#     jak dzwignia,
#   - dwa RAMIONA (bark + przedramie, zawiasy wokol osi y): bark laczy deske
#     z ramieniem, lokiec ramie z przedramieniem,
#   - DLONIE na kontakcie z podlozem (SilaKontaktu): reakcja podloza + tarcie
#     trzyma dlonie w miejscu, a ramiona pchaja o nie cialo w gore.
#
# Sterowanie: aktuatory PID. Barki trzymaja ramie ustawione ku podlozu, a
# lokcie sa napedzane sinusoidalnie miedzy wyprostem (gora - cialo uniesione)
# a zgieciem (dol - klatka nisko). Wyprost lokci pcha o dlonie -> tulow rotuje
# w gore wokol stop = pompka.
#
# Calkowanie: sym3 polniejawny (RATTLE + niejawne tlumienie/sprezyny), dt=1e-3
# - sztywny kontakt dloni bez maleńkiego kroku.
#
# Wynik: web/dane_pompka.js do wizualizacji Three.js (web/pompka.html).

import json
import os

import numpy as np

from uw_dyn import (Uklad, Czlon, Polaczenie_Obr, MomentWzgledny, SilaKontaktu,
                    wektor, u2p, R, wektor_p, mnoz_kwaterniony)

# ----- wymiary (metry, kilogramy) -----
L_TUL = 1.35              # dlugosc deski (palce stop -> barki)
M_TUL = 58.0             # masa tulowia+nog+glowy [kg]
R_TUL = 0.13
L_UA, L_FA = 0.30, 0.28  # ramie (bark-lokiec), przedramie (lokiec-dlon)
M_UA, M_FA = 2.6, 1.6
R_ARM = 0.05
Y_BARK = 0.20            # polowa rozstawu barkow (i dloni)
Z_TOE = 0.06             # wysokosc palcow stop (przegub) [m]

TUL, UAL, FAL, UAP, FAP = 1, 2, 3, 4, 5
N = 5
OS_Y = np.array([0.0, 1.0, 0.0])
# zawias wokol osi y: (vi, wi, uj) tak, ze vi x wi = y
OSIE_Y = (wektor(0, 0, 1), wektor(1, 0, 0), wektor(0, 1, 0))

# ----- sterownik -----
SEGMENT = 0.01
DT = 1.0e-3
K_SH, C_SH, KI_SH, CMAX_SH = 220.0, 20.0, 120.0, 200.0   # bark (trzyma ramie)
K_ELB, C_ELB, KI_ELB, CMAX_ELB = 200.0, 18.0, 150.0, 200.0  # lokiec (naped)
LICZBA_POMPEK = 3
OKRES = 1.4              # czas jednej pompki [s]
CZAS = LICZBA_POMPEK*OKRES + 0.4
# katy lokcia (wzgledny, 0 = wyprost): gora ~wyprost, dol ~zgiecie
ELB_GORA = 0.15
ELB_DOL = 1.15


def tensor_preta(m, L, r):
    Jo = m*r**2/2
    Jp = m*(3*r**2 + L**2)/12
    return np.diag([Jp, Jp, Jo])


def _os_z(p):
    return R(wektor_p(*p))[:, 2]


def q_startowe(theta, elb):
    """FK pozy pompki: theta = nachylenie deski od poziomu, elb = kat lokcia.
    Zwraca (q, dane) gdzie dane zawiera pozycje pomocnicze."""
    q = np.zeros(7*N)
    b = 3*N
    # deska: palce w (0,0,Z_TOE), lokalna z wzdluz ciala (palce -> glowa)
    p_toe = np.array([0.0, 0.0, Z_TOE])
    # deska pozioma: lokalna z (palce->glowa) nachylona o theta nad poziom
    p_tul = u2p(OS_Y, np.pi/2 - theta)
    Rt = R(wektor_p(*p_tul))
    os_ciala = Rt[:, 2]                        # lokalna z deski (swiat)
    srodek_tul = p_toe + os_ciala*(L_TUL/2)
    q[3*(TUL-1):3*(TUL-1)+3] = srodek_tul
    q[b + 4*(TUL-1):b + 4*(TUL-1)+4] = p_tul

    # barki: przy koncu deski (+z), rozstaw +-Y_BARK
    for strona, (ua, fa) in (('L', (UAL, FAL)), ('P', (UAP, FAP))):
        y = Y_BARK if strona == 'L' else -Y_BARK
        bark = srodek_tul + Rt.dot(np.array([0, y, L_TUL/2 - 0.05]))
        # ramie i przedramie pionowo w dol; lokiec zgiety o elb (do tylu)
        p_ua = u2p(OS_Y, np.pi)               # lokalna z w dol
        p_fa = mnoz_kwaterniony(p_ua, u2p(OS_Y, elb))
        srodek_ua = bark + _os_z(p_ua)*(L_UA/2)
        lokiec = bark + _os_z(p_ua)*L_UA
        srodek_fa = lokiec + _os_z(p_fa)*(L_FA/2)
        q[3*(ua-1):3*(ua-1)+3] = srodek_ua
        q[3*(fa-1):3*(fa-1)+3] = srodek_fa
        q[b + 4*(ua-1):b + 4*(ua-1)+4] = p_ua
        q[b + 4*(fa-1):b + 4*(fa-1)+4] = p_fa
    return q


def _dlon(q, fa):
    r = q[3*(fa-1):3*(fa-1)+3]
    p = q[3*N + 4*(fa-1):3*N + 4*(fa-1)+4]
    return r + _os_z(p)*(L_FA/2)


def zbuduj(theta0, elb0):
    ukl = Uklad()
    ukl.dodajCzlon(Czlon(TUL, M_TUL, tensor_preta(M_TUL, L_TUL, R_TUL)))
    for nr in (UAL, FAL, UAP, FAP):
        m, L = (M_UA, L_UA) if nr in (UAL, UAP) else (M_FA, L_FA)
        ukl.dodajCzlon(Czlon(nr, m, tensor_preta(m, L, R_ARM)))

    # przegub palcow stop: deska obraca sie wokol y w punkcie (0,0,Z_TOE)
    ukl.dodajWiez(Polaczenie_Obr(0, TUL, wektor(0, 0, Z_TOE),
                                 wektor(0, 0, -L_TUL/2), *OSIE_Y))
    akt = {}
    for strona, (ua, fa) in (('L', (UAL, FAL)), ('P', (UAP, FAP))):
        y = Y_BARK if strona == 'L' else -Y_BARK
        ukl.dodajWiez(Polaczenie_Obr(TUL, ua, wektor(0, y, L_TUL/2 - 0.05),
                                     wektor(0, 0, -L_UA/2), *OSIE_Y))
        akt['bark_'+strona] = MomentWzgledny(TUL, ua, wektor(0, 1, 0),
                                             wektor(0, 0, 1), K_SH, 0.0, C_SH,
                                             ki=KI_SH, calka_max=CMAX_SH)
        ukl.dodajWiez(Polaczenie_Obr(ua, fa, wektor(0, 0, L_UA/2),
                                     wektor(0, 0, -L_FA/2), *OSIE_Y))
        akt['lokiec_'+strona] = MomentWzgledny(ua, fa, wektor(0, 1, 0),
                                               wektor(0, 0, 1), K_ELB, elb0,
                                               C_ELB, ki=KI_ELB,
                                               calka_max=CMAX_ELB)
    for a in akt.values():
        ukl.dodajSileWewn(a)
    # kontakt dloni z podlozem (koniec przedramienia)
    for fa in (FAL, FAP):
        ukl.dodajSileWewn(SilaKontaktu(fa, wektor(0, 0, L_FA/2),
                                       k=3.0e4, c=300.0, mu=1.2, eps=0.003))
    ukl.grawitacja = True

    # cel barku = kat w pozie startowej (punkt staly)
    q0 = q_startowe(theta0, elb0)
    for strona in ('L', 'P'):
        akt['bark_'+strona].theta_cel = akt['bark_'+strona].kat(q0, N)
    return ukl, akt, q0


def _theta_dla_wyprostu(elb):
    """Dobiera nachylenie deski tak, by przy danym kacie lokcia dlon byla na
    podlozu (z=0)."""
    lo, hi = 0.05, 0.9
    for _ in range(40):
        th = 0.5*(lo + hi)
        z = _dlon(q_startowe(th, elb), FAL)[2]
        if z > 0:
            hi = th
        else:
            lo = th
    return 0.5*(lo + hi)


def sterowanie(t, akt, bazy):
    """Naped lokci: sinusoida miedzy wyprostem (gora) a zgieciem (dol)."""
    faza = (1 - np.cos(2*np.pi*min(t, LICZBA_POMPEK*OKRES)/OKRES))/2  # 0->1->0
    elb = ELB_GORA + (ELB_DOL - ELB_GORA)*faza
    for strona in ('L', 'P'):
        akt['lokiec_'+strona].theta_cel = elb


def symuluj():
    theta0 = _theta_dla_wyprostu(ELB_GORA)
    ukl, akt, q0 = zbuduj(theta0, ELB_GORA)
    q = ukl.projekcja_polozen(q0)
    dq = np.zeros(7*N)

    klatki = []
    kat_min, kat_max = np.inf, -np.inf
    t = 0.0
    while t < CZAS:
        sterowanie(t, akt, None)
        ukl.sym3(np.concatenate((q, dq)), 0.0, SEGMENT, DT, polniejawne=True)
        Y = ukl.Y
        for w in Y[:-1]:
            klatki.append(w[0:7*N].copy())
        q = Y[-1][0:7*N].copy()
        dq = Y[-1][7*N:14*N].copy()
        # nachylenie deski (kat osi ciala od poziomu)
        os_c = _os_z(q[3*N + 4*(TUL-1):3*N + 4*(TUL-1)+4])
        kat = np.degrees(np.arctan2(os_c[2], os_c[0]))
        kat_min, kat_max = min(kat_min, kat), max(kat_max, kat)
        t += SEGMENT
    return ukl, klatki, kat_min, kat_max


def eksportuj(klatki, co_ile=5, plik='web/dane_pompka.js'):
    dane_klatki = []
    for q in klatki[::co_ile]:
        czlony = []
        for k in range(N):
            r = q[3*k:3*k+3]
            p = q[3*N + 4*k:3*N + 4*k + 4]
            czlony.append({'r': [round(float(v), 4) for v in r],
                           'p': [round(float(v), 5) for v in p]})
        dane_klatki.append(czlony)
    dane = {
        'dt': DT*co_ile,
        'wymiary': {'tul': L_TUL, 'r_tul': R_TUL, 'ramie': L_UA,
                    'przedramie': L_FA, 'r_arm': R_ARM},
        'indeksy': {'tul': TUL-1, 'ua_L': UAL-1, 'fa_L': FAL-1,
                    'ua_P': UAP-1, 'fa_P': FAP-1},
        'klatki': dane_klatki,
    }
    katalog = os.path.dirname(os.path.abspath(plik))
    os.makedirs(katalog, exist_ok=True)
    with open(plik, 'w') as f:
        f.write('const DANE = ' + json.dumps(dane) + ';\n')
    print(f'zapisano {len(dane_klatki)} klatek do {plik}')


if __name__ == '__main__':
    ukl, klatki, kat_min, kat_max = symuluj()
    print(f'pompki: {LICZBA_POMPEK}, nachylenie deski od {kat_min:.1f} do '
          f'{kat_max:.1f} st. (zakres {kat_max - kat_min:.1f} st.)')
    print('brak NaN:', not np.isnan(np.array(klatki)).any())
    sciezka = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           '..', 'web', 'dane_pompka.js')
    eksportuj(klatki, plik=os.path.normpath(sciezka))
