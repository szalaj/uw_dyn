# -*- coding: utf-8 -*-
# Przyklad: kickboxing - walka z cieniem. Prawy sierpowy.
#
# Bokser w pozycji wyjsciowej (garda) rzuca prawy sierpowy (hook): cios
# lukiem w plaszczyznie poziomej, napedzany obrotem tulowia (praca bioder
# i skretu korpusu). Model gornej czesci ciala:
#   - tulow (przypiety do podloza przegubem obrotowym wokol pionu w pasie:
#     skret korpusu to zrodlo mocy ciosu),
#   - prawe i lewe ramie: ramie (upper) + przedramie (fore), przeguby
#     obrotowe wokol pionu (bark, lokiec).
# Stawy napedzane aktuatorami MomentWzgledny (sprezyna-tlumik z celem
# katowym). Sekwencja katow: garda -> wyprowadzenie sierpowego -> powrot.
#
# Sierpowy to ruch poziomy, dlatego wszystkie osie obrotu sa pionowe (z),
# a czlony ramion rozciagaja sie wzdluz lokalnej osi x (grawitacja nie daje
# momentu wzgledem pionu, wiec ramiona trzyma sam wiez, bez oklapniecia).
#
# Nogi i glowa sa rysowane pogladowo w wizualizacji (nieruchoma postawa).
#
# Wynik: web/dane_bokser.js do wizualizacji Three.js (web/bokser.html).

import json
import os

import numpy as np

from uw_dyn import Uklad, Czlon, Polaczenie_Obr, MomentWzgledny, wektor, u2p

# ----- wymiary (metry, kilogramy) -----
TUL_DOL, TUL_GORA = 0.95, 1.50   # pas i barki (wysokosc) [m]
M_TUL = 34.0
Z_BARK = 1.42                    # wysokosc barkow [m]
ROZSTAW_BARK = 0.19             # polowa rozstawu barkow [m]
L_RAMIE, L_PRZEDRAMIE = 0.30, 0.28
M_RAMIE, M_PRZEDRAMIE = 2.6, 1.6
PROMIEN = 0.05

TUL_H = TUL_GORA - TUL_DOL
Z_TUL = (TUL_GORA + TUL_DOL)/2   # srodek masy tulowia
BARK_W_TUL = Z_BARK - Z_TUL      # bark w ukladzie tulowia (skladowa z)

OS_Z = np.array([0.0, 0.0, 1.0])
# przegub obrotowy wokol pionu z: os uj=z na czlonie j, prostopadla do x,y na i
OSIE_Z = (wektor(1, 0, 0), wektor(0, 1, 0), wektor(0, 0, 1))
REF_X = wektor(1, 0, 0)          # wektor odniesienia aktuatora (w plaszczyznie obrotu)

# numery czlonow
TUL, RA_G, RA_P, LA_G, LA_P = 1, 2, 3, 4, 5
N_CZLONOW = 5

# ----- katy (yaw wokol z), mierzone od +x (przod) ku +y (lewo) -----
# garda: postawa lekko zbladowana, prawa (tylna) piesc przy brodzie
YAW_TUL_G = 0.30      # skret tulowia (prawy bark cofniety)
RA_G_REL_G = -0.63    # bark: ramie wzgledem tulowia
RA_P_REL_G = 2.60     # lokiec mocno zgiety (piesc schowana przy twarzy)
# lewe ramie (przednia reka): lustrzana garda
LA_G_REL = 0.63
LA_P_REL = -2.60

# prawy sierpowy: mocny skret tulowia + zamach barku, lokiec sie otwiera
YAW_TUL_H = -0.45
RA_G_REL_H = 0.39
RA_P_REL_H = 0.90

# ----- sterownik (czasy faz) -----
SEGMENT = 0.004
DT = 0.0005
K_STAWU, C_STAWU = 220.0, 5.0     # bark i lokiec (snappy)
K_TUL, C_TUL = 400.0, 16.0        # skret korpusu (ciezki, ale mocny)
T_GARDA = 0.6                     # postawa wyjsciowa
T_ZAMACH = 0.14                   # wyprowadzenie ciosu (szybkie)
T_TRAFIENIE = 0.10               # pelne wyprostowanie
T_POWROT = 0.34                  # sciagniecie do gardy
T_KONIEC = 0.7
CZAS = T_GARDA + T_ZAMACH + T_TRAFIENIE + T_POWROT + T_KONIEC


def tensor_ramienia(m, L):
    Jo = m*PROMIEN**2/2                 # wzdluz wlasnej osi (x)
    Jp = m*(3*PROMIEN**2 + L**2)/12     # poprzecznie
    return np.diag([Jo, Jp, Jp])


def tensor_tulowia():
    Jx = M_TUL*(0.32**2 + TUL_H**2)/12
    Jy = M_TUL*(0.22**2 + TUL_H**2)/12
    Jz = M_TUL*(0.32**2 + 0.22**2)/12
    return np.diag([Jx, Jy, Jz])


def zbuduj():
    ukl = Uklad()
    ukl.dodajCzlon(Czlon(TUL, M_TUL, tensor_tulowia()))
    ukl.dodajCzlon(Czlon(RA_G, M_RAMIE, tensor_ramienia(M_RAMIE, L_RAMIE)))
    ukl.dodajCzlon(Czlon(RA_P, M_PRZEDRAMIE, tensor_ramienia(M_PRZEDRAMIE, L_PRZEDRAMIE)))
    ukl.dodajCzlon(Czlon(LA_G, M_RAMIE, tensor_ramienia(M_RAMIE, L_RAMIE)))
    ukl.dodajCzlon(Czlon(LA_P, M_PRZEDRAMIE, tensor_ramienia(M_PRZEDRAMIE, L_PRZEDRAMIE)))

    # tulow przypiety do podloza przegubem obrotowym wokol pionu w pasie
    ukl.dodajWiez(Polaczenie_Obr(0, TUL, wektor(0, 0, TUL_DOL),
                                 wektor(0, 0, -TUL_H/2), *OSIE_Z))
    # barki (na tulowie) -> ramiona (koniec ramienia wzdluz lokalnej osi x)
    ukl.dodajWiez(Polaczenie_Obr(TUL, RA_G, wektor(0, -ROZSTAW_BARK, BARK_W_TUL),
                                 wektor(-L_RAMIE/2, 0, 0), *OSIE_Z))
    ukl.dodajWiez(Polaczenie_Obr(RA_G, RA_P, wektor(L_RAMIE/2, 0, 0),
                                 wektor(-L_PRZEDRAMIE/2, 0, 0), *OSIE_Z))
    ukl.dodajWiez(Polaczenie_Obr(TUL, LA_G, wektor(0, ROZSTAW_BARK, BARK_W_TUL),
                                 wektor(-L_RAMIE/2, 0, 0), *OSIE_Z))
    ukl.dodajWiez(Polaczenie_Obr(LA_G, LA_P, wektor(L_RAMIE/2, 0, 0),
                                 wektor(-L_PRZEDRAMIE/2, 0, 0), *OSIE_Z))

    # aktuatory obrotowe (wokol z) w kazdym przegubie
    akt = {}
    akt['tul'] = MomentWzgledny(0, TUL, wektor(0, 0, 1), REF_X, K_TUL, YAW_TUL_G, C_TUL)
    akt['ra_g'] = MomentWzgledny(TUL, RA_G, wektor(0, 0, 1), REF_X, K_STAWU, RA_G_REL_G, C_STAWU)
    akt['ra_p'] = MomentWzgledny(RA_G, RA_P, wektor(0, 0, 1), REF_X, K_STAWU, RA_P_REL_G, C_STAWU)
    akt['la_g'] = MomentWzgledny(TUL, LA_G, wektor(0, 0, 1), REF_X, K_STAWU, LA_G_REL, C_STAWU)
    akt['la_p'] = MomentWzgledny(LA_G, LA_P, wektor(0, 0, 1), REF_X, K_STAWU, LA_P_REL, C_STAWU)
    for a in akt.values():
        ukl.dodajSileWewn(a)

    ukl.grawitacja = True
    return ukl, akt


def _rz(yaw, v):
    c, s = np.cos(yaw), np.sin(yaw)
    x, y, z = v
    return np.array([c*x - s*y, s*x + c*y, z])


def q_startowe():
    """Poza gardy: pozycje srodkow mas i kwaterniony (obrot wokol z)."""
    yaw_tul = YAW_TUL_G
    yaw_ra_g = yaw_tul + RA_G_REL_G
    yaw_ra_p = yaw_ra_g + RA_P_REL_G
    yaw_la_g = yaw_tul + LA_G_REL
    yaw_la_p = yaw_la_g + LA_P_REL

    q = np.zeros(7*N_CZLONOW)

    # tulow
    q[3*(TUL-1):3*(TUL-1)+3] = [0, 0, Z_TUL]

    def ramie(bark_tul, yaw_g, yaw_p):
        bark = np.array([0, 0, Z_TUL]) + _rz(yaw_tul, bark_tul)
        srodek_g = bark + _rz(yaw_g, np.array([L_RAMIE/2, 0, 0]))
        lokiec = bark + _rz(yaw_g, np.array([L_RAMIE, 0, 0]))
        srodek_p = lokiec + _rz(yaw_p, np.array([L_PRZEDRAMIE/2, 0, 0]))
        return srodek_g, srodek_p

    sg, sp = ramie(np.array([0, -ROZSTAW_BARK, BARK_W_TUL]), yaw_ra_g, yaw_ra_p)
    q[3*(RA_G-1):3*(RA_G-1)+3] = sg
    q[3*(RA_P-1):3*(RA_P-1)+3] = sp
    lg, lp = ramie(np.array([0, ROZSTAW_BARK, BARK_W_TUL]), yaw_la_g, yaw_la_p)
    q[3*(LA_G-1):3*(LA_G-1)+3] = lg
    q[3*(LA_P-1):3*(LA_P-1)+3] = lp

    b = 3*N_CZLONOW
    for nr, yaw in ((TUL, yaw_tul), (RA_G, yaw_ra_g), (RA_P, yaw_ra_p),
                    (LA_G, yaw_la_g), (LA_P, yaw_la_p)):
        q[b + 4*(nr-1):b + 4*(nr-1)+4] = u2p(OS_Z, yaw)
    return q


def _ramp(t, t0, dt_faza):
    """Gladkie przejscie 0->1 na odcinku [t0, t0+dt_faza] (kosinusoida)."""
    if t <= t0:
        return 0.0
    if t >= t0 + dt_faza:
        return 1.0
    return (1 - np.cos(np.pi*(t - t0)/dt_faza))/2


def cele(t):
    """Cele katowe w chwili t: garda, wyprowadzenie, trafienie, powrot."""
    t1 = T_GARDA
    t2 = t1 + T_ZAMACH + T_TRAFIENIE
    # a: 0 = garda, 1 = pelny sierpowy
    if t < t1:
        a = 0.0
    elif t < t2:
        a = _ramp(t, t1, T_ZAMACH + T_TRAFIENIE)
    else:
        a = 1.0 - _ramp(t, t2, T_POWROT)

    return {
        'tul': YAW_TUL_G + (YAW_TUL_H - YAW_TUL_G)*a,
        'ra_g': RA_G_REL_G + (RA_G_REL_H - RA_G_REL_G)*a,
        'ra_p': RA_P_REL_G + (RA_P_REL_H - RA_P_REL_G)*a,
        'la_g': LA_G_REL,
        'la_p': LA_P_REL,
    }


def piesc_prawa(q):
    """Pozycja prawej piesci (koniec prawego przedramienia) w konfiguracji q."""
    from uw_dyn import R, wektor_p
    r = q[3*(RA_P-1):3*(RA_P-1)+3].reshape(3, 1)
    p = q[3*N_CZLONOW + 4*(RA_P-1):3*N_CZLONOW + 4*(RA_P-1)+4]
    return (r + R(wektor_p(*p)).dot(wektor(L_PRZEDRAMIE/2, 0, 0))).ravel()


def symuluj():
    ukl, akt = zbuduj()
    q = ukl.projekcja_polozen(q_startowe())
    dq = np.zeros(7*N_CZLONOW)

    klatki = []
    t = 0.0
    n_seg = int(CZAS/SEGMENT)
    for seg in range(n_seg):
        c = cele(t)
        for nazwa, wartosc in c.items():
            akt[nazwa].theta_cel = wartosc

        ukl.sym2(np.concatenate((q, dq)), 0.0, SEGMENT, DT)
        Y = ukl.Y
        for w in Y[:-1]:
            klatki.append(w[0:7*N_CZLONOW].copy())
        q = Y[-1][0:7*N_CZLONOW].copy()
        dq = Y[-1][7*N_CZLONOW:14*N_CZLONOW].copy()
        t += SEGMENT

    return ukl, klatki


def eksportuj(klatki, co_ile=8, plik='web/dane_bokser.js'):
    dane_klatki = []
    predkosci = []
    for i, q in enumerate(klatki[::co_ile]):
        czlony = []
        for k in range(N_CZLONOW):
            r = q[3*k:3*k+3]
            p = q[3*N_CZLONOW + 4*k:3*N_CZLONOW + 4*k + 4]
            czlony.append({'r': [round(float(v), 4) for v in r],
                           'p': [round(float(v), 5) for v in p]})
        dane_klatki.append(czlony)

    dane = {
        'dt': DT*co_ile,
        'wymiary': {'ramie': L_RAMIE, 'przedramie': L_PRZEDRAMIE,
                    'tul_dol': TUL_DOL, 'tul_gora': TUL_GORA,
                    'z_bark': Z_BARK, 'rozstaw_bark': ROZSTAW_BARK},
        'indeksy': {'tulow': TUL-1, 'ra_g': RA_G-1, 'ra_p': RA_P-1,
                    'la_g': LA_G-1, 'la_p': LA_P-1},
        'klatki': dane_klatki,
    }
    katalog = os.path.dirname(os.path.abspath(plik))
    os.makedirs(katalog, exist_ok=True)
    with open(plik, 'w') as f:
        f.write('const DANE = ' + json.dumps(dane) + ';\n')
    print(f'zapisano {len(dane_klatki)} klatek do {plik}')


if __name__ == '__main__':
    ukl, klatki = symuluj()
    # predkosc piesci (szczyt = moment ciosu)
    dt_kl = DT
    piesci = np.array([piesc_prawa(q) for q in klatki])
    predk = np.linalg.norm(np.diff(piesci, axis=0), axis=1)/dt_kl
    print(f'liczba klatek: {len(klatki)}')
    print(f'zasieg prawej piesci: x od {piesci[:,0].min():.2f} do {piesci[:,0].max():.2f} m')
    print(f'szczytowa predkosc piesci: {predk.max():.1f} m/s')

    sciezka = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           '..', 'web', 'dane_bokser.js')
    eksportuj(klatki, plik=os.path.normpath(sciezka))
