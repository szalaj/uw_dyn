# -*- coding: utf-8 -*-
# Przyklad: kickboxing - walka z cieniem. Prawy sierpowy (model 3D).
#
# Wersja anatomiczna na stawach kulistych (Etap A): barki to napedzane stawy
# kuliste 3 DOF (Para_Sferyczna + MomentSferyczny), lokcie to zawiasy 1 DOF
# (Polaczenie_Obr + MomentWzgledny), a skret korpusu w pasie to przegub
# obrotowy wokol pionu. Dzieki stawom kulistym mozliwa jest naturalna garda
# (lokcie w dol, piesci przy twarzy) i pelny trojwymiarowy sierpowy, czego
# nie dawala poprzednia wersja w plaszczyznie poziomej.
#
# Czlony ramion to prety wzdluz lokalnej osi z (bark -> lokiec -> piesc).
# Orientacje ramion zadaje sie kierunkiem (funkcja `orientacja`), a cele
# sterowania interpoluje slerpem miedzy garda a sierpowym.
#
# Nogi: udo + podudzie z biodrem i STAWEM KOLANOWYM (zawiasy wokol osi y),
# w ugietej postawie atletycznej (masy z antropometrii). Miednica przypieta,
# wiec nogi na razie trzymaja postawe (baza pod przyszle kopniecia/kolana).
#
# Sekwencja: garda -> wyprowadzenie ciosu -> powrot. Glowa rysowana pogladowo.
#
# Wynik: web/dane_bokser.js do wizualizacji Three.js (web/bokser.html).

import json
import os

import numpy as np

from uw_dyn import (Uklad, Czlon, Polaczenie_Obr, Para_Sferyczna,
                    MomentWzgledny, MomentSferyczny, wektor, u2p, R, wektor_p,
                    mnoz_kwaterniony, sprzezenie_kwaternionu,
                    kwaternion_na_wektor_obrotu, macierz_na_kwaternion)

# ----- wymiary (metry, kilogramy) -----
TUL_DOL, TUL_GORA = 0.95, 1.50
Z_TUL = (TUL_GORA + TUL_DOL)/2
TUL_H = TUL_GORA - TUL_DOL
M_TUL = 34.0
BARK_Z = 1.42 - Z_TUL           # bark w ukladzie tulowia (skladowa z)
BARK_Y = 0.19                   # polowa rozstawu barkow
L_UA, L_FA = 0.30, 0.28         # ramie (bark-lokiec), przedramie (lokiec-piesc)
M_UA, M_FA = 2.6, 1.6
PROM = 0.05
# nogi: udo (biodro-kolano) i podudzie (kolano-stopa); atletyczna postawa
L_UD, L_PD = 0.45, 0.45         # dlugosci czlonow nogi [m]
M_UD, M_PD = 7.8, 3.6           # masy (antropometria, ~78 kg)
BIODRO_Y = 0.12                 # polowa rozstawu bioder
HIP_FLEX = 0.20                 # ugiecie biodra (udo lekko do przodu) [rad]
KNEE_FLEX = 0.45                # ugiecie kolana w postawie [rad]

TUL, RA_G, RA_P, LA_G, LA_P = 1, 2, 3, 4, 5
UD_L, PD_L, UD_P, PD_P = 6, 7, 8, 9
N_CZLONOW = 9

OS_Z = np.array([0.0, 0.0, 1.0])
OS_Y = np.array([0.0, 1.0, 0.0])
OSIE_Z = (wektor(1, 0, 0), wektor(0, 1, 0), wektor(0, 0, 1))   # przegub wokol pionu (pas)
# zawias lokcia: os = lokalna x ramienia (vi=y, wi=z -> vi x wi = x); uj = x przedramienia
OSIE_LOKIEC = (wektor(0, 1, 0), wektor(0, 0, 1), wektor(1, 0, 0))
# zawias biodra/kolana: os = lokalna y (zgiecie w plaszczyznie strzalkowej)
OSIE_NOGA = (wektor(0, 0, 1), wektor(1, 0, 0), wektor(0, 1, 0))
REF_X = wektor(1, 0, 0)
REF_Z = wektor(0, 0, 1)

# ----- pozy (w ukladzie tulowia): kierunek ramienia dz, plaszczyzna lokcia dx, zgiecie -----
# garda: lokcie w dol, piesci przy policzkach
G = {
    'R': dict(dz=(0.30, 0.20, -1.0), dx=(0.0, -0.5, 0.6), flex=2.60),
    # lewa reka to lustro prawej; odbicie zmienia skretnosc, wiec zgiecie
    # lokcia ma znak ujemny (inaczej przedramie wygina sie nienaturalnie w tyl)
    'L': dict(dz=(0.30, -0.20, -1.0), dx=(0.0, 0.5, 0.6), flex=-2.60),
}
# prawy sierpowy: ramie uniesione (lokiec do wysokosci barku), piesc lukiem
# w przod-w poprzek; korpus obraca sie CCW (prawy bark do przodu), yaw rosnie
H_R = dict(dz=(0.50, -0.20, 0.10), dx=(0.0, 0.3, 1.0), flex=1.70)

YAW_TUL_G = 0.28
YAW_TUL_H = 0.60

# ----- sterownik -----
SEGMENT = 0.004
DT = 0.0005
K_BARK, C_BARK = 90.0, 7.0
K_LOKIEC, C_LOKIEC = 60.0, 3.0
K_TUL, C_TUL = 400.0, 16.0
T_GARDA = 0.6
T_CIOS = 0.26        # wyprowadzenie + trafienie
T_POWROT = 0.34
T_KONIEC = 0.7
CZAS = T_GARDA + T_CIOS + T_POWROT + T_KONIEC


# ---------- pomocnicze: orientacje i slerp ----------
def _norm(v):
    v = np.asarray(v, dtype=float)
    return v/np.linalg.norm(v)


def orientacja(dz, dx_hint):
    """Parametry Eulera obrotu, ktorego lokalna os z -> dz, a os x ~ dx_hint."""
    z = _norm(dz)
    x = np.asarray(dx_hint, float) - np.dot(dx_hint, z)*z
    x = _norm(x)
    y = np.cross(z, x)
    return macierz_na_kwaternion(np.column_stack([x, y, z]))


def slerp(pa, pb, t):
    q_rel = mnoz_kwaterniony(sprzezenie_kwaternionu(pa), pb)
    v = kwaternion_na_wektor_obrotu(q_rel)
    kat = np.linalg.norm(v)
    if kat < 1e-9:
        return np.asarray(pa, float)
    return mnoz_kwaterniony(pa, u2p(v/kat, kat*t))


def _os_z(p):
    return R(wektor_p(*p))[:, 2]


# ---------- model ----------
def tensor_preta(m, L):
    Jo = m*PROM**2/2
    Jp = m*(3*PROM**2 + L**2)/12
    return np.diag([Jp, Jp, Jo])   # pret wzdluz lokalnej z


def tensor_tulowia():
    Jx = M_TUL*(0.32**2 + TUL_H**2)/12
    Jy = M_TUL*(0.22**2 + TUL_H**2)/12
    Jz = M_TUL*(0.32**2 + 0.22**2)/12
    return np.diag([Jx, Jy, Jz])


def bark_w_tulowiu(strona):
    return wektor(0, BARK_Y if strona == 'L' else -BARK_Y, BARK_Z)


def biodro_w_tulowiu(strona):
    return wektor(0, BIODRO_Y if strona == 'L' else -BIODRO_Y, -TUL_H/2)


# orientacje nog w postawie (obrot wokol y; lokalna os z w dol)
def _p_udo():
    return u2p(OS_Y, np.pi - HIP_FLEX)       # udo lekko do przodu


def _p_podudzie():
    return u2p(OS_Y, np.pi - HIP_FLEX + KNEE_FLEX)  # podudzie cofniete (kolano zgiete)


def zbuduj():
    ukl = Uklad()
    ukl.dodajCzlon(Czlon(TUL, M_TUL, tensor_tulowia()))
    ukl.dodajCzlon(Czlon(RA_G, M_UA, tensor_preta(M_UA, L_UA)))
    ukl.dodajCzlon(Czlon(RA_P, M_FA, tensor_preta(M_FA, L_FA)))
    ukl.dodajCzlon(Czlon(LA_G, M_UA, tensor_preta(M_UA, L_UA)))
    ukl.dodajCzlon(Czlon(LA_P, M_FA, tensor_preta(M_FA, L_FA)))
    for nr, m, L in ((UD_L, M_UD, L_UD), (PD_L, M_PD, L_PD),
                     (UD_P, M_UD, L_UD), (PD_P, M_PD, L_PD)):
        ukl.dodajCzlon(Czlon(nr, m, tensor_preta(m, L)))

    # pas: przegub obrotowy wokol pionu (skret korpusu)
    ukl.dodajWiez(Polaczenie_Obr(0, TUL, wektor(0, 0, TUL_DOL),
                                 wektor(0, 0, -TUL_H/2), *OSIE_Z))
    akt = {}
    akt['tul'] = MomentWzgledny(0, TUL, wektor(0, 0, 1), REF_X, K_TUL, YAW_TUL_G, C_TUL)

    for strona, (ua, fa) in (('R', (RA_G, RA_P)), ('L', (LA_G, LA_P))):
        # bark: staw kulisty (polozenie) + moment 3D (orientacja ramienia)
        ukl.dodajWiez(Para_Sferyczna(TUL, ua, bark_w_tulowiu(strona), wektor(0, 0, -L_UA/2)))
        p_bark = orientacja(G[strona]['dz'], G[strona]['dx'])
        akt['bark_'+strona] = MomentSferyczny(TUL, ua, K_BARK, C_BARK, p_cel=p_bark)
        # lokiec: zawias (os = lokalna x ramienia) + moment 1 DOF (zgiecie)
        ukl.dodajWiez(Polaczenie_Obr(ua, fa, wektor(0, 0, L_UA/2),
                                     wektor(0, 0, -L_FA/2), *OSIE_LOKIEC))
        akt['lokiec_'+strona] = MomentWzgledny(ua, fa, wektor(1, 0, 0), wektor(0, 0, 1),
                                               K_LOKIEC, G[strona]['flex'], C_LOKIEC)

    # nogi: biodro (zawias, os y) tulow->udo, kolano (zawias, os y) udo->podudzie
    for strona, (ud, pd) in (('L', (UD_L, PD_L)), ('P', (UD_P, PD_P))):
        ukl.dodajWiez(Polaczenie_Obr(TUL, ud, biodro_w_tulowiu(strona),
                                     wektor(0, 0, -L_UD/2), *OSIE_NOGA))
        akt['biodro_'+strona] = MomentWzgledny(TUL, ud, wektor(0, 1, 0), REF_Z,
                                               K_BARK, 0.0, C_BARK)
        ukl.dodajWiez(Polaczenie_Obr(ud, pd, wektor(0, 0, L_UD/2),
                                     wektor(0, 0, -L_PD/2), *OSIE_NOGA))
        akt['kolano_'+strona] = MomentWzgledny(ud, pd, wektor(0, 1, 0), REF_Z,
                                               K_BARK, 0.0, C_BARK)

    for a in akt.values():
        ukl.dodajSileWewn(a)
    ukl.grawitacja = True

    # cele bioder/kolan = katy w postawie neutralnej (punkt staly, bez zgadywania)
    q_neu = q_startowe(YAW_TUL_G, {s: (orientacja(G[s]['dz'], G[s]['dx']),
                                       G[s]['flex']) for s in ('R', 'L')})
    for nazwa in ('biodro_L', 'biodro_P', 'kolano_L', 'kolano_P'):
        akt[nazwa].theta_cel = akt[nazwa].kat(q_neu, N_CZLONOW)
    return ukl, akt


def q_startowe(yaw_tul, pozy):
    """FK pozy: pozycje srodkow mas i kwaterniony. pozy[strona] = (p_bark_rel, flex)."""
    q = np.zeros(7*N_CZLONOW)
    p_tul = u2p(OS_Z, yaw_tul)
    Rt = R(wektor_p(*p_tul))
    q[3*(TUL-1):3*(TUL-1)+3] = [0, 0, Z_TUL]
    q[3*N_CZLONOW + 4*(TUL-1):3*N_CZLONOW + 4*(TUL-1)+4] = p_tul

    for strona, (ua, fa) in (('R', (RA_G, RA_P)), ('L', (LA_G, LA_P))):
        p_bark_rel, flex = pozy[strona]
        bark = np.array([0, 0, Z_TUL]) + Rt.dot(bark_w_tulowiu(strona).ravel())
        p_ua = mnoz_kwaterniony(p_tul, p_bark_rel)
        srodek_ua = bark + _os_z(p_ua)*L_UA/2
        lokiec = bark + _os_z(p_ua)*L_UA
        p_fa = mnoz_kwaterniony(p_ua, u2p(np.array([1.0, 0, 0]), flex))
        srodek_fa = lokiec + _os_z(p_fa)*L_FA/2

        q[3*(ua-1):3*(ua-1)+3] = srodek_ua
        q[3*(fa-1):3*(fa-1)+3] = srodek_fa
        b = 3*N_CZLONOW
        q[b + 4*(ua-1):b + 4*(ua-1)+4] = p_ua
        q[b + 4*(fa-1):b + 4*(fa-1)+4] = p_fa

    # nogi w postawie atletycznej (biodro lekko ugiete, kolano zgiete)
    b = 3*N_CZLONOW
    p_ud = mnoz_kwaterniony(p_tul, _p_udo())
    p_pd = mnoz_kwaterniony(p_tul, _p_podudzie())
    for strona, (ud, pd) in (('L', (UD_L, PD_L)), ('P', (UD_P, PD_P))):
        biodro = np.array([0, 0, Z_TUL]) + Rt.dot(biodro_w_tulowiu(strona).ravel())
        srodek_ud = biodro + _os_z(p_ud)*L_UD/2
        kolano = biodro + _os_z(p_ud)*L_UD
        srodek_pd = kolano + _os_z(p_pd)*L_PD/2
        q[3*(ud-1):3*(ud-1)+3] = srodek_ud
        q[3*(pd-1):3*(pd-1)+3] = srodek_pd
        q[b + 4*(ud-1):b + 4*(ud-1)+4] = p_ud
        q[b + 4*(pd-1):b + 4*(pd-1)+4] = p_pd
    return q


def _ramp(t, t0, dt_faza):
    if t <= t0:
        return 0.0
    if t >= t0 + dt_faza:
        return 1.0
    return (1 - np.cos(np.pi*(t - t0)/dt_faza))/2


def sterowanie(t, akt, p_bark_R_g, p_bark_R_h):
    """Aktualizacja celow: garda -> sierpowy -> powrot (prawa reka + korpus)."""
    if t < T_GARDA:
        a = 0.0
    elif t < T_GARDA + T_CIOS:
        a = _ramp(t, T_GARDA, T_CIOS)
    else:
        a = 1.0 - _ramp(t, T_GARDA + T_CIOS, T_POWROT)

    akt['tul'].theta_cel = YAW_TUL_G + (YAW_TUL_H - YAW_TUL_G)*a
    akt['bark_R'].p_cel = slerp(p_bark_R_g, p_bark_R_h, a)
    akt['lokiec_R'].theta_cel = G['R']['flex'] + (H_R['flex'] - G['R']['flex'])*a
    # lewa reka i tak trzyma garde (cele ustawione przy budowie)


def piesc_prawa(q):
    r = q[3*(RA_P-1):3*(RA_P-1)+3]
    p = q[3*N_CZLONOW + 4*(RA_P-1):3*N_CZLONOW + 4*(RA_P-1)+4]
    return r + _os_z(p)*L_FA/2


def symuluj():
    ukl, akt = zbuduj()
    pozy_g = {s: (orientacja(G[s]['dz'], G[s]['dx']), G[s]['flex']) for s in ('R', 'L')}
    q = ukl.projekcja_polozen(q_startowe(YAW_TUL_G, pozy_g))
    dq = np.zeros(7*N_CZLONOW)

    p_bark_R_g = orientacja(G['R']['dz'], G['R']['dx'])
    p_bark_R_h = orientacja(H_R['dz'], H_R['dx'])

    klatki = []
    t = 0.0
    for _ in range(int(CZAS/SEGMENT)):
        sterowanie(t, akt, p_bark_R_g, p_bark_R_h)
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
        'wymiary': {'ramie': L_UA, 'przedramie': L_FA,
                    'tul_dol': TUL_DOL, 'tul_gora': TUL_GORA,
                    'udo': L_UD, 'podudzie': L_PD},
        'indeksy': {'tulow': TUL-1, 'ra_g': RA_G-1, 'ra_p': RA_P-1,
                    'la_g': LA_G-1, 'la_p': LA_P-1,
                    'ud_L': UD_L-1, 'pd_L': PD_L-1, 'ud_P': UD_P-1, 'pd_P': PD_P-1},
        'klatki': dane_klatki,
    }
    katalog = os.path.dirname(os.path.abspath(plik))
    os.makedirs(katalog, exist_ok=True)
    with open(plik, 'w') as f:
        f.write('const DANE = ' + json.dumps(dane) + ';\n')
    print(f'zapisano {len(dane_klatki)} klatek do {plik}')


if __name__ == '__main__':
    ukl, klatki = symuluj()
    piesci = np.array([piesc_prawa(q) for q in klatki])
    predk = np.linalg.norm(np.diff(piesci, axis=0), axis=1)/DT
    i_g = int(round((T_GARDA*0.5)/DT))
    i_h = int(np.argmax(piesci[:, 0]))
    print('garda:  piesc R = (%.2f, %.2f, %.2f)' % tuple(piesci[i_g]))
    print('szczyt: piesc R = (%.2f, %.2f, %.2f) w t=%.2f s' % (*piesci[i_h], i_h*DT))
    print('szczytowa predkosc piesci: %.1f m/s' % predk.max())
    print('brak NaN:', not np.isnan(np.array(klatki)).any())

    sciezka = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           '..', 'web', 'dane_bokser.js')
    eksportuj(klatki, plik=os.path.normpath(sciezka))
