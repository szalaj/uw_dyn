# -*- coding: utf-8 -*-
# Przyklad: kickboxing na worku. Kombinacja jab-prosty-sierpowy-kopniecie (3D).
#
# Worek bokserski jest PRAWDZIWYM CIALEM: cieznym wahadlem (22 kg)
# podwieszonym stawem kulistym z TLUMIENIEM przegubu (MomentSferyczny k=0,
# c>0), wiec po uderzeniu wychyla sie i realistycznie wraca do pionu, a nie
# buja bez konca. Ciosy trafiaja przez kontakt penalty bryla-bryla
# (SilaUderzenia: punkt piesci/stopy vs kapsula worka), uderzenie przekazuje
# ped i mierzy sile/impuls. Dodatkowo tulow vs worek jako kontakt kapsula-
# kapsula, dzieki czemu odbity worek nie przenika przez boksera. To realne
# trafianie w cel z bezwladnoscia, nie walka z cieniem.
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
# Nogi: udo + podudzie z biodrem i stawem kolanowym (zawiasy wokol osi y),
# w ugietej postawie atletycznej (masy z antropometrii). Miednica przypieta,
# noga podporowa trzyma postawe podczas kopniecia.
#
# Sekwencja (kombinacja treningowa): jab (lewa prosta) -> prosty (prawa) ->
# sierpowy (prawa) -> front kick tylna noga. Ciosy szybkie (mocne aktuatory,
# krotkie fazy T_OUT/T_BACK), kazdy z osobna metryka uderzenia (szczytowa sila
# i impuls, mierzone w oknie ciosu). Glowa rysowana pogladowo.
#
# Wynik: web/dane_bokser.js do wizualizacji Three.js (web/bokser.html).

import json
import os

import numpy as np

from uw_dyn import (Uklad, Czlon, Polaczenie_Obr, Para_Sferyczna,
                    MomentWzgledny, MomentSferyczny, SilaUderzenia,
                    wektor, u2p, R, wektor_p,
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
WOREK = 10                       # worek bokserski jako prawdziwe cialo (wahadlo)
N_CZLONOW = 10

# ----- worek bokserski (wahadlo z masa, podwieszone stawem kulistym) -----
M_WOREK = 22.0                   # masa worka [kg] (realistyczny worek treningowy)
R_WOREK = 0.22                   # promien worka [m] (szeroki, by hak i kopniecie sie zmiescily)
H_WOREK = 0.90                   # wysokosc worka [m]
X_WOREK = 0.50                   # odleglosc worka przed bokserem (os x) [m]
Z_ZACZEP = 1.55                  # przegub kulisty (gora worka) [m]
Z_WOREK = Z_ZACZEP - H_WOREK/2   # srodek masy worka w spoczynku
PROM_KONTAKT = R_WOREK + PROM    # promien kontaktu piesc/stopa - worek
C_PRZEGUB = 30.0                 # tlumienie przegubu (wahania zanikaja, powrot do pionu)

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
# proste ciosy: piesc leci do przodu (+x), lokiec niemal wyprostowany (flex~0)
J_L = dict(dz=(1.0, 0.12, 0.05), dx=(0.0, 0.4, 1.0), flex=-0.25)   # jab (lewa)
C_R = dict(dz=(1.0, -0.12, 0.05), dx=(0.0, -0.4, 1.0), flex=0.25)  # prosty (prawa)

YAW_TUL_G = 0.28
YAW_TUL_H = 0.60     # skret korpusu przy ciosach prawa reka (moc z bioder)
YAW_TUL_J = 0.10     # lekki kontr-skret przy jabie (lewa do przodu)

# ----- sterownik -----
SEGMENT = 0.004
DT = 0.0005
K_BARK, C_BARK = 170.0, 11.0    # mocny bark -> szybki, twardy cios
K_LOKIEC, C_LOKIEC = 130.0, 6.0
K_TUL, C_TUL = 500.0, 20.0
K_NOGA, C_NOGA = 450.0, 22.0    # biodro/kolano: ciezka noga, mocny naped
T_GARDA = 0.45       # ustawienie gardy na starcie
T_OUT = 0.09         # wyprowadzenie ciosu (szybkie)
T_BACK = 0.15        # powrot piesci do gardy
# kombinacja: jab (L) -> prosty (R) -> sierpowy (R) -> front kick
COMBO = [
    dict(nazwa='jab',      reka='L', poza=J_L, yaw=YAW_TUL_J, t=T_GARDA),
    dict(nazwa='prosty',   reka='R', poza=C_R, yaw=YAW_TUL_H, t=T_GARDA + 0.32),
    dict(nazwa='sierpowy', reka='R', poza=H_R, yaw=YAW_TUL_H, t=T_GARDA + 0.66),
]
for _c in COMBO:
    _c['t_end'] = _c['t'] + T_OUT + T_BACK
# front kick (tylna noga): kolano w gore, wyprost, zlozenie, powrot
T_CHAMBER = 0.18
T_STRIKE = 0.10
T_ZLOZ = 0.10
T_KICK_POW = 0.34
T_KONIEC = 0.6
T_KICK0 = COMBO[-1]['t_end'] + 0.15   # kopniecie po ciosach recznych
CZAS = T_KICK0 + T_CHAMBER + T_STRIKE + T_ZLOZ + T_KICK_POW + T_KONIEC


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


def tensor_worka():
    """Tensor walca (worek) wzdluz lokalnej z: promien R_WOREK, wys. H_WOREK."""
    Jo = M_WOREK*R_WOREK**2/2
    Jp = M_WOREK*(3*R_WOREK**2 + H_WOREK**2)/12
    return np.diag([Jp, Jp, Jo])


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
                                               K_NOGA, 0.0, C_NOGA)
        ukl.dodajWiez(Polaczenie_Obr(ud, pd, wektor(0, 0, L_UD/2),
                                     wektor(0, 0, -L_PD/2), *OSIE_NOGA))
        akt['kolano_'+strona] = MomentWzgledny(ud, pd, wektor(0, 1, 0), REF_Z,
                                               K_NOGA, 0.0, C_NOGA)

    for a in akt.values():
        ukl.dodajSileWewn(a)

    # worek bokserski: cialo z masa podwieszone przegubem kulistym do zaczepu,
    # z tlumieniem obrotowym w przegubie (wahania zanikaja -> realistyczny powrot)
    ukl.dodajCzlon(Czlon(WOREK, M_WOREK, tensor_worka()))
    ukl.dodajWiez(Para_Sferyczna(0, WOREK, wektor(X_WOREK, 0, Z_ZACZEP),
                                 wektor(0, 0, H_WOREK/2)))
    ukl.dodajSileWewn(MomentSferyczny(0, WOREK, 0.0, C_PRZEGUB))  # tlumik przegubu
    # kontakt penalty bryla-bryla: obie piesci i stopa (P) uderzaja w worek
    kontakty = {
        'piesc_L': SilaUderzenia(LA_P, wektor(0, 0, L_FA/2), WOREK,
                                 PROM_KONTAKT, H_WOREK/2),
        'piesc_R': SilaUderzenia(RA_P, wektor(0, 0, L_FA/2), WOREK,
                                 PROM_KONTAKT, H_WOREK/2),
        'stopa': SilaUderzenia(PD_P, wektor(0, 0, L_PD/2), WOREK,
                               PROM_KONTAKT, H_WOREK/2),
    }
    for k in kontakty.values():
        ukl.dodajSileWewn(k)
    # kontakt tulow-worek (kapsula-kapsula): worek odpycha sie od ciala i nie
    # przenika przez boksera przy odbiciu
    ukl.dodajSileWewn(SilaUderzenia(TUL, wektor(0, 0, 0), WOREK,
                                    R_WOREK + 0.16, H_WOREK/2,
                                    polowa_wys_i=TUL_H/2, k=2.0e4, c=200.0))
    ukl.grawitacja = True

    # cele bioder/kolan = katy w postawie neutralnej (punkt staly, bez zgadywania)
    q_neu = q_startowe(YAW_TUL_G, {s: (orientacja(G[s]['dz'], G[s]['dx']),
                                       G[s]['flex']) for s in ('R', 'L')})
    for nazwa in ('biodro_L', 'biodro_P', 'kolano_L', 'kolano_P'):
        akt[nazwa].theta_cel = akt[nazwa].kat(q_neu, N_CZLONOW)
    return ukl, akt, kontakty


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

    # worek: srodek masy pod zaczepem, kwaternion jednostkowy (wisi pionowo)
    q[3*(WOREK-1):3*(WOREK-1)+3] = [X_WOREK, 0, Z_WOREK]
    q[b + 4*(WOREK-1)] = 1.0
    return q


def _ramp(t, t0, dt_faza):
    if t <= t0:
        return 0.0
    if t >= t0 + dt_faza:
        return 1.0
    return (1 - np.cos(np.pi*(t - t0)/dt_faza))/2


def _akt(t, t0):
    """Aktywacja ciosu 0->1->0: szybkie wyprowadzenie T_OUT, powrot T_BACK."""
    if t <= t0 or t >= t0 + T_OUT + T_BACK:
        return 0.0
    if t < t0 + T_OUT:
        return _ramp(t, t0, T_OUT)
    return 1.0 - _ramp(t, t0 + T_OUT, T_BACK)


def _pozy_cache():
    """Prekalkulacja orientacji gardy i poz ciosow (parametry Eulera)."""
    c = {'g': {}}
    for s in ('L', 'R'):
        c['g'][s] = (orientacja(G[s]['dz'], G[s]['dx']), G[s]['flex'])
    for strike in COMBO:
        c[strike['nazwa']] = (orientacja(strike['poza']['dz'], strike['poza']['dx']),
                              strike['poza']['flex'])
    return c


def sterowanie(t, akt, cache):
    """Kombinacja ciosow: garda + aktywny cios (jab/prosty/sierpowy)."""
    yaw = YAW_TUL_G
    cele = {s: cache['g'][s] for s in ('L', 'R')}
    for strike in COMBO:
        a = _akt(t, strike['t'])
        if a <= 0.0:
            continue
        reka = strike['reka']
        cele[reka] = (slerp(cache['g'][reka][0], cache[strike['nazwa']][0], a),
                      _lerp(cache['g'][reka][1], cache[strike['nazwa']][1], a))
        yaw = YAW_TUL_G + (strike['yaw'] - YAW_TUL_G)*a
    akt['tul'].theta_cel = yaw
    for s in ('L', 'R'):
        akt['bark_'+s].p_cel = cele[s][0]
        akt['lokiec_'+s].theta_cel = cele[s][1]


def _lerp(a, b, f):
    return a + (b - a)*f


def sterowanie_kopniecie(t, akt, kick):
    """Front kick tylna noga (P): chamber -> wyprost -> zloz -> powrot.
    kick = {'stance':(hip,knee), 'chamber':..., 'strike':...}."""
    st, ch, sr = kick['stance'], kick['chamber'], kick['strike']
    tl = t - T_KICK0
    if tl < 0:
        para = st
    elif tl < T_CHAMBER:                      # stance -> chamber
        f = _ramp(tl, 0, T_CHAMBER)
        para = (_lerp(st[0], ch[0], f), _lerp(st[1], ch[1], f))
    elif tl < T_CHAMBER + T_STRIKE:           # chamber -> wyprost (uderzenie)
        f = _ramp(tl, T_CHAMBER, T_STRIKE)
        para = (_lerp(ch[0], sr[0], f), _lerp(ch[1], sr[1], f))
    elif tl < T_CHAMBER + T_STRIKE + T_ZLOZ:  # wyprost -> chamber (zloz)
        f = _ramp(tl, T_CHAMBER + T_STRIKE, T_ZLOZ)
        para = (_lerp(sr[0], ch[0], f), _lerp(sr[1], ch[1], f))
    else:                                     # chamber -> stance (powrot)
        f = _ramp(tl, T_CHAMBER + T_STRIKE + T_ZLOZ, T_KICK_POW)
        para = (_lerp(ch[0], st[0], f), _lerp(ch[1], st[1], f))
    akt['biodro_P'].theta_cel, akt['kolano_P'].theta_cel = para


def cele_kopniecia(akt):
    """Cele hip/kolano tylnej nogi (P) dla klatek: stance, chamber, uderzenie.
    Liczone jako kat zmierzony w danej pozie (punkt staly, bez zgadywania)."""
    def q_noga(lean, knee):
        q = q_startowe(YAW_TUL_G, {s: (orientacja(G[s]['dz'], G[s]['dx']),
                                       G[s]['flex']) for s in ('R', 'L')})
        b = 3*N_CZLONOW
        p_tul = q[b + 4*(TUL-1):b + 4*(TUL-1)+4]
        q[b + 4*(UD_P-1):b + 4*(UD_P-1)+4] = mnoz_kwaterniony(p_tul, u2p(OS_Y, np.pi-lean))
        q[b + 4*(PD_P-1):b + 4*(PD_P-1)+4] = mnoz_kwaterniony(p_tul, u2p(OS_Y, np.pi-lean+knee))
        return q
    def para(lean, knee):
        q = q_noga(lean, knee)
        return (akt['biodro_P'].kat(q, N_CZLONOW), akt['kolano_P'].kat(q, N_CZLONOW))
    return {'stance': para(HIP_FLEX, KNEE_FLEX),
            'chamber': para(1.5, 1.7),
            'strike': para(1.5, 0.15)}


def stopa_prawa(q):
    """Pozycja prawej stopy (koniec podudzia P)."""
    r = q[3*(PD_P-1):3*(PD_P-1)+3]
    p = q[3*N_CZLONOW + 4*(PD_P-1):3*N_CZLONOW + 4*(PD_P-1)+4]
    return r + _os_z(p)*L_PD/2


def piesc(q, fa):
    r = q[3*(fa-1):3*(fa-1)+3]
    p = q[3*N_CZLONOW + 4*(fa-1):3*N_CZLONOW + 4*(fa-1)+4]
    return r + _os_z(p)*L_FA/2


def piesc_prawa(q):
    return piesc(q, RA_P)


def symuluj():
    ukl, akt, kontakty = zbuduj()
    pozy_g = {s: (orientacja(G[s]['dz'], G[s]['dx']), G[s]['flex']) for s in ('R', 'L')}
    q = ukl.projekcja_polozen(q_startowe(YAW_TUL_G, pozy_g))
    dq = np.zeros(7*N_CZLONOW)
    cache = _pozy_cache()
    kick = cele_kopniecia(akt)

    # okna ciosow do metryk per cios: (nazwa, kontakt, t_start, t_end)
    okna = [(c['nazwa'], 'piesc_L' if c['reka'] == 'L' else 'piesc_R',
             c['t'], c['t_end']) for c in COMBO]
    okna.append(('kopniecie', 'stopa', T_KICK0,
                 T_KICK0 + T_CHAMBER + T_STRIKE + 0.06))
    wyniki = {}
    zres, zzap = set(), set()

    klatki = []
    t = 0.0
    for _ in range(int(CZAS/SEGMENT)):
        sterowanie(t, akt, cache)
        sterowanie_kopniecie(t, akt, kick)
        for nazwa, kkey, ts, te in okna:      # reset na starcie, zapis na koncu okna
            if nazwa not in zres and t >= ts:
                kontakty[kkey].zeruj_metryki()
                zres.add(nazwa)
            if nazwa not in zzap and t >= te:
                k = kontakty[kkey]
                wyniki[nazwa] = (k.F_szczyt, k.impuls(DT))
                zzap.add(nazwa)
        ukl.sym2(np.concatenate((q, dq)), 0.0, SEGMENT, DT)
        Y = ukl.Y
        for w in Y[:-1]:
            klatki.append(w[0:7*N_CZLONOW].copy())
        q = Y[-1][0:7*N_CZLONOW].copy()
        dq = Y[-1][7*N_CZLONOW:14*N_CZLONOW].copy()
        t += SEGMENT
    for nazwa, kkey, ts, te in okna:          # okna niedomkniete przed koncem
        if nazwa not in zzap:
            k = kontakty[kkey]
            wyniki[nazwa] = (k.F_szczyt, k.impuls(DT))
    return ukl, klatki, wyniki


def eksportuj(klatki, wyniki=None, co_ile=8, plik='web/dane_bokser.js'):
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
                    'ud_L': UD_L-1, 'pd_L': PD_L-1, 'ud_P': UD_P-1, 'pd_P': PD_P-1,
                    'worek': WOREK-1},
        'worek': {'r': R_WOREK, 'h': H_WOREK, 'zaczep': [X_WOREK, 0, Z_ZACZEP]},
        'klatki': dane_klatki,
    }
    if wyniki is not None:
        dane['metryki'] = {nazwa: {'F': round(F, 0), 'imp': round(imp, 2)}
                           for nazwa, (F, imp) in wyniki.items()}
    katalog = os.path.dirname(os.path.abspath(plik))
    os.makedirs(katalog, exist_ok=True)
    with open(plik, 'w') as f:
        f.write('const DANE = ' + json.dumps(dane) + ';\n')
    print(f'zapisano {len(dane_klatki)} klatek do {plik}')


def _szczyt_predkosci(klatki, punkt_fn, t0, t1):
    """Szczytowa predkosc wskazanego punktu w oknie czasu [t0, t1]."""
    i0, i1 = int(t0/DT), min(int(t1/DT), len(klatki)-1)
    P = np.array([punkt_fn(klatki[i]) for i in range(i0, i1+1)])
    if len(P) < 2:
        return 0.0
    return float(np.linalg.norm(np.diff(P, axis=0), axis=1).max()/DT)


if __name__ == '__main__':
    ukl, klatki, wyniki = symuluj()

    punkt = {'jab': lambda q: piesc(q, LA_P), 'prosty': lambda q: piesc(q, RA_P),
             'sierpowy': lambda q: piesc(q, RA_P), 'kopniecie': stopa_prawa}
    print('kombinacja na worku (jab-prosty-sierpowy-kopniecie):')
    for c in COMBO:
        F, imp = wyniki[c['nazwa']]
        v = _szczyt_predkosci(klatki, punkt[c['nazwa']], c['t'], c['t_end'])
        print('  %-9s v=%.1f m/s, sila %6.0f N, impuls %5.2f N s'
              % (c['nazwa'], v, F, imp))
    Fk, impk = wyniki['kopniecie']
    vk = _szczyt_predkosci(klatki, stopa_prawa, T_KICK0, CZAS)
    print('  %-9s v=%.1f m/s, sila %6.0f N, impuls %5.2f N s'
          % ('kopniecie', vk, Fk, impk))

    worek = np.array([q[3*(WOREK-1):3*(WOREK-1)+3] for q in klatki])
    wych = np.hypot(worek[:, 0] - X_WOREK, worek[:, 1])
    print('worek: maks. wychylenie CoM %.3f m (koncowe %.3f m, wraca do pionu)'
          % (wych.max(), wych[-1]))
    print('brak NaN:', not np.isnan(np.array(klatki)).any())

    sciezka = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           '..', 'web', 'dane_bokser.js')
    eksportuj(klatki, wyniki, plik=os.path.normpath(sciezka))
