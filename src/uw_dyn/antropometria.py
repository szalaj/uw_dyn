# -*- coding: utf-8 -*-
# uw_dyn: antropometria - parametry segmentow ciala czlowieka.
# Etap A: pelna sylwetka skalowana wzrostem i masa (tablice Wintera).
#
# Zrodlo: D. A. Winter, "Biomechanics and Motor Control of Human Movement",
# tabela parametrow segmentow. Dla kazdego segmentu:
#   masa   - ulamek calkowitej masy ciala,
#   dlug   - ulamek wzrostu (dlugosc segmentu miedzy stawami),
#   com    - polozenie srodka masy od stawu bliższego (ulamek dlugosci),
#   rg     - promien bezwladnosci wzgledem SM (ulamek dlugosci, os poprzeczna),
#   prom   - promien poprzeczny segmentu (ulamek dlugosci; do osi podluznej
#            i wizualizacji).
# Dlon zostala scalona z przedramieniem (masa + dlugosc), stad "przedramie".
# Tensor bezwladnosci: pret/walec wzdluz lokalnej osi z (jak w reszcie
# biblioteki); osie poprzeczne z tablicy rg, os podluzna z przyblizenia walca.

from dataclasses import dataclass

import numpy as np

# nazwa: (masa, dlug, com, rg, prom)  -- ulamki (Winter, dlon w przedramieniu)
TABELA_WINTER = {
    'glowa':      (0.081, 0.200, 0.500, 0.495, 0.45),
    'tulow':      (0.497, 0.300, 0.500, 0.320, 0.45),
    'ramie':      (0.028, 0.186, 0.436, 0.322, 0.09),
    'przedramie': (0.022, 0.160, 0.450, 0.300, 0.08),
    'udo':        (0.100, 0.245, 0.433, 0.323, 0.10),
    'podudzie':   (0.0465, 0.246, 0.433, 0.302, 0.08),
    'stopa':      (0.0145, 0.152, 0.500, 0.475, 0.15),
}

# segmenty parzyste (lewy/prawy) budowane po obu stronach
SEGMENTY_PARZYSTE = ('ramie', 'przedramie', 'udo', 'podudzie', 'stopa')


@dataclass
class Segment:
    """Parametry pojedynczego segmentu ciala (jednostki SI)."""
    nazwa: str
    masa: float          # [kg]
    dlugosc: float       # [m]
    com: float           # polozenie SM od stawu bliższego (ulamek dlugosci)
    promien: float       # promien poprzeczny [m]
    tensor: np.ndarray   # tensor bezwladnosci wzgledem SM (3x3), pret wzdluz z


def tensor_segmentu(masa, dlugosc, rg, promien):
    """Tensor bezwladnosci segmentu (pret/walec wzdluz lokalnej osi z).

    Osie poprzeczne (x, y): I = masa*(rg*dlugosc)^2 (promien bezwladnosci z
    tablicy). Os podluzna (z): przyblizenie walca I = 0.5*masa*promien^2."""
    I_poprz = masa*(rg*dlugosc)**2
    I_podl = 0.5*masa*promien**2
    return np.diag([I_poprz, I_poprz, I_podl])


def segment(nazwa, masa_ciala, wzrost):
    """Segment o zadanej nazwie dla osoby (masa_ciala [kg], wzrost [m])."""
    frac_m, frac_l, com, rg, frac_p = TABELA_WINTER[nazwa]
    masa = frac_m*masa_ciala
    dlugosc = frac_l*wzrost
    promien = frac_p*dlugosc
    return Segment(nazwa, masa, dlugosc, com, promien,
                   tensor_segmentu(masa, dlugosc, rg, promien))


def segmenty(masa_ciala, wzrost):
    """Wszystkie segmenty ciala (parzyste z sufiksem _L / _P).

    Zwraca slownik nazwa -> Segment. Suma mas segmentow = masa_ciala."""
    wynik = {}
    for nazwa in TABELA_WINTER:
        s = segment(nazwa, masa_ciala, wzrost)
        if nazwa in SEGMENTY_PARZYSTE:
            wynik[nazwa+'_L'] = s
            wynik[nazwa+'_P'] = Segment(nazwa+'_P', s.masa, s.dlugosc, s.com,
                                        s.promien, s.tensor)
            wynik[nazwa+'_L'] = Segment(nazwa+'_L', s.masa, s.dlugosc, s.com,
                                        s.promien, s.tensor)
        else:
            wynik[nazwa] = s
    return wynik


def masa_calkowita(masa_ciala, wzrost):
    """Suma mas wszystkich segmentow (kontrola: ~ masa_ciala)."""
    return sum(s.masa for s in segmenty(masa_ciala, wzrost).values())


# ---------------------------------------------------------------------------
# Builder pelnej sylwetki (stojaca postac) jako Uklad wieloczlonowy
# ---------------------------------------------------------------------------

def zbuduj_postac(masa_ciala=75.0, wzrost=1.80, k_staw=300.0, c_staw=15.0):
    """Buduje stojaca postac czlowieka jako Uklad wieloczlonowy.

    Segmenty: tulow (podstawa, przypieta do podloza), glowa, po obu stronach
    ramie+przedramie oraz udo+podudzie+stopa (11 czlonow). Stawy: bark, biodro,
    szyja i przypiecie tulowia jako kuliste (Para_Sferyczna + MomentSferyczny);
    lokiec, kolano, kostka jako zawiasy (Polaczenie_Obr + MomentWzgledny).
    Aktuatory trzymaja poze neutralna (stojaca).

    Zwraca (uklad, nry, q0, aktuatory): nry mapuje nazwe segmentu na numer
    czlonu, q0 to wektor wspolrzednych pozy neutralnej, aktuatory to slownik
    nazwa_stawu -> aktuator (cele mozna podmieniac w sterowaniu)."""
    from uw_dyn.uklad import Uklad
    from uw_dyn.czlony import Czlon
    from uw_dyn.wiezy import Para_Sferyczna, Polaczenie_Obr
    from uw_dyn.sily import MomentSferyczny, MomentWzgledny
    from uw_dyn.algebra import (wektor, u2p, R, wektor_p,
                                mnoz_kwaterniony, sprzezenie_kwaternionu)

    S = segmenty(masa_ciala, wzrost)
    OSX = np.array([1.0, 0.0, 0.0])
    OSY = np.array([0.0, 1.0, 0.0])
    p_dol = u2p(OSX, np.pi)         # lokalna os z w dol (ramiona, nogi)
    p_id = np.array([1.0, 0.0, 0.0, 0.0])
    p_stopa = u2p(OSY, np.pi/2)     # lokalna os z do przodu (stopa)

    bark_y = 0.13*wzrost
    biodro_y = 0.09*wzrost
    L_t = S['tulow'].dlugosc
    com_t = S['tulow'].com
    P = S['udo_L'].dlugosc + S['podudzie_L'].dlugosc + 0.08   # wysokosc miednicy

    # numeracja czlonow
    kolejnosc = ['tulow', 'glowa',
                 'ramie_L', 'przedramie_L', 'ramie_P', 'przedramie_P',
                 'udo_L', 'podudzie_L', 'stopa_L',
                 'udo_P', 'podudzie_P', 'stopa_P']
    nry = {nazwa: i+1 for i, nazwa in enumerate(kolejnosc)}
    N = len(kolejnosc)

    def prox(seg):      # punkt stawu bliższego w ukladzie ciala (0,0,-com*L)
        return wektor(0, 0, -seg.com*seg.dlugosc)

    def dist(seg):      # punkt stawu dalszego (0,0,(1-com)*L)
        return wektor(0, 0, (1-seg.com)*seg.dlugosc)

    # --- kinematyka pozy neutralnej: pozycje SM i orientacje ---
    orient = {}
    poz = {}
    Z_tul = P + com_t*L_t                     # SM tulowia
    poz['tulow'] = np.array([0, 0, Z_tul]);   orient['tulow'] = p_id

    def osadz(nazwa, staw_world, p):
        s = S[nazwa]
        com_w = staw_world - R(wektor_p(*p)).dot(prox(s).ravel())
        poz[nazwa] = com_w
        orient[nazwa] = p
        return com_w + R(wektor_p(*p)).dot(dist(s).ravel())   # staw dalszy

    szyja = poz['tulow'] + R(wektor_p(*p_id)).dot(dist(S['tulow']).ravel())
    osadz('glowa', szyja, p_id)

    for bok, sy in (('L', +1), ('P', -1)):
        bark = poz['tulow'] + R(wektor_p(*p_id)).dot(
            np.array([0, sy*bark_y, (1-com_t)*L_t]))
        lokiec = osadz('ramie_'+bok, bark, p_dol)
        osadz('przedramie_'+bok, lokiec, p_dol)

        biodro = poz['tulow'] + R(wektor_p(*p_id)).dot(
            np.array([0, sy*biodro_y, -com_t*L_t]))
        kolano = osadz('udo_'+bok, biodro, p_dol)
        kostka = osadz('podudzie_'+bok, kolano, p_dol)
        osadz('stopa_'+bok, kostka, p_stopa)

    # --- wektor q0 ---
    q0 = np.zeros(7*N)
    for nazwa, nr in nry.items():
        q0[3*(nr-1):3*(nr-1)+3] = poz[nazwa]
        q0[3*N+4*(nr-1):3*N+4*(nr-1)+4] = orient[nazwa]

    # --- budowa ukladu ---
    ukl = Uklad()
    for nazwa in kolejnosc:
        ukl.dodajCzlon(Czlon(nry[nazwa], S[nazwa].masa, S[nazwa].tensor))

    OSIE_ZAWIAS = (wektor(1, 0, 0), wektor(0, 0, 1), wektor(0, 1, 0))  # os y ciala
    aktuatory = {}

    def staw_kulisty(nazwa, i, j, staw_i, staw_j, k, c):
        ukl.dodajWiez(Para_Sferyczna(i, j, staw_i, staw_j))
        a = MomentSferyczny(i, j, k, c)
        ukl.dodajSileWewn(a)
        aktuatory[nazwa] = a

    def staw_zawias(nazwa, i, j, staw_i, staw_j, k, c):
        ukl.dodajWiez(Polaczenie_Obr(i, j, staw_i, staw_j, *OSIE_ZAWIAS))
        a = MomentWzgledny(i, j, wektor(0, 1, 0), wektor(0, 0, 1), k, 0.0, c)
        ukl.dodajSileWewn(a)
        aktuatory[nazwa] = a

    t = nry['tulow']
    # przypiecie tulowia do podloza (staw kulisty w miednicy) + pozycja
    pelvis_w = poz['tulow'] + R(wektor_p(*p_id)).dot(prox(S['tulow']).ravel())
    staw_kulisty('miednica', 0, t, wektor(*pelvis_w), prox(S['tulow']),
                 8*k_staw, 4*c_staw)
    staw_kulisty('szyja', t, nry['glowa'],
                 wektor(0, 0, (1-com_t)*L_t), prox(S['glowa']), k_staw, c_staw)

    for bok, sy in (('L', +1), ('P', -1)):
        staw_kulisty('bark_'+bok, t, nry['ramie_'+bok],
                     wektor(0, sy*bark_y, (1-com_t)*L_t), prox(S['ramie_'+bok]),
                     k_staw, c_staw)
        staw_zawias('lokiec_'+bok, nry['ramie_'+bok], nry['przedramie_'+bok],
                    dist(S['ramie_'+bok]), prox(S['przedramie_'+bok]), k_staw, c_staw)
        staw_kulisty('biodro_'+bok, t, nry['udo_'+bok],
                     wektor(0, sy*biodro_y, -com_t*L_t), prox(S['udo_'+bok]),
                     2*k_staw, 2*c_staw)
        staw_zawias('kolano_'+bok, nry['udo_'+bok], nry['podudzie_'+bok],
                    dist(S['udo_'+bok]), prox(S['podudzie_'+bok]), 2*k_staw, 2*c_staw)
        staw_zawias('kostka_'+bok, nry['podudzie_'+bok], nry['stopa_'+bok],
                    dist(S['podudzie_'+bok]), prox(S['stopa_'+bok]), k_staw, c_staw)

    ukl.grawitacja = True

    # --- cele aktuatorow = poza neutralna (punkt staly) ---
    for nazwa, a in aktuatory.items():
        if isinstance(a, MomentSferyczny):
            pj = orient[_czlon_j(nazwa, nry, kolejnosc)]
            pi = (p_id if a.i == 0 else orient[kolejnosc[a.i-1]])
            a.p_cel = mnoz_kwaterniony(sprzezenie_kwaternionu(pi), pj)
        else:
            a.theta_cel = a.kat(q0, N)

    return ukl, nry, q0, aktuatory


def _czlon_j(nazwa_stawu, nry, kolejnosc):
    """Nazwa segmentu bedacego czlonem j danego stawu (do ustawienia celu)."""
    mapowanie = {
        'miednica': 'tulow', 'szyja': 'glowa',
        'bark_L': 'ramie_L', 'bark_P': 'ramie_P',
        'biodro_L': 'udo_L', 'biodro_P': 'udo_P',
    }
    return mapowanie[nazwa_stawu]
