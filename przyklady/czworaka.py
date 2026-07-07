# -*- coding: utf-8 -*-
# Przyklad: czlowiek chodzacy na czworaka (pelzanie / crawl).
#
# Pelna sylwetka o proporcjach antropometrycznych (tablice Wintera) ustawiona
# poziomo: tulow to swobodna podstawa (floating base) na wysokosci ~0.42 m,
# glowa z przodu, cztery dwuczlonowe konczyny siegajace do podloza. Konczyny
# przednie to ramiona (ramie + przedramie, kontakt dloni), tylne to nogi
# (udo + podudzie, kontakt stopy).
#
# Bark i biodro to STAWY KULISTE (Para_Sferyczna + MomentSferyczny), wiec
# konczyny moga wychodzic w bok - to konieczne, by przenosic ciezar w bok
# przy unoszeniu konczyny (bez tego postac tylko z zawiasami strzalkowymi
# nieuchronnie przewraca sie na bok, bo srodek masy laduje na krawedzi
# trojkata podparcia). Lokiec i kolano to zawiasy (Polaczenie_Obr +
# MomentWzgledny) zginajace sie w plaszczyznie konczyny.
#
# Chod: pelzanie (sekwencja boczna) - w danej chwili jedna konczyna wykonuje
# wymach, trzy pozostale podpieraja (statycznie stabilnie). Rownolegle
# tulow przenosi ciezar na strone przeciwna do unoszonej konczyny (waddle).
# W fazie podporu stopa jest "przyklejona" tarciem, a aktuatory przesuwaja
# ja do tylu wzgledem ciala (chodnik ruchomy), co pcha tulow do przodu.
#
# Sterowanie liczy w kazdym kroku pozadane polozenia stop (swiat), rozwiazuje
# przestrzenna odwrotna kinematyke kazdej konczyny (staw kulisty + zawias) i
# ustawia cele aktuatorow: orientacja wzgledna barku/biodra (p_cel) oraz kat
# lokcia/kolana (odczyt z konfiguracji docelowej przez a.kat).
#
# Wynik: web/dane_czworaka.js do wizualizacji Three.js (web/czworaka.html).

import json
import os

import numpy as np

from uw_dyn import (Uklad, Czlon, Para_Sferyczna, Polaczenie_Obr, SilaKontaktu,
                    MomentSferyczny, MomentWzgledny, segmenty,
                    wektor, u2p, R, wektor_p, mnoz_kwaterniony,
                    sprzezenie_kwaternionu, macierz_na_kwaternion)

# ----- osoba -----
MASA, WZROST = 75.0, 1.80
S = segmenty(MASA, WZROST)

OSY = np.array([0.0, 1.0, 0.0])
OSIE = (wektor(1, 0, 0), wektor(0, 0, 1), wektor(0, 1, 0))  # zawias wokol y ciala
P_TUL = u2p(OSY, np.pi/2)          # tulow poziomo: lokalna z (kregoslup) -> +x

BARK_Y = 0.13*WZROST
BIODRO_Y = 0.09*WZROST
L_TUL = S['tulow'].dlugosc
COM_TUL = S['tulow'].com

# konczyny: segment gorny/dolny, przod?, y-offset zaczepienia na tulowiu
KONCZYNY = {
    'PL': dict(gorny='ramie_L', dolny='przedramie_L', przod=True,  y=+BARK_Y),
    'PP': dict(gorny='ramie_P', dolny='przedramie_P', przod=True,  y=-BARK_Y),
    'TL': dict(gorny='udo_L',   dolny='podudzie_L',   przod=False, y=+BIODRO_Y),
    'TP': dict(gorny='udo_P',   dolny='podudzie_P',   przod=False, y=-BIODRO_Y),
}

# ----- postura i chod -----
H_TULOW = 0.40        # wysokosc srodka tulowia [m]
FLEX_KOLANO = 2.2     # podwiniecie kolana (tylne konczyny na kolanach) [rad]
BAZA_Y = 0.27         # szerokosc bazy: podpory na +/- BAZA_Y [m] (stabilnosc boczna)
STRIDE = 0.11         # dlugosc kroku (przesuw stopy wzgledem ciala) [m]
KLIRENS = 0.045       # wznios stopy w wymachu [m]
SWAY = 0.06           # przenoszenie ciezaru w bok na strone podporu [m]
RAMPA_LEAN = 0.8      # maks. tempo zmiany przechylenia [m/s] (gladki waddle)
OKRES = 2.2           # okres cyklu chodu [s]
DUTY = 0.82           # udzial fazy podporu (pelzanie: jedna konczyna w wymachu)
LICZBA_CYKLI = 4

_LEAN = [0.0]         # wygladzone przechylenie boczne (stan sterownika)

# pelzanie (sekwencja boczna): jedna konczyna w wymachu na raz
FAZA = {'PL': 0.0, 'TP': 0.25, 'PP': 0.5, 'TL': 0.75}

# ----- solver -----
DT = 3.0e-4
SEGMENT = 0.005
CZAS_START = 0.6      # ustabilizowanie postury przed chodem [s]
CZAS = CZAS_START + LICZBA_CYKLI*OKRES + 0.4

# ----- aktuatory / kontakt -----
# umiarkowana sztywnosc + mocna calka (PID) trzyma poze bez sagu; ograniczenie
# momentu (realne dla stawu czlowieka) chroni przed gwaltownymi kopnieciami
K_KUL, C_KUL, KI_KUL, MOM_KUL = 2000.0, 55.0, 1200.0, 300.0   # bark / biodro
K_ZAW, C_ZAW, KI_ZAW, MOM_ZAW = 1400.0, 40.0, 800.0, 220.0    # lokiec / kolano
CALKA_KUL, CALKA_ZAW = 250.0, 180.0
K_KONT, C_KONT, MU_KONT, EPS_KONT = 5.0e4, 500.0, 1.2, 0.003

KOLEJNOSC = (['tulow', 'glowa']
             + [k+'_'+cz for k in KONCZYNY for cz in ('gorny', 'dolny')])
NR = {n: i+1 for i, n in enumerate(KOLEJNOSC)}
N = len(KOLEJNOSC)


def seg(n):
    """Segment anatomiczny dla nazwy czlonu (np. 'PL_gorny' -> 'ramie_L')."""
    if n in S:
        return S[n]
    k, cz = n.rsplit('_', 1)
    return S[KONCZYNY[k][cz]]


def com_lok(s):       # od SM do stawu blizszego (lokalnie, wzdluz z)
    return -s.com*s.dlugosc


def dist_lok(s):      # od SM do stawu dalszego
    return (1 - s.com)*s.dlugosc


def zaczep_lok(k):    # punkt zaczepienia konczyny w ukladzie tulowia
    cfg = KONCZYNY[k]
    z = (1 - COM_TUL)*L_TUL if cfg['przod'] else -COM_TUL*L_TUL
    return wektor(0, cfg['y'], z)


def rodrigues(v, os, kat):
    """Obrot wektora v o kat wokol jednostkowej osi os (wzor Rodriguesa)."""
    c, s = np.cos(kat), np.sin(kat)
    return v*c + np.cross(os, v)*s + os*np.dot(os, v)*(1 - c)


def orient_z_y(zdir, ynorm):
    """Parametry Eulera czlonu, ktorego lokalna os z = zdir (os dluga segmentu),
    a lokalna os y = ynorm (os zawiasu lokcia/kolana, normalna plaszczyzny)."""
    zb = zdir/np.linalg.norm(zdir)
    yb = ynorm - zb*np.dot(zb, ynorm)
    yb = yb/np.linalg.norm(yb)
    xb = np.cross(yb, zb)
    return macierz_na_kwaternion(np.column_stack([xb, yb, zb]))


def _normalna_pionowa(u):
    """Normalna pionowej plaszczyzny zawierajacej kierunek u (do osi zawiasu)."""
    nrm = np.cross(u, np.array([0.0, 0.0, 1.0]))
    if np.linalg.norm(nrm) < 1e-6:
        nrm = np.array([0.0, 1.0, 0.0])
    return nrm/np.linalg.norm(nrm)


def ik_dwuczlon(A, F, L1, L2):
    """IK konczyny 2-czlonowej (przednia, dlon na podlozu): staw kulisty w A,
    zawias (lokiec) w plaszczyznie pionowej, dlon (koniec dolnego) w F."""
    d = np.asarray(F, float) - np.asarray(A, float)
    r = float(np.linalg.norm(d))
    r = min(r, L1 + L2 - 1e-4)
    r = max(r, abs(L1 - L2) + 1e-4)
    dh = d/np.linalg.norm(d)
    nrm = _normalna_pionowa(d)
    a = np.arccos(np.clip((r*r + L1*L1 - L2*L2)/(2*L1*r), -1, 1))
    best = None
    for s in (+1.0, -1.0):
        u1 = rodrigues(dh, nrm, s*a)
        kolano = A + L1*u1
        u2 = F - kolano
        n2 = np.linalg.norm(u2)
        if n2 < 1e-9:
            continue
        if best is None or kolano[2] > best[0]:       # lokiec wyzej
            best = (kolano[2], u1, u2/n2)
    _, u1, u2 = best
    return orient_z_y(u1, nrm), orient_z_y(u2, nrm)


def ik_kolano(A, F, L1, L2):
    """IK konczyny tylnej na KOLANIE: udo (gorny) jako podpora hip->kolano na
    podlozu, podudzie (dolny) podwiniete o staly kat FLEX_KOLANO. F wyznacza
    tylko kierunek poziomy kolana; kolano lezy na okregu (udo siega podlogi),
    wiec udo stoi prawie pionowo i przenosi ciezar osiowo (maly moment)."""
    A = np.asarray(A, float)
    hz = A[2]
    rad = np.sqrt(max(L1*L1 - hz*hz, 1e-4))           # poziomy zasieg kolana
    dxy = np.asarray(F, float)[:2] - A[:2]
    m = np.linalg.norm(dxy)
    dxy = dxy*(rad/m) if m > 1e-6 else np.array([rad, 0.0])
    kolano = np.array([A[0]+dxy[0], A[1]+dxy[1], 0.0])
    u1 = (kolano - A)/L1                               # udo hip->kolano (|.|=L1)
    nrm = _normalna_pionowa(u1)
    u2 = rodrigues(u1, nrm, FLEX_KOLANO)              # podudzie podwiniete
    return orient_z_y(u1, nrm), orient_z_y(u2, nrm)


def ik_konczyna(A, F, L1, L2, dwuczlon=True):
    return (ik_dwuczlon(A, F, L1, L2) if dwuczlon
            else ik_kolano(A, F, L1, L2))


def konfiguracja(r_tul, p_tul, cele):
    """Pelna kinematyka: pozycje SM i orientacje czlonow dla pozy tulowia
    (r_tul, p_tul) i celow stop cele={k: F_swiat}."""
    orient = {'tulow': np.asarray(p_tul, float), 'glowa': np.asarray(p_tul, float)}
    poz = {'tulow': np.asarray(r_tul, float)}
    Rt = R(wektor_p(*p_tul))
    szyja = poz['tulow'] + Rt.dot(wektor(0, 0, (1-COM_TUL)*L_TUL).ravel())
    poz['glowa'] = szyja - Rt.dot(wektor(0, 0, com_lok(S['glowa'])).ravel())
    for k, cfg in KONCZYNY.items():
        A = poz['tulow'] + Rt.dot(zaczep_lok(k).ravel())
        sg, sd = S[cfg['gorny']], S[cfg['dolny']]
        pu, pl = ik_konczyna(A, cele[k], sg.dlugosc, sd.dlugosc, cfg['przod'])
        Ru, Rd = R(wektor_p(*pu)), R(wektor_p(*pl))
        poz[k+'_gorny'] = A + Ru.dot(wektor(0, 0, sg.com*sg.dlugosc).ravel())
        kolano = A + Ru.dot(wektor(0, 0, sg.dlugosc).ravel())
        poz[k+'_dolny'] = kolano + Rd.dot(wektor(0, 0, sd.com*sd.dlugosc).ravel())
        orient[k+'_gorny'], orient[k+'_dolny'] = pu, pl
    return poz, orient


def zloz_q(poz, orient):
    q = np.zeros(7*N)
    for n, nr in NR.items():
        q[3*(nr-1):3*(nr-1)+3] = poz[n]
        q[3*N+4*(nr-1):3*N+4*(nr-1)+4] = orient[n]
    return q


def lean_boczny(faza_glob):
    """Przechylenie ciezaru: strona przeciwna do konczyny wchodzacej w wymach
    (z wyprzedzeniem). Zwraca przesuniecie boczne stop [m] (waddle)."""
    WYPRZEDZENIE = 0.12
    for k in KONCZYNY:
        p = (faza_glob + FAZA[k] + WYPRZEDZENIE) % 1.0
        if p >= DUTY:                                 # ta konczyna wnet w wymachu
            ss = 1.0 if KONCZYNY[k]['y'] > 0 else -1.0
            # stopa przyklejona: przesuniecie celu stopy pcha CIALO w strone
            # przeciwna (reakcja), wiec by przeniesc ciezar na strone podporu
            # (od unoszonej konczyny), przesuwamy stopy w strone konczyny
            return ss*SWAY
    return 0.0


def stopa_cel(k, faza_glob, A, lean=0.0):
    """Cel stopy (swiat): x wzgledem zaczepu A (chodnik ruchomy -> pcha do przodu),
    y bezwzgledne (baza + przechylenie) -> kontrola boczna ciala."""
    p = (faza_glob + FAZA[k]) % 1.0
    if p < DUTY:                                      # podpor: stopa w tyl
        u = p/DUTY
        dx, dz = STRIDE*(0.5 - u), 0.0
    else:                                             # wymach: powrot do przodu
        u = (p - DUTY)/(1 - DUTY)
        dx, dz = STRIDE*(-0.5 + u), KLIRENS*np.sin(np.pi*u)
    y = (BAZA_Y if KONCZYNY[k]['y'] > 0 else -BAZA_Y) + lean
    return np.array([A[0] + dx, y, dz])


def cele_stop(r_tul, p_tul, faza_glob):
    """Cele stop wszystkich konczyn dla biezacej pozy tulowia i fazy chodu."""
    Rt = R(wektor_p(*p_tul))
    cele = {}
    if faza_glob < 0:                                 # poza neutralna (start)
        for k in KONCZYNY:
            A = np.asarray(r_tul, float) + Rt.dot(zaczep_lok(k).ravel())
            y = BAZA_Y if KONCZYNY[k]['y'] > 0 else -BAZA_Y
            cele[k] = np.array([A[0], y, 0.0])        # pod zaczepem, na podlozu
        return cele
    # wygladzenie przechylenia: rampa do celu skokowego (bez impulsow momentu)
    cel_lean = lean_boczny(faza_glob)
    krok = RAMPA_LEAN*SEGMENT
    _LEAN[0] += np.clip(cel_lean - _LEAN[0], -krok, krok)
    for k in KONCZYNY:
        A = np.asarray(r_tul, float) + Rt.dot(zaczep_lok(k).ravel())
        cele[k] = stopa_cel(k, faza_glob, A, _LEAN[0])
    return cele


def zbuduj():
    ukl = Uklad()
    for n in KOLEJNOSC:
        ukl.dodajCzlon(Czlon(NR[n], seg(n).masa, seg(n).tensor))

    aktuatory = {}
    t = NR['tulow']
    # glowa: zawias w szyi
    ukl.dodajWiez(Polaczenie_Obr(t, NR['glowa'],
                                 wektor(0, 0, (1-COM_TUL)*L_TUL),
                                 wektor(0, 0, com_lok(S['glowa'])), *OSIE))
    a_glowa = MomentWzgledny(t, NR['glowa'], wektor(0, 1, 0), wektor(0, 0, 1),
                             300.0, 0.0, 10.0)
    ukl.dodajSileWewn(a_glowa)
    aktuatory['glowa'] = a_glowa

    for k, cfg in KONCZYNY.items():
        g, d = NR[k+'_gorny'], NR[k+'_dolny']
        sg, sd = S[cfg['gorny']], S[cfg['dolny']]
        # bark/biodro: staw kulisty tulow -> gorny
        ukl.dodajWiez(Para_Sferyczna(t, g, zaczep_lok(k),
                                     wektor(0, 0, com_lok(sg))))
        a_kul = MomentSferyczny(t, g, K_KUL, C_KUL, moment_max=MOM_KUL,
                                ki=KI_KUL, calka_max=CALKA_KUL)
        ukl.dodajSileWewn(a_kul)
        # lokiec/kolano: zawias gorny -> dolny
        ukl.dodajWiez(Polaczenie_Obr(g, d, wektor(0, 0, dist_lok(sg)),
                                     wektor(0, 0, com_lok(sd)), *OSIE))
        a_zaw = MomentWzgledny(g, d, wektor(0, 1, 0), wektor(0, 0, 1),
                               K_ZAW, 0.0, C_ZAW, moment_max=MOM_ZAW,
                               ki=KI_ZAW, calka_max=CALKA_ZAW)
        ukl.dodajSileWewn(a_zaw)
        # kontakt: przednie na dloni (koniec przedramienia), tylne na kolanie
        # (koniec uda) - podudzie podwiniete nad podlogą
        if cfg['przod']:
            ukl.dodajSileWewn(SilaKontaktu(d, wektor(0, 0, dist_lok(sd)),
                                           k=K_KONT, c=C_KONT, mu=MU_KONT, eps=EPS_KONT))
        else:
            ukl.dodajSileWewn(SilaKontaktu(g, wektor(0, 0, dist_lok(sg)),
                                           k=K_KONT, c=C_KONT, mu=MU_KONT, eps=EPS_KONT))
        aktuatory[k] = (a_kul, a_zaw)

    ukl.grawitacja = True
    return ukl, aktuatory


def ustaw_cele(aktuatory, q, faza_glob):
    """Ustawia cele aktuatorow wzgledem tulowia ODNIESIENIA (idealna wysokosc
    H_TULOW i poziom), nie wzgledem aktualnego. Dzieki temu barki serwowane sa
    do orientacji bezwzglednej: gdy tulow opada albo przechyla sie, konczyny
    odpychaja go z powrotem do pozy odniesienia (regulacja wysokosci/postawy
    przez kontakt). Postep do przodu przez sledzenie x tulowia."""
    p_akt = q[3*N:3*N+4]                               # aktualna orientacja tulowia
    r_ref = np.array([q[0], 0.0, H_TULOW])            # tulow odniesienia
    cele = cele_stop(r_ref, P_TUL, faza_glob)
    _, orient = konfiguracja(r_ref, P_TUL, cele)
    q_cel = q.copy()
    for n in KOLEJNOSC:
        q_cel[3*N+4*(NR[n]-1):3*N+4*(NR[n]-1)+4] = orient[n]

    aktuatory['glowa'].theta_cel = 0.0
    for k in KONCZYNY:
        a_kul, a_zaw = aktuatory[k]
        # bark/biodro: orientacja BEZWZGLEDNA gornego (serwo postawy)
        a_kul.p_cel = mnoz_kwaterniony(sprzezenie_kwaternionu(p_akt),
                                       orient[k+'_gorny'])
        # lokiec/kolano: kat z konfiguracji odniesienia
        a_zaw.theta_cel = a_zaw.kat(q_cel, N)


def symuluj():
    ukl, aktuatory = zbuduj()
    _LEAN[0] = 0.0
    cele0 = cele_stop([0.0, 0.0, H_TULOW], P_TUL, -1.0)
    poz0, orient0 = konfiguracja([0.0, 0.0, H_TULOW], P_TUL, cele0)
    q = ukl.projekcja_polozen(zloz_q(poz0, orient0))
    dq = np.zeros(7*N)
    ustaw_cele(aktuatory, q, -1.0)

    klatki = []
    t = 0.0
    n_seg = int(CZAS/SEGMENT)
    for _ in range(n_seg):
        faza = (t - CZAS_START)/OKRES if t >= CZAS_START else -1.0
        ustaw_cele(aktuatory, q, faza)
        ukl.sym2(np.concatenate((q, dq)), 0.0, SEGMENT, DT)
        Y = ukl.Y
        for w in Y[:-1]:
            klatki.append(w[0:7*N].copy())
        q = Y[-1][0:7*N].copy()
        dq = Y[-1][7*N:14*N].copy()
        t += SEGMENT
        if not np.all(np.isfinite(q)):
            print(f'NaN w t={t:.3f}')
            break
    return ukl, klatki


def eksportuj(klatki, co_ile=36, plik='web/dane_czworaka.js'):
    dane_klatki = []
    for q in klatki[::co_ile]:
        czlony = []
        for i in range(N):
            r = q[3*i:3*i+3]
            p = q[3*N+4*i:3*N+4*i+4]
            czlony.append({'r': [round(float(v), 4) for v in r],
                           'p': [round(float(v), 5) for v in p]})
        dane_klatki.append(czlony)

    segmenty_op = []
    for n in KOLEJNOSC:
        s = seg(n)
        segmenty_op.append({'nazwa': n, 'dlugosc': round(s.dlugosc, 4),
                            'com': s.com, 'promien': round(s.promien, 4)})

    dane = {'dt': DT*co_ile, 'segmenty': segmenty_op, 'klatki': dane_klatki}
    katalog = os.path.dirname(os.path.abspath(plik))
    os.makedirs(katalog, exist_ok=True)
    with open(plik, 'w') as f:
        f.write('const DANE = ' + json.dumps(dane) + ';\n')
    print(f'zapisano {len(dane_klatki)} klatek do {plik}')


if __name__ == '__main__':
    ukl, klatki = symuluj()
    if klatki:
        x0, xk = klatki[0][0], klatki[-1][0]
        zt = [q[2] for q in klatki]
        print(f'tulow: start x={x0:.3f}, koniec x={xk:.3f} '
              f'(do przodu {xk-x0:+.3f} m przez {LICZBA_CYKLI} cykli)')
        print(f'wysokosc tulowia: {min(zt):.3f}..{max(zt):.3f} m (nom {H_TULOW})')
        sciezka = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               '..', 'web', 'dane_czworaka.js')
        eksportuj(klatki, plik=os.path.normpath(sciezka))
