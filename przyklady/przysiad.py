# -*- coding: utf-8 -*-
# Przyklad: staw kolanowy podczas przysiadu.
#
# Model czlowieka w plaszczyznie strzalkowej (x: przod, z: gora):
# trzy czlony (podudzie, udo, tulow) polaczone przegubami obrotowymi
# o osi y (staw skokowy zakotwiczony w podlozu, staw kolanowy,
# staw biodrowy). Miesnie zamodelowane jako elementy sprezysto-tlumiace
# (SilaWewnProst) o dlugosciach swobodnych dobranych do pozycji przysiadu:
# postac startuje ze stania i "miesnie" sciagaja ja do przysiadu.
#
# Wynik: web/dane_przysiad.js (pozycje i kwaterniony klatka po klatce)
# do wizualizacji w Three.js (web/przysiad.html).

import json
import os

import numpy as np

from uw_dyn import (Uklad, Czlon, Polaczenie_Obr, SilaWewnProst,
                    wektor, u2p, R, wektor_p)

# ----- antropometria (obie nogi zlaczone w jeden lancuch) -----
L1, M1 = 0.45, 7.0    # podudzie: dlugosc [m], masa [kg]
L2, M2 = 0.45, 16.0   # udo
L3, M3 = 0.80, 45.0   # tulow z glowa i ramionami

# tensory bezwladnosci jak dla preta wokol srodka (os podluzna z)
def tensor_preta(m, L, promien=0.06):
    Jxx = m * (3 * promien ** 2 + L ** 2) / 12
    Jzz = m * promien ** 2 / 2
    return np.diag([Jxx, Jxx, Jzz])

# ----- pozycja docelowa przysiadu (katy absolutne wokol osi y) -----
# dodatni kat pochyla wierzcholek czlonu do przodu (+x)
KAT_PODUDZIE = 0.44   # ~25 stopni do przodu
KAT_UDO = -1.31       # ~-75 stopni (kolano z przodu, biodro z tylu)
KAT_TULOW = 0.52      # ~30 stopni pochylenia tulowia

OS_Y = np.array([0.0, 1.0, 0.0])


def kinematyka(katy):
    """Polozenia srodkow mas i kwaterniony czlonow dla katow (t1,t2,t3)."""
    t1, t2, t3 = katy
    R1 = R(wektor_p(*u2p(OS_Y, t1)))
    R2 = R(wektor_p(*u2p(OS_Y, t2)))
    R3 = R(wektor_p(*u2p(OS_Y, t3)))

    kostka = np.zeros((3, 1))
    kolano = kostka + R1.dot(wektor(0, 0, L1))
    biodro = kolano + R2.dot(wektor(0, 0, L2))

    r1 = kostka + R1.dot(wektor(0, 0, L1 / 2))
    r2 = kolano + R2.dot(wektor(0, 0, L2 / 2))
    r3 = biodro + R3.dot(wektor(0, 0, L3 / 2))

    p1 = u2p(OS_Y, t1)
    p2 = u2p(OS_Y, t2)
    p3 = u2p(OS_Y, t3)
    return (r1, r2, r3), (p1, p2, p3)


def q_z_katow(katy):
    """Wektor wspolrzednych q (3N pozycji + 4N kwaternionow) z katow."""
    (r1, r2, r3), (p1, p2, p3) = kinematyka(katy)
    q = np.zeros(7 * 3)
    q[0:3], q[3:6], q[6:9] = r1.ravel(), r2.ravel(), r3.ravel()
    q[9:13], q[13:17], q[17:21] = p1, p2, p3
    return q


def punkt_na_czlonie(nr_czlonu, s_lokalny, q):
    """Punkt zaczepu (wsp. globalne) dla czlonu w konfiguracji q."""
    N = 3
    r = q[3 * (nr_czlonu - 1):3 * (nr_czlonu - 1) + 3].reshape(3, 1)
    p = q[3 * N + 4 * (nr_czlonu - 1):3 * N + 4 * (nr_czlonu - 1) + 4]
    return r + R(wektor_p(*p)).dot(s_lokalny)


def zbuduj_model():
    ukl = Uklad()
    ukl.dodajCzlon(Czlon(1, M1, tensor_preta(M1, L1)))
    ukl.dodajCzlon(Czlon(2, M2, tensor_preta(M2, L2)))
    ukl.dodajCzlon(Czlon(3, M3, tensor_preta(M3, L3)))

    OSIE = (wektor(1, 0, 0), wektor(0, 0, 1), wektor(0, 1, 0))
    # staw skokowy: podloze-podudzie, w poczatku ukladu
    ukl.dodajWiez(Polaczenie_Obr(0, 1, wektor(0, 0, 0), wektor(0, 0, -L1 / 2), *OSIE))
    # staw kolanowy: podudzie-udo
    ukl.dodajWiez(Polaczenie_Obr(1, 2, wektor(0, 0, L1 / 2), wektor(0, 0, -L2 / 2), *OSIE))
    # staw biodrowy: udo-tulow
    ukl.dodajWiez(Polaczenie_Obr(2, 3, wektor(0, 0, L2 / 2), wektor(0, 0, -L3 / 2), *OSIE))

    ukl.grawitacja = True
    return ukl


# zaczepy "miesni": jeden element na staw, rozpiety miedzy srodkami
# sasiednich segmentow (odleglosc srodkow jest monotoniczna w kacie stawu,
# wiec sprezyna stabilizuje staw w obu kierunkach);
# format: (czlon_i, s_i, czlon_j, s_j)
MIESNIE = [
    # staw skokowy: podloze (za pieta) - srodek podudzia
    (0, wektor(-0.35, 0, 0.0), 1, wektor(0, 0, 0)),
    # staw kolanowy: srodek podudzia - srodek uda
    (1, wektor(0, 0, 0), 2, wektor(0, 0, 0)),
    # staw biodrowy: srodek uda - srodek tulowia
    (2, wektor(0, 0, 0), 3, wektor(0, 0, 0)),
]

SZTYWNOSC = [1.2e5, 1.0e5, 1.0e5]   # k [N/m]
TLUMIENIE = [1.2e4, 1.0e4, 1.0e4]   # c [N*s/m]


def dlugosci_w_pozycji(katy):
    """Dlugosci miesni w zadanej pozycji (punkt startowy strojenia)."""
    q = q_z_katow(katy)
    dl = []
    for (ci, si, cj, sj) in MIESNIE:
        pA = si if ci == 0 else punkt_na_czlonie(ci, si, q)
        pB = punkt_na_czlonie(cj, sj, q)
        dl.append(float(np.linalg.norm(pB - pA)))
    return np.array(dl)


def zbuduj_uklad_z_miesniami(l0_lista):
    """Uklad z miesniami o zadanych dlugosciach swobodnych."""
    ukl = zbuduj_model()
    for (ci, si, cj, sj), k, c, l0 in zip(MIESNIE, SZTYWNOSC, TLUMIENIE, l0_lista):
        ukl.dodajSileWewn(SilaWewnProst(ci, cj, si, sj, k, float(l0), c, 0))
    return ukl


def katy_z_q(q):
    """Katy absolutne czlonow wokol osi y z kwaternionow (obrot czysto wokol y)."""
    N = 3
    katy = []
    for k in range(N):
        p = q[3 * N + 4 * k:3 * N + 4 * k + 4]
        katy.append(2 * np.arctan2(p[2], p[0]))
    return np.array(katy)


GRAW = 9.80665


def _energia_potencjalna(katy):
    (r1, r2, r3), _ = kinematyka(katy)
    return GRAW * (M1 * r1[2, 0] + M2 * r2[2, 0] + M3 * r3[2, 0])


def _gradient(fun, katy, eps=1e-6):
    """Gradient funkcji skalarnej lub wektorowej po katach (roznice centralne)."""
    kolumny = []
    for j in range(3):
        kp = np.array(katy, dtype=float)
        km = kp.copy()
        kp[j] += eps
        km[j] -= eps
        kolumny.append((np.asarray(fun(kp)) - np.asarray(fun(km))) / (2 * eps))
    return np.column_stack(kolumny)


def l0_rownowagi(katy_docelowe):
    """Dlugosci swobodne dajace dokladna rownowage statyczna w pozycji docelowej.

    Warunek rownowagi: dV_graw/dkat + A f = 0, gdzie A[j,m] = dl_m/dkat_j
    (ramiona dzialania miesni), f = k (l - l0) to sily miesni.
    Stad f = -A^-1 grad(V), a l0 = l_cel - f/k."""
    l_cel = dlugosci_w_pozycji(katy_docelowe)
    grad_V = _gradient(_energia_potencjalna, katy_docelowe).ravel()
    A = _gradient(dlugosci_w_pozycji, katy_docelowe).T  # A[j,m] = dl_m/dkat_j
    f = np.linalg.solve(A, -grad_V)
    l0 = l_cel - f / np.array(SZTYWNOSC)

    # kontrola statecznosci: hesjan energii calkowitej musi byc dodatnio okreslony
    def V_calk(katy):
        dl = dlugosci_w_pozycji(katy)
        return (_energia_potencjalna(katy)
                + 0.5 * np.sum(np.array(SZTYWNOSC) * (dl - l0) ** 2))
    H = _gradient(lambda k: _gradient(V_calk, k).ravel(), katy_docelowe, eps=1e-4)
    wart_wlasne = np.linalg.eigvalsh(0.5 * (H + H.T))
    if wart_wlasne.min() <= 0:
        raise RuntimeError(f'pozycja docelowa niestateczna: {wart_wlasne}')

    return l0


def symuluj(t_kon=5.0, dt=0.001, l0_lista=None):
    if l0_lista is None:
        l0_lista = l0_rownowagi((KAT_PODUDZIE, KAT_UDO, KAT_TULOW))
    ukl = zbuduj_uklad_z_miesniami(l0_lista)
    # start: stanie z lekkim ugieciem inicjujacym ruch w dol
    q0 = q_z_katow((0.05, -0.10, 0.05))
    y0 = np.concatenate((q0, np.zeros(21)))
    ukl.sym2(y0, 0.0, t_kon, dt, 5, 5)
    return ukl


def eksportuj(ukl, dt, co_ile=20, plik='web/dane_przysiad.js'):
    """Zapis klatek animacji jako plik JS (const DANE = {...})."""
    N = 3
    klatki = []
    for wiersz in ukl.Y[::co_ile]:
        q = wiersz[0:7 * N]
        czlony = []
        for k in range(N):
            r = q[3 * k:3 * k + 3]
            p = q[3 * N + 4 * k:3 * N + 4 * k + 4]
            czlony.append({'r': [round(v, 5) for v in r],
                           'p': [round(v, 6) for v in p]})
        klatki.append(czlony)

    dane = {
        'dt': dt * co_ile,
        'dlugosci': [L1, L2, L3],
        'nazwy': ['podudzie', 'udo', 'tulow'],
        'klatki': klatki,
    }
    katalog = os.path.dirname(os.path.abspath(plik))
    os.makedirs(katalog, exist_ok=True)
    with open(plik, 'w') as f:
        f.write('const DANE = ' + json.dumps(dane) + ';\n')
    print(f'zapisano {len(klatki)} klatek do {plik}')


if __name__ == '__main__':
    dt = 0.001
    ukl = symuluj(t_kon=3.0, dt=dt)

    katy0 = np.degrees(katy_z_q(ukl.Y[0][0:21]))
    katyK = np.degrees(katy_z_q(ukl.Y[-1][0:21]))
    cel = np.degrees([KAT_PODUDZIE, KAT_UDO, KAT_TULOW])
    print('katy [stopnie]   start     koniec    cel')
    for n, a, b, c in zip(['podudzie', 'udo     ', 'tulow   '], katy0, katyK, cel):
        print(f'  {n}    {a:8.1f}  {b:8.1f}  {c:8.1f}')
    zgiecie_kolana = katyK[0] - katyK[1]
    print(f'zgiecie kolana na koncu: {zgiecie_kolana:.1f} stopni')

    sciezka = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           '..', 'web', 'dane_przysiad.js')
    eksportuj(ukl, dt, co_ile=20, plik=os.path.normpath(sciezka))
