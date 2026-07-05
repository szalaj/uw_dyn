# -*- coding: utf-8 -*-
# Przyklad dla projektu logistyka: urealnienie kosztu transportu.
#
# W grze "Zasoby" dowoz kosztuje `d * koszt_transportu + w * paliwo_za_wysokosc`
# (w = poziom wysokosci zony, 0..3). Ten przyklad wyprowadza te wspolczynniki
# z fizyki jazdy po terenie o rosnacej szorstkosci i nachyleniu:
#
#  - model pol-samochodu (nadwozie na dwoch resorach z tlumikami), jak
#    w klasycznej analizie zawieszen: pojazd "stoi", a profil drogi
#    przesuwa sie pod nim z predkoscia jazdy (wymuszenie kinematyczne
#    pionowe pod kolami),
#  - straty paliwa ~ energia rozproszona w amortyzatorach (nierownosci)
#    + praca wspinaczki m*g*(suma podjazdow),
#  - rms przyspieszen pionowych nadwozia ~ ryzyko uszkodzenia ladunku.
#
# Biblioteka nie ma wiezow zaleznych od czasu, wiec jazda jest symulowana
# odcinkami: co segment kotwice resorow dostaja nowa wysokosc terenu
# (zero-order hold), a stan (q, dq) przechodzi dalej.
#
# Wynik: tabela wspolczynnikow + web/dane_transport.js (animacja przejazdu
# po terenie poziomu 3).

import json
import os

import numpy as np

from uw_dyn import Uklad, Czlon, SilaWewnProst, wektor
from uw_dyn.uklad import GRAWITACJA

# ----- pojazd (pol-samochod) -----
MASA = 800.0          # nadwozie [kg]
DLUGOSC = 3.0         # dlugosc nadwozia [m]
ROZSTAW = 2.4         # rozstaw resorow [m]
K_RESORU = 8.0e4      # sztywnosc resoru [N/m]
C_RESORU = 8.0e3      # tlumienie amortyzatora [N*s/m]
L0_RESORU = 0.6       # dlugosc swobodna resoru [m]
PREDKOSC = 8.0        # predkosc jazdy [m/s]
DYSTANS = 100.0       # dlugosc przejazdu (1 heks) [m]
DT = 0.002
SEGMENT = 0.004       # czas segmentu (aktualizacja wysokosci kotwic) [s]

J_NADWOZIA = np.diag([MASA * (1.0**2 + 1.2**2) / 12,
                      MASA * (DLUGOSC**2 + 1.2**2) / 12,
                      MASA * (DLUGOSC**2 + 1.0**2) / 12])

# ----- profile terenu odpowiadajace poziomom wysokosci gry -----
POZIOMY = {
    0: {'nazwa': 'rownina', 'amplituda': 0.03, 'nachylenie': 0.00},
    1: {'nazwa': 'pagorki', 'amplituda': 0.12, 'nachylenie': 0.03},
    2: {'nazwa': 'wyzyna',  'amplituda': 0.25, 'nachylenie': 0.06},
    3: {'nazwa': 'gory',    'amplituda': 0.40, 'nachylenie': 0.10},
}


def teren(x, poziom):
    """Wysokosc profilu drogi w punkcie x (bez skladnika nachylenia,
    ktory rozliczany jest analitycznie jako praca wspinaczki)."""
    a = POZIOMY[poziom]['amplituda']
    return (a * np.sin(2 * np.pi * x / 11.0)
            + 0.6 * a * np.sin(2 * np.pi * x / 4.7 + 1.3)
            + 0.3 * a * np.sin(2 * np.pi * x / 2.3 + 2.1))


# statyczne ugiecie resoru pod polowa masy nadwozia
UGIECIE = MASA * GRAWITACJA / (2 * K_RESORU)
# zaczepy resorow w ukladzie nadwozia (naroznik dolny)
ZACZEPY = (ROZSTAW / 2, -ROZSTAW / 2)


def zbuduj_pojazd(x_jazdy, poziom):
    """Uklad: nadwozie na dwoch pionowych resorach; kotwice na wysokosci
    terenu pod kolami dla biezacego polozenia x_jazdy."""
    ukl = Uklad()
    ukl.dodajCzlon(Czlon(1, MASA, J_NADWOZIA))
    for poz in ZACZEPY:
        z_terenu = teren(x_jazdy + poz, poziom)
        ukl.dodajSileWewn(SilaWewnProst(
            0, 1, wektor(poz, 0, z_terenu), wektor(poz, 0, -0.3),
            K_RESORU, L0_RESORU, C_RESORU, 0))
    ukl.grawitacja = True
    return ukl


def stan_poczatkowy():
    """Nadwozie w rownowadze statycznej na plaskim odcinku startowym."""
    q = np.zeros(7)
    q[2] = 0.3 + L0_RESORU + UGIECIE
    q[3] = 1.0
    return np.concatenate((q, np.zeros(7)))


def przejazd(poziom, zapis_klatek=False):
    """Symulacja przejazdu; zwraca metryki energetyczne (i klatki)."""
    y = stan_poczatkowy()
    q, dq = y[0:7], y[7:14]
    klatki = []
    E_tlumikow = 0.0
    rms_akum, rms_n = 0.0, 0
    dz_poprz = 0.0

    n_seg = int(DYSTANS / PREDKOSC / SEGMENT)
    for seg in range(n_seg):
        x_jazdy = seg * SEGMENT * PREDKOSC
        ukl = zbuduj_pojazd(x_jazdy, poziom)
        ukl.sym2(np.concatenate((q, dq)), 0.0, SEGMENT, DT)
        Y = ukl.Y

        # straty w tlumikach: c * (dl/dt)^2 sumowane po krokach
        for s in ukl.silyWewn:
            dl_poprz = None
            for w in Y:
                l = s.dlugosc(w[0:7], 1)
                if dl_poprz is not None:
                    E_tlumikow += s.c * ((l - dl_poprz) / DT) ** 2 * DT
                dl_poprz = l

        # rms pionowego przyspieszenia nadwozia (komfort ladunku)
        if seg > n_seg // 10:  # pomin poczatkowy stan nieustalony
            a_z = (Y[-1][9] - dz_poprz) / SEGMENT
            rms_akum += a_z ** 2
            rms_n += 1
        dz_poprz = Y[-1][9]

        if zapis_klatek:
            for w in Y[:-1]:
                klatki.append({'x': round(x_jazdy, 3),
                               'z': round(float(w[2]), 4),
                               'p': [round(v, 5) for v in w[3:7]]})
        q = Y[-1][0:7].copy()
        dq = Y[-1][7:14].copy()

    # praca wspinaczki: nachylenie + dodatnie przyrosty falek
    xx = np.linspace(0, DYSTANS, 1000)
    zz = teren(xx, poziom) + POZIOMY[poziom]['nachylenie'] * xx
    E_wspinaczki = MASA * GRAWITACJA * float(np.sum(np.maximum(np.diff(zz), 0)))

    return {
        'poziom': poziom,
        'nazwa': POZIOMY[poziom]['nazwa'],
        'E_tlumikow_kJ': round(E_tlumikow / 1000, 1),
        'E_wspinaczki_kJ': round(E_wspinaczki / 1000, 1),
        'E_razem_kJ': round((E_tlumikow + E_wspinaczki) / 1000, 1),
        'rms_przysp_ms2': round(float(np.sqrt(rms_akum / max(rms_n, 1))), 2),
    }, klatki


if __name__ == '__main__':
    wyniki = []
    klatki_gory = []
    for w in range(4):
        wynik, klatki = przejazd(w, zapis_klatek=(w == 3))
        wyniki.append(wynik)
        if w == 3:
            klatki_gory = klatki
        print(wynik)

    E0 = wyniki[0]['E_razem_kJ']
    print('\npoziom  nazwa     E [kJ/heks]  mnoznik paliwa')
    for wnk in wyniki:
        print(f"  {wnk['poziom']}    {wnk['nazwa']:9s} {wnk['E_razem_kJ']:10.1f}  "
              f"{wnk['E_razem_kJ'] / E0:8.2f}")

    # dopasowanie liniowe: mnoznik(w) ~ 1 + w * paliwo_za_wysokosc
    mnozniki = np.array([wnk['E_razem_kJ'] / E0 for wnk in wyniki])
    wsp = float(np.polyfit(np.arange(4), mnozniki, 1)[0])
    print("\nregula gry: koszt = d * koszt_transportu + w * paliwo_za_wysokosc")
    print(f"fizyczne 'paliwo za wysokosc' = {wsp:.2f} x koszt na rowninie")

    dane = {
        'wyniki': wyniki,
        'paliwo_za_wysokosc': round(wsp, 2),
        'dt': DT,
        'predkosc': PREDKOSC,
        'teren': [[round(float(x), 2), round(float(teren(x, 3)), 3)]
                  for x in np.linspace(0, DYSTANS, 500)],
        'pojazd': {'dlugosc': DLUGOSC, 'rozstaw': ROZSTAW},
        'klatki': klatki_gory[::10],
    }
    sciezka = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           '..', 'web', 'dane_transport.js')
    with open(os.path.normpath(sciezka), 'w') as f:
        f.write('const DANE = ' + json.dumps(dane) + ';\n')
    print(f"\nzapisano {len(dane['klatki'])} klatek do web/dane_transport.js")
