# -*- coding: utf-8 -*-
# Przyklad: najprostszy robot kroczacy (chod cyrklowy, compass gait).
#
# Dwie nogi (prety) polaczone przegubem biodrowym. Noga podporowa jest
# przypieta do podloza przegubem obrotowym (os y); robot przewaza sie
# nad nia jak odwrocone wahadlo, a najprostszy silownik biodra (stala
# para momentow: +tau na noge podporowa, -tau na wymachowa) przerzuca
# noge wymachowa do przodu i uzupelnia energie tracona w uderzeniach.
# Gdy noga podporowa przejdzie do przodu o kat kroku, stopa nogi
# wymachowej laduje: podpora przechodzi na druga noge.
#
# Zmiana podpory wykorzystuje nowe metody biblioteki:
#  - projekcja_polozen: osadzenie stopy dokladnie na wiezach nowej podpory,
#  - projekcja_predkosci: uderzenie plastyczne (predkosci po lądowaniu).
#
# Wynik: web/dane_robot.js do wizualizacji Three.js (web/robot.html).

import json
import os

import numpy as np

from uw_dyn import (Uklad, Czlon, Polaczenie_Obr, SilaZewn, SilaWewnProst,
                    wektor, u2p, R, G, wektor_p)

# ----- parametry robota -----
L = 1.0          # dlugosc nogi [m]
MASA = 5.0       # masa nogi [kg]
ALFA = 0.25      # pol kata miedzy nogami przy ladowaniu [rad]
COM = 0.8        # polozenie srodka masy nogi (ulamek dlugosci od stopy;
                 # masa skupiona przy biodrze jak w chodziarzach McGeera)
OMEGA0 = 1.8      # poczatkowa predkosc katowa (pchniecie startowe) [rad/s]
TAU_POD = 9.0     # moment na nodze podporowej (napedza przewazanie) [N*m]
TAU_WYM = 10.0    # moment na nodze wymachowej (przerzuca ja do przodu) [N*m]
C_BIODRA = 3.0    # tlumik miedzy nogami (gasi wymach przy celu) [N*s/m]
DT = 0.002
MAKS_KROKOW = 8
OS_Y = np.array([0.0, 1.0, 0.0])


def tensor_preta(m, L, promien=0.04):
    Jxx = m * (3 * promien ** 2 + L ** 2) / 12
    Jzz = m * promien ** 2 / 2
    return np.diag([Jxx, Jxx, Jzz])


def kat_czlonu(q, nr):
    """Kat absolutny czlonu wokol osi y (obrot czysto wokol y)."""
    p = q[6 + 4 * (nr - 1):6 + 4 * (nr - 1) + 4]
    return 2 * np.arctan2(p[2], p[0])


def stopa_czlonu(q, nr):
    """Polozenie stopy (dolnego konca) nogi nr w konfiguracji q."""
    r = q[3 * (nr - 1):3 * (nr - 1) + 3].reshape(3, 1)
    p = q[6 + 4 * (nr - 1):6 + 4 * (nr - 1) + 4]
    return (r + R(wektor_p(*p)).dot(wektor(0, 0, -COM * L))).ravel()


def zbuduj_faze(noga_podporowa, stopa_xyz):
    """Uklad dla fazy z podpora na wskazanej nodze przypieta w stopa_xyz."""
    ukl = Uklad()
    ukl.dodajCzlon(Czlon(1, MASA, tensor_preta(MASA, L)))
    ukl.dodajCzlon(Czlon(2, MASA, tensor_preta(MASA, L)))

    OSIE = (wektor(1, 0, 0), wektor(0, 0, 1), wektor(0, 1, 0))
    # przegub podporowy: podloze - stopa nogi podporowej
    ukl.dodajWiez(Polaczenie_Obr(0, noga_podporowa,
                                 wektor(*stopa_xyz), wektor(0, 0, -COM * L), *OSIE))
    # biodro: gorne konce obu nog
    ukl.dodajWiez(Polaczenie_Obr(1, 2, wektor(0, 0, (1 - COM) * L),
                                 wektor(0, 0, (1 - COM) * L), *OSIE))

    # silownik biodra: stale momenty wokol osi y
    # (+tau_pod napedza przewazanie, -tau_wym przerzuca wymach do przodu)
    noga_wymachowa = 2 if noga_podporowa == 1 else 1
    ukl.dodajSileZewn(SilaZewn(noga_podporowa, 'ny', TAU_POD))
    ukl.dodajSileZewn(SilaZewn(noga_wymachowa, 'ny', -TAU_WYM))
    # tlumik miedzy srodkami nog: gasi predkosc wzgledna wymachu,
    # zeby noga nie przekrecala sie pod robotem
    ukl.dodajSileWewn(SilaWewnProst(1, 2, wektor(0, 0, 0), wektor(0, 0, 0),
                                    0, 0, C_BIODRA, 0))
    ukl.grawitacja = True
    return ukl


def q_startowe(stopa_x, kat1, kat2):
    """Konfiguracja: noga 1 podporowa (stopa w stopa_x), katy absolutne nog."""
    R1 = R(wektor_p(*u2p(OS_Y, kat1)))
    R2 = R(wektor_p(*u2p(OS_Y, kat2)))
    stopa = wektor(stopa_x, 0, 0)
    biodro = stopa + R1.dot(wektor(0, 0, L))
    r1 = stopa + R1.dot(wektor(0, 0, COM * L))
    r2 = biodro - R2.dot(wektor(0, 0, (1 - COM) * L))
    q = np.zeros(14)
    q[0:3], q[3:6] = r1.ravel(), r2.ravel()
    q[6:10] = u2p(OS_Y, kat1)
    q[10:14] = u2p(OS_Y, kat2)
    return q


def dq_obrotu(q, omega):
    """Predkosci odpowiadajace sztywnemu obrotowi calosci wokol stopy nogi 1."""
    stopa = stopa_czlonu(q, 1)
    dq = np.zeros(14)
    for nr in (1, 2):
        r = q[3 * (nr - 1):3 * (nr - 1) + 3]
        wzgl = r - stopa
        dq[3 * (nr - 1):3 * (nr - 1) + 3] = omega * np.array([wzgl[2], 0.0, -wzgl[0]])
        p = wektor_p(*q[6 + 4 * (nr - 1):6 + 4 * (nr - 1) + 4])
        dq[6 + 4 * (nr - 1):6 + 4 * (nr - 1) + 4] = \
            (0.5 * G(p).T.dot(wektor(0, omega, 0))).ravel()
    return dq


def symuluj_chod():
    """Kolejne kroki chodu. Zwraca (klatki, opisy_krokow)."""
    stopa_x = 0.0
    noga_podporowa = 1
    q = q_startowe(stopa_x, -ALFA, ALFA)
    dq = dq_obrotu(q, OMEGA0)

    klatki = []      # (q, noga_podporowa) dla kazdego zapisu
    kroki = []

    for krok in range(MAKS_KROKOW):
        ukl = zbuduj_faze(noga_podporowa, (stopa_x, 0, 0))
        # osadzenie stanu na wiezach fazy (przy zmianie podpory: uderzenie)
        q = ukl.projekcja_polozen(q)
        dq = ukl.projekcja_predkosci(q, dq)

        ukl.sym2(np.concatenate((q, dq)), 0.0, 2.5, DT)
        Y = ukl.Y

        # zdarzenie ladowania (jak w klasycznym chodzie cyrklowym): kat
        # podporowej + kat wymachowej przecina zero przy podporowej z przodu,
        # czyli nogi sa symetryczne wzgledem pionu, a stopa wymachowa jest
        # na podlozu w odleglosci 2L sin(kat podporowej) przed podpora
        # (szuranie stopy w polowie wymachu jest ignorowane, jak w modelu
        # McGeera)
        noga_wymachowa = 2 if noga_podporowa == 1 else 1
        koniec = None
        for i in range(1, len(Y)):
            q_i = Y[i][0:14]
            kat_pod = kat_czlonu(q_i, noga_podporowa)
            kat_wym = kat_czlonu(q_i, noga_wymachowa)
            if kat_pod >= 0.02 and kat_pod + kat_wym <= 0.0:
                koniec = i
                break
        if koniec is None:
            kroki.append({'krok': krok, 'status': 'zatrzymanie'})
            klatki.extend((Y[i][0:14], noga_podporowa) for i in range(len(Y)))
            break

        klatki.extend((Y[i][0:14], noga_podporowa) for i in range(koniec))

        q = Y[koniec][0:14].copy()
        dq = Y[koniec][14:28].copy()

        # nowa podpora: stopa nogi wymachowej (dociagnieta do podloza)
        nowa_stopa = stopa_czlonu(q, noga_wymachowa)
        dlugosc_kroku = nowa_stopa[0] - stopa_x
        kroki.append({'krok': krok, 'status': 'ok',
                      'dlugosc': round(float(dlugosc_kroku), 3),
                      'stopa_x': round(float(nowa_stopa[0]), 3)})
        stopa_x = float(nowa_stopa[0])
        noga_podporowa = noga_wymachowa

    return klatki, kroki


def eksportuj(klatki, co_ile=10, plik='web/dane_robot.js'):
    dane_klatki = []
    for q, podpora in klatki[::co_ile]:
        czlony = []
        for k in range(2):
            r = q[3 * k:3 * k + 3]
            p = q[6 + 4 * k:6 + 4 * k + 4]
            czlony.append({'r': [round(v, 5) for v in r],
                           'p': [round(v, 6) for v in p]})
        dane_klatki.append({'czlony': czlony, 'podpora': podpora})

    dane = {'dt': DT * co_ile, 'dlugosc_nogi': L, 'com': COM,
            'klatki': dane_klatki}
    katalog = os.path.dirname(os.path.abspath(plik))
    os.makedirs(katalog, exist_ok=True)
    with open(plik, 'w') as f:
        f.write('const DANE = ' + json.dumps(dane) + ';\n')
    print(f'zapisano {len(dane_klatki)} klatek do {plik}')


if __name__ == '__main__':
    klatki, kroki = symuluj_chod()
    print('kroki:')
    for k in kroki:
        print('  ', k)
    udane = [k for k in kroki if k['status'] == 'ok']
    print(f'{len(udane)} udanych krokow, dystans: '
          f'{udane[-1]["stopa_x"] if udane else 0} m')

    sciezka = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           '..', 'web', 'dane_robot.js')
    eksportuj(klatki, co_ile=10, plik=os.path.normpath(sciezka))
