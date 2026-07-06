# -*- coding: utf-8 -*-
# testy integratora sym3 (RATTLE: Verlet predkosciowy z wiezami w kroku)
#  - rzad 2 (blad maleje jak dt^2; sym2 ma rzad 1)
#  - wiezy na poziomie maszynowym
#  - energia na ukladach z wiezami o rzedy wielkosci dokladniejsza niz sym2

import numpy as np
import pytest

from uw_dyn import (Uklad, Czlon, Polaczenie_Obr, SilaWewnProst, wektor, u2p,
                    R, G, wektor_p)


def _bryla(met, dt):
    """Bryla swobodna (os posrednia): wzgledny blad kretu w ukladzie swiata."""
    Jb = np.diag([0.3, 0.7, 1.1])
    rng = np.random.default_rng(3)
    p0 = rng.normal(size=4)
    p0 /= np.linalg.norm(p0)
    om0 = np.array([2.0, 0.1, 1.5])
    ukl = Uklad()
    ukl.dodajCzlon(Czlon(1, 2.0, Jb))
    ukl.grawitacja = False
    q0 = np.zeros(7)
    q0[3:7] = p0
    dp0 = 0.5*G(wektor_p(*p0)).T.dot(om0.reshape(3, 1)).ravel()
    getattr(ukl, met)(np.concatenate((q0, np.zeros(3), dp0)), 0.0, 1.0, dt)

    def kret(w):
        p = wektor_p(*w[3:7])
        dp = wektor_p(*w[10:14])
        return R(p).dot(Jb.dot(2*G(p).dot(dp))).ravel()
    L0 = kret(ukl.Y[0])
    return np.linalg.norm(kret(ukl.Y[-1]) - L0)/np.linalg.norm(L0)


def _wahadlo(met, dt, T=4.0):
    """Wahadlo na przegubie: (maks. blad energii, maks. naruszenie wiezow)."""
    ukl = Uklad()
    ukl.dodajCzlon(Czlon(1, 1.0, np.diag([1.0, 1.0, 1.0])))
    ukl.dodajWiez(Polaczenie_Obr(0, 1, wektor(0, 0, 0), wektor(0, 0, 2),
                                 wektor(1, 0, 0), wektor(0, 0, 1), wektor(0, 1, 0)))
    ukl.grawitacja = True
    q0 = np.zeros(7)
    q0[0] = -2*np.sin(0.5)
    q0[2] = -2*np.cos(0.5)
    q0[3:7] = u2p(np.array([0.0, 1.0, 0.0]), 0.5)
    getattr(ukl, met)(np.concatenate((q0, np.zeros(7))), 0.0, T, dt)
    E = [ukl.energia(ukl.Y[i]) for i in range(0, len(ukl.Y), max(1, len(ukl.Y)//40))]
    wz = max(np.linalg.norm(ukl.wiezyKP(w[:7])) for w in ukl.Y[::100])
    return max(abs(e - E[0]) for e in E), wz


def test_sym3_rzad_2():
    """Blad kretu bryly swobodnej maleje jak dt^2 (stosunek ~4 przy dt/2)."""
    b1 = _bryla('sym3', 8e-4)
    b2 = _bryla('sym3', 4e-4)
    assert b1/b2 == pytest.approx(4.0, rel=0.2)


def test_sym3_dokladniejszy_od_sym2():
    """Przy tym samym dt sym3 ma blad kretu i energii >=50x mniejszy."""
    assert _bryla('sym2', 4e-4)/_bryla('sym3', 4e-4) > 50
    e2, _ = _wahadlo('sym2', 2e-3)
    e3, _ = _wahadlo('sym3', 2e-3)
    assert e2/e3 > 50


def test_sym3_wiezy_maszynowe():
    """Wiezy (SHAKE/RATTLE wewnatrz kroku) spelnione na poziomie maszynowym."""
    _, wz = _wahadlo('sym3', 4e-3)
    assert wz < 1e-10


def test_sym3_duzy_krok_nadal_dokladny():
    """sym3 przy dt 10x wiekszym jest wciaz dokladniejszy niz sym2."""
    e2, _ = _wahadlo('sym2', 2e-3)
    e3, _ = _wahadlo('sym3', 2e-2)
    assert e3 < e2


def _damper_light(polniejawne, dt, c=200.0):
    """Lekki czlon (m=0.2) z mocnym tlumikiem liniowym do podstawy (c/m=1000,
    prog jawny dt<2m/c=2e-3). Zwraca norme predkosci koncowej (0 = wygaszone)
    albo None gdy symulacja wybucha."""
    ukl = Uklad()
    ukl.dodajCzlon(Czlon(1, 0.2, np.diag([0.02, 0.02, 0.02])))
    ukl.dodajSileWewn(SilaWewnProst(0, 1, wektor(0, 0, 0), wektor(0, 0, 0),
                                    5.0, 0.0, c, 0))
    ukl.grawitacja = False
    q0 = np.zeros(7)
    q0[3] = 1.0
    dq0 = np.zeros(7)
    dq0[0] = 3.0
    try:
        ukl.sym3(np.concatenate((q0, dq0)), 0.0, 1.0, dt, polniejawne=polniejawne)
        return np.linalg.norm(ukl.Y[-1][7:10])
    except RuntimeError:
        return None


def test_sym3_polniejawny_stabilny_przy_sztywnym_tlumieniu():
    """Przy dt > 2m/c jawny sym3 wybucha (albo nie wygasza), a polniejawny
    poprawnie tlumi predkosc do zera (stabilnosc bezwarunkowa w tlumieniu)."""
    # dt=4e-3 (c/m*dt/2=2): jawny wybucha, polniejawny wygasza
    assert _damper_light(False, 4e-3) is None
    vp = _damper_light(True, 4e-3)
    assert vp is not None and vp < 1e-2


def test_sym3_polniejawny_rzad_2():
    """Tryb polniejawny zachowuje rzad 2 na ukladzie niesztywnym (bryla)."""
    b1 = _bryla('sym3', 8e-4)      # jawny (bryla bez tlumienia -> zgodny)
    ukl = Uklad()
    # policz blad kretu dla polniejawnego przy dwoch dt
    def blad(dt):
        Jb = np.diag([0.3, 0.7, 1.1])
        rng = np.random.default_rng(3)
        p0 = rng.normal(size=4)
        p0 /= np.linalg.norm(p0)
        om0 = np.array([2.0, 0.1, 1.5])
        u = Uklad()
        u.dodajCzlon(Czlon(1, 2.0, Jb))
        u.grawitacja = False
        q0 = np.zeros(7)
        q0[3:7] = p0
        dp0 = 0.5*G(wektor_p(*p0)).T.dot(om0.reshape(3, 1)).ravel()
        u.sym3(np.concatenate((q0, np.zeros(3), dp0)), 0.0, 1.0, dt, polniejawne=True)

        def kret(w):
            p = wektor_p(*w[3:7])
            dp = wektor_p(*w[10:14])
            return R(p).dot(Jb.dot(2*G(p).dot(dp))).ravel()
        L0 = kret(u.Y[0])
        return np.linalg.norm(kret(u.Y[-1]) - L0)/np.linalg.norm(L0)
    assert blad(8e-4)/blad(4e-4) == pytest.approx(4.0, rel=0.2)


def test_sym3_pid_dziala():
    """Hook czlonu calkujacego (PID) dziala takze w sym3."""
    from uw_dyn import MomentWzgledny
    ukl = Uklad()
    ukl.dodajCzlon(Czlon(1, 1.0, np.diag([1.0, 1.0, 1.0])))
    ukl.dodajWiez(Polaczenie_Obr(0, 1, wektor(0, 0, 0), wektor(0, 0, 2),
                                 wektor(1, 0, 0), wektor(0, 0, 1), wektor(0, 1, 0)))
    ukl.dodajSileWewn(MomentWzgledny(0, 1, wektor(0, 1, 0), wektor(1, 0, 0),
                                     40.0, 0.8, 8.0, ki=60.0, calka_max=200.0))
    ukl.grawitacja = True
    q0 = np.zeros(7)
    q0[2] = -2
    q0[3] = 1
    ukl.sym3(np.concatenate((q0, np.zeros(7))), 0.0, 12.0, 2e-3)
    p = ukl.Y[-1][3:7]
    kat = 2*np.arctan2(p[2], p[0])
    assert abs(0.8 - kat) < 0.06     # PID znosi sag jak w sym2
