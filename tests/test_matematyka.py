# -*- coding: utf-8 -*-
# Audyt matematyczny biblioteki (weryfikacja numeryczna sformulowania):
#  1. jakobiany wiezow = pochodne kierunkowe Phi wzdluz kierunkow STYCZNYCH
#     (dp _|_ p; wzdluz kierunku radialnego jakobiany sa z definicji niedokladne,
#     patrz CLAUDE.md - dlatego iteracje Newtona normalizuja kwaterniony),
#  2. gammaK = -(dJ/dt) v (tozsamosc Phi'' = J a + Jdot v => J a = gamma),
#  3. sily potencjalne: Q = -dV/dq wzdluz kierunkow stycznych (zawiasy: na
#     rozmaitosci przegubu),
#  4. bryla swobodna: zachowanie kretu w ukladzie swiata; blad maleje liniowo
#     z dt (rzad 1 poljawnego Eulera) => sformulowanie rownan egzaktne.

import numpy as np
import pytest

from uw_dyn import (Para_Sferyczna, Polaczenie_Obr, Polaczenie_Cyl,
                    Polaczenie_Przes, Para_Prostopadla, Para_Prostopadla_D,
                    SilaWewnProst, MomentWzgledny, MomentSferyczny,
                    OgranicznikKata, OgranicznikStozka, SilaKontaktu,
                    Uklad, Czlon, wektor, u2p, R, G, wektor_p, mnoz_kwaterniony)

N = 2
EPS = 1e-6
rng = np.random.default_rng(7)


def _kw():
    p = rng.normal(size=4)
    return p/np.linalg.norm(p)


def _q():
    q = np.zeros(7*N)
    q[0:3*N] = rng.normal(size=3*N)
    for k in range(N):
        q[3*N+4*k:3*N+4*k+4] = _kw()
    return q


def _styczny(q):
    d = np.zeros(7*N)
    d[0:3*N] = rng.normal(size=3*N)
    for k in range(N):
        p = q[3*N+4*k:3*N+4*k+4]
        dp = rng.normal(size=4)
        dp -= p*np.dot(p, dp)
        d[3*N+4*k:3*N+4*k+4] = dp
    return d/np.linalg.norm(d)


def _q_do_jak(x):
    out = np.zeros_like(x)
    for k in range(N):
        out[7*k:7*k+3] = x[3*k:3*k+3]
        out[7*k+3:7*k+7] = x[3*N+4*k:3*N+4*k+4]
    return out


def _pelny_jak(w, q):
    Fqi, Fqj = w.jakobianK(q, N)
    J = np.zeros((Fqi.shape[0], 7*N))
    if w.i != 0:
        J[:, 7*(w.i-1):7*(w.i-1)+7] = Fqi
    J[:, 7*(w.j-1):7*(w.j-1)+7] = Fqj
    return J


def _wiezy_testowe():
    v = lambda: wektor(*rng.normal(size=3))
    return [
        Para_Sferyczna(0, 1, v(), v()),
        Para_Sferyczna(1, 2, v(), v()),
        Para_Prostopadla(0, 1, v(), v()),
        Para_Prostopadla(1, 2, v(), v()),
        Para_Prostopadla_D(0, 1, v(), v(), v()),
        Para_Prostopadla_D(1, 2, v(), v(), v()),
        Polaczenie_Obr(1, 2, wektor(0, 0, 1), wektor(0, 0, -1),
                       wektor(1, 0, 0), wektor(0, 0, 1), wektor(0, 1, 0)),
        Polaczenie_Cyl(1, 2, wektor(0, 0, 1), wektor(0, 0, -1),
                       wektor(1, 0, 0), wektor(0, 1, 0), wektor(0, 0, 1)),
        Polaczenie_Przes(1, 2, wektor(0, 0, 1), wektor(0, 0, -1),
                         wektor(1, 0, 0), wektor(0, 1, 0), wektor(0, 1, 0),
                         wektor(0, 0, 1)),
    ]


def test_jakobiany_styczne():
    """J * d == pochodna kierunkowa Phi wzdluz kazdego kierunku stycznego."""
    for w in _wiezy_testowe():
        for _ in range(3):
            q = _q()
            Ja = _pelny_jak(w, q)
            for _ in range(5):
                d = _styczny(q)
                num = (np.asarray(w.wiezyK(q+EPS*d, N))
                       - np.asarray(w.wiezyK(q-EPS*d, N)))/(2*EPS)
                an = Ja @ _q_do_jak(d)
                assert np.abs(num.ravel() - an.ravel()).max() < 1e-6, type(w).__name__


def test_gamma_rowna_minus_Jdot_v():
    """gammaK == -(dJ/dt) v dla losowych stanow (dp styczne do sfery)."""
    for w in _wiezy_testowe():
        for _ in range(3):
            q = _q()
            dq = _styczny(q)*rng.uniform(0.5, 2.0)
            v = _q_do_jak(dq)
            Jdot = (_pelny_jak(w, q+EPS*dq) - _pelny_jak(w, q-EPS*dq))/(2*EPS)
            gamma_num = -(Jdot @ v)
            gamma_an = np.atleast_2d(w.gammaK(q, dq, N)).ravel()
            assert np.abs(gamma_an - gamma_num.ravel()).max() < 1e-6, type(w).__name__


def _sprawdz_potencjal(s, q, kierunki):
    Qr_i, Qp_i, Qr_j, Qp_j = s.sila(q, np.zeros(7*N), N)
    Q = np.zeros(7*N)
    if s.i != 0:
        Q[3*(s.i-1):3*(s.i-1)+3] = np.ravel(Qr_i)
        Q[3*N+4*(s.i-1):3*N+4*(s.i-1)+4] = np.ravel(Qp_i)
    Q[3*(s.j-1):3*(s.j-1)+3] += np.ravel(Qr_j)
    Q[3*N+4*(s.j-1):3*N+4*(s.j-1)+4] += np.ravel(Qp_j)
    for d in kierunki:
        dV = (s.energia_potencjalna(q+EPS*d, N)
              - s.energia_potencjalna(q-EPS*d, N))/(2*EPS)
        assert abs(dV + Q @ d) < 1e-5, type(s).__name__


def test_sily_potencjalne_ogolne():
    """Q == -dV/dq (styczne) dla sil o dowolnej konfiguracji."""
    sily = [
        SilaWewnProst(1, 2, wektor(0.1, 0, 0.3), wektor(-0.2, 0.1, 0), 37.0, 0.4, 0.0, 0.0),
        SilaWewnProst(0, 1, wektor(0.5, 0.2, 0), wektor(-0.2, 0.1, 0), 37.0, 0.4, 0.0, 0.0),
        MomentSferyczny(1, 2, 25.0, 0.0, p_cel=u2p(np.array([0.0, 0.6, 0.8]), 0.6)),
        OgranicznikStozka(1, 2, wektor(0, 0, 1), 0.1, 300.0, 0.0),
        SilaKontaktu(1, wektor(0.05, 0.02, -0.3), k=1000.0, c=0.0, mu=0.0),
    ]
    for s in sily:
        for _ in range(3):
            q = _q()
            _sprawdz_potencjal(s, q, [_styczny(q) for _ in range(5)])


def test_sily_zawiasowe_na_rozmaitosci():
    """Zawiasy (MomentWzgledny, OgranicznikKata): Q == -dV/d(kat) na
    rozmaitosci przegubu (obrot wzgledny czysto wokol osi zawiasu)."""
    OS_L = np.array([0.0, 1.0, 0.0])
    for cls, kw in ((MomentWzgledny, dict(k=25.0, theta_cel=0.4, c=0.0)),
                    (OgranicznikKata, dict(kat_min=-0.2, kat_max=0.1, k=500.0, c=0.0))):
        s = cls(1, 2, wektor(0, 1, 0), wektor(0, 0, 1), **kw)
        for _ in range(6):
            q = _q()
            p_i = q[3*N:3*N+4]
            p_j = mnoz_kwaterniony(p_i, u2p(OS_L, rng.uniform(-2.5, 2.5)))
            q[3*N+4:3*N+8] = p_j
            a_g = R(wektor_p(*p_i)).dot(OS_L.reshape(3, 1))
            dpj = (0.5*G(wektor_p(*p_j)).T.dot(a_g)).ravel()
            d = np.zeros(7*N)
            d[3*N+4:3*N+8] = dpj
            _sprawdz_potencjal(s, q, [d/np.linalg.norm(d)])


def test_bryla_swobodna_zbieznosc_kretu():
    """Bryla swobodna (asymetryczny tensor, os posrednia): kret w ukladzie
    swiata zachowany, a blad maleje liniowo z dt (rzad 1 integratora) =>
    sformulowanie rownan obrotu (4G^T J G, 8Gdot^T J Gdot p) jest egzaktne."""
    Jb = np.diag([0.3, 0.7, 1.1])
    p0 = _kw()
    om0 = np.array([2.0, 0.1, 1.5])

    def blad(dt):
        ukl = Uklad()
        ukl.dodajCzlon(Czlon(1, 2.0, Jb))
        ukl.grawitacja = False
        q0 = np.zeros(7)
        q0[3:7] = p0
        dp0 = 0.5*G(wektor_p(*p0)).T.dot(om0.reshape(3, 1)).ravel()
        ukl.sym2(np.concatenate((q0, np.zeros(3), dp0)), 0.0, 1.0, dt)

        def kret(w):
            p = wektor_p(*w[3:7])
            dp = wektor_p(*w[10:14])
            return R(p).dot(Jb.dot(2*G(p).dot(dp))).ravel()
        L0 = kret(ukl.Y[0])
        return np.linalg.norm(kret(ukl.Y[-1]) - L0)/np.linalg.norm(L0)

    b1, b2 = blad(4e-4), blad(2e-4)
    assert b1 < 1e-3                      # kret zachowany (maly blad)
    assert b1/b2 == pytest.approx(2.0, rel=0.15)   # zbieznosc rzedu 1
