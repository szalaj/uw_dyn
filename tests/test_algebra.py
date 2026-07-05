# -*- coding: utf-8 -*-
# testy funkcji algebraicznych: wektory, kwaterniony, macierze obrotu

import numpy as np
import pytest

from uw_dyn import wektor, wektor_p, skew, R, G, u2p, EA_to_EP


def test_wektor_ksztalt():
    a = wektor(1, 2, 3)
    assert a.shape == (3, 1)
    assert np.allclose(a.ravel(), [1, 2, 3])


def test_wektor_p_ksztalt():
    p = wektor_p(1, 0, 0, 0)
    assert p.shape == (4, 1)


def test_skew_antysymetria():
    a = wektor(1.5, -2.0, 3.0)
    A = skew(a)
    assert A.shape == (3, 3)
    assert np.allclose(A, -A.T)


def test_skew_iloczyn_wektorowy():
    a = wektor(1.0, 2.0, 3.0)
    b = wektor(-4.0, 5.0, 0.5)
    assert np.allclose(skew(a).dot(b).ravel(),
                       np.cross(a.ravel(), b.ravel()))


def test_R_kwaternion_jednostkowy():
    p = wektor_p(1, 0, 0, 0)
    assert np.allclose(R(p), np.eye(3))


def test_R_ortonormalna():
    rng = np.random.default_rng(42)
    for _ in range(10):
        p = rng.normal(size=4)
        p = p / np.linalg.norm(p)
        Rot = R(wektor_p(*p))
        assert np.allclose(Rot.T.dot(Rot), np.eye(3), atol=1e-12)
        assert np.isclose(np.linalg.det(Rot), 1.0)


def test_u2p_norma_i_obrot():
    kat = 0.7
    p = u2p(np.array([0.0, 0.0, 1.0]), kat)
    assert np.isclose(np.linalg.norm(p), 1.0)
    # obrot wokol osi z przeksztalca os x na (cos, sin, 0)
    Rot = R(wektor_p(*p))
    assert np.allclose(Rot.dot(wektor(1, 0, 0)).ravel(),
                       [np.cos(kat), np.sin(kat), 0.0])


def test_G_wlasnosci():
    rng = np.random.default_rng(7)
    p = rng.normal(size=4)
    p = p / np.linalg.norm(p)
    pp = wektor_p(*p)
    Gp = G(pp)
    # G(p) p = 0
    assert np.allclose(Gp.dot(pp), np.zeros((3, 1)), atol=1e-12)
    # G G^T = I (dla kwaternionu jednostkowego)
    assert np.allclose(Gp.dot(Gp.T), np.eye(3), atol=1e-12)


def test_EA_to_EP_obrot_wokol_z():
    kat = 0.9
    e0, e1, e2, e3 = EA_to_EP(kat, 0.0, 0.0)
    oczekiwane = u2p(np.array([0.0, 0.0, 1.0]), kat)
    assert np.allclose([e0, e1, e2, e3], oczekiwane)
