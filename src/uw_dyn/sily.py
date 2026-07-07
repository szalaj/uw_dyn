# -*- coding: utf-8 -*-
# uw_dyn: dynamika 3D ukladow wieloczlonowych
# autor: Marcin Szalajski (2016; pakiet od 2026)

"""Sily: elementy sprezysto-tlumiace i sily zewnetrzne."""

import numpy as np

from uw_dyn.algebra import (r_i, dr_i, p_i, dp_i, R, G, skew,
                            mnoz_kwaterniony, sprzezenie_kwaternionu,
                            kwaternion_na_wektor_obrotu)


class SilaWewnProst:
    """Element sprezysto-tlumiacy z sila stala (sprezyna k o dlugosci
    swobodnej l0, tlumik c, sila F) miedzy punktem A czlonu i oraz punktem B
    czlonu j. Z flaga tylko_rozciaganie=True dziala jak lina: sila znika,
    gdy element jest krotszy niz dlugosc swobodna."""
    def __init__(self, i, j, sA_i, sB_j, k,l0, c, F, tylko_rozciaganie=False):
            self.i = i
            self.j = j
            self.sA_i = sA_i
            self.sB_j = sB_j
            self.k = k
            self.l0 = l0
            self.c = c
            self.F = F
            self.tylko_rozciaganie = tylko_rozciaganie

    def energia_potencjalna(self,q,N):
        """Energia sprezysta elementu (czlon tlumika i sily stalej pominiety)."""
        if self.k == 0:
            return 0.0
        l = self.dlugosc(q,N)
        if self.tylko_rozciaganie and l < self.l0:
            return 0.0
        return 0.5*self.k*(l - self.l0)**2
            
    def dlugosc(self,q,N):
        """Aktualna dlugosc elementu (odleglosc punktow zaczepienia)."""
        rj = r_i(self.j,q)
        pj = p_i(self.j,q,N)
        Rj = R(pj)
        if self.i == 0:
            dij = rj + Rj.dot(self.sB_j) - self.sA_i
        else:
            ri = r_i(self.i,q)
            pi = p_i(self.i,q,N)
            Ri = R(pi)
            dij = rj + Rj.dot(self.sB_j) - ri - Ri.dot(self.sA_i)
        return float(np.linalg.norm(dij))

    def sila(self,q,dq,N):
        rj = r_i(self.j,q)
        drj = dr_i(self.j,dq)
        
        pj = p_i(self.j,q,N)
        dpj = dp_i(self.j,dq,N)
        Rj = R(pj)
        Gj = G(pj)

        om_j = 2*Gj.dot(dpj)

        Qr_i = np.zeros((3,1))
        Qp_i = np.zeros((4,1))         
        Qr_j = np.zeros((3,1))
        Qp_j = np.zeros((4,1))  
              
        if self.i==0:

            dij = rj + Rj.dot(self.sB_j) - self.sA_i
            l = np.sqrt( dij.transpose().dot(dij) )

            if self.tylko_rozciaganie and l < self.l0:
                return Qr_i, Qp_i, Qr_j, Qp_j

            if l>0.01:
                pom = drj - Rj.dot(skew(self.sB_j)).dot(om_j)

                dl = (dij/l).transpose().dot(pom)
                f = self.k*(l-self.l0) + self.c*dl + self.F

                Qr_j = -(f/l)*dij
                Qp_j = -(f/l)*2*Gj.transpose().dot(skew(self.sB_j)).dot(Rj.transpose()).dot(dij)
             
        else:
            ri = r_i(self.i,q)
            dri = dr_i(self.i,dq)
            
            pi = p_i(self.i,q,N)
            dpi = dp_i(self.i,dq,N)
            Ri = R(pi)
            Gi = G(pi)


            om_i = 2*Gi.dot(dpi)
            
            dij = rj + Rj.dot(self.sB_j) - ri - Ri.dot(self.sA_i)
            l = np.sqrt( dij.transpose().dot(dij) )

            if self.tylko_rozciaganie and l < self.l0:
                return Qr_i, Qp_i, Qr_j, Qp_j

            if l>0.01:
                pom = drj - Rj.dot(skew(self.sB_j)).dot(om_j) - dri + Ri.dot(skew(self.sA_i)).dot(om_i)

                dl = (dij/l).transpose().dot(pom)
                f = self.k*(l-self.l0) + self.c*dl + self.F
                #print('i: ',self.i,' j: ', self.j, ' f: ',f)

                Qr_i = (f/l)*dij
                Qp_i = (f/l)*2*Gi.transpose().dot(skew(self.sA_i)).dot(Ri.transpose()).dot(dij)
                
                Qr_j = -(f/l)*dij
                Qp_j = -(f/l)*2*Gj.transpose().dot(skew(self.sB_j)).dot(Rj.transpose()).dot(dij)
            
        return Qr_i, Qp_i, Qr_j, Qp_j
        
        
class SilaWPunkcie:
    """Sila zaczepiona w punkcie ciala, zadana w ukladzie ciala
    (follower force), np. ciag wirnika drona.

    czlon: numer czlonu; s_punkt: punkt zaczepienia w ukladzie ciala;
    f_lokalna: wektor sily w ukladzie ciala (podazy za orientacja).
    Wektor sily mozna podmieniac miedzy segmentami symulacji
    (dyskretny regulator)."""

    def __init__(self, czlon, s_punkt, f_lokalna):
        self.i = 0
        self.j = czlon
        self.s_punkt = s_punkt
        self.f_lokalna = f_lokalna

    def energia_potencjalna(self,q,N):
        return 0.0

    def sila(self,q,dq,N):
        pj = p_i(self.j,q,N)
        Rj = R(pj)
        Gj = G(pj)

        Qr_i = np.zeros((3,1))
        Qp_i = np.zeros((4,1))
        # sila globalna: R f'; moment w ukladzie ciala: s' x f'
        Qr_j = Rj.dot(self.f_lokalna)
        n_lok = skew(self.s_punkt).dot(self.f_lokalna)
        Qp_j = 2*Gj.transpose().dot(n_lok)
        return Qr_i, Qp_i, Qr_j, Qp_j


class SilaKontaktu:
    """Jednostronny kontakt punktu ciala z podlozem z = 0 (model penalty).

    Gdy punkt zaczepienia jest pod powierzchnia: sila normalna
    N = max(0, k*wnikanie - c*predkosc_pionowa) oraz tarcie styczne
    (regularyzowany Coulomb): T = -mu*N * v_t/(|v_t| + eps).

    czlon: numer czlonu; s_punkt: punkt stopy w ukladzie ciala;
    k, c: sztywnosc i tlumienie podloza; mu: wspolczynnik tarcia."""

    def __init__(self, czlon, s_punkt, k=2.0e4, c=200.0, mu=0.8, eps=0.01):
        self.i = 0
        self.j = czlon
        self.s_punkt = s_punkt
        self.k = k
        self.c = c
        self.mu = mu
        self.eps = eps

    def _punkt(self,q,N):
        rj = r_i(self.j,q)
        pj = p_i(self.j,q,N)
        Rj = R(pj)
        return rj + Rj.dot(self.s_punkt), Rj

    def energia_potencjalna(self,q,N):
        pw, _ = self._punkt(q,N)
        wnikanie = -float(pw[2,0])
        if wnikanie <= 0:
            return 0.0
        return 0.5*self.k*wnikanie**2

    def sila(self,q,dq,N):
        Qr_i = np.zeros((3,1))
        Qp_i = np.zeros((4,1))
        Qr_j = np.zeros((3,1))
        Qp_j = np.zeros((4,1))

        pw, Rj = self._punkt(q,N)
        wnikanie = -float(pw[2,0])
        if wnikanie <= 0:
            return Qr_i, Qp_i, Qr_j, Qp_j

        # predkosc punktu stopy: dr + R (om' x s')
        drj = dr_i(self.j,dq)
        pj = p_i(self.j,q,N)
        dpj = dp_i(self.j,dq,N)
        Gj = G(pj)
        om_j = 2*Gj.dot(dpj)
        v = drj + Rj.dot(skew(om_j)).dot(self.s_punkt)

        Fn = max(0.0, self.k*wnikanie - self.c*float(v[2,0]))
        vt = np.array([float(v[0,0]), float(v[1,0])])
        Ft = -self.mu*Fn*vt/(np.linalg.norm(vt) + self.eps)

        F_glob = np.array([[Ft[0]], [Ft[1]], [Fn]])
        Qr_j = F_glob
        # moment w ukladzie ciala: s' x (R^T F)
        n_lok = skew(self.s_punkt).dot(Rj.transpose().dot(F_glob))
        Qp_j = 2*Gj.transpose().dot(n_lok)
        return Qr_i, Qp_i, Qr_j, Qp_j


def _najblizsze_odcinki(Ci, ui, hi, Cj, uj, hj):
    """Najblizsze punkty dwoch odcinkow (kapsuly): odcinek i o srodku Ci,
    kierunku ui (jedn.), polowie dlugosci hi; analogicznie j. Zwraca (Pi, Pj).
    Standardowy algorytm z klamrowaniem parametrow do [-h, h]."""
    r = Ci - Cj
    b = ui.T.dot(uj).item()
    d = ui.T.dot(r).item()
    e = uj.T.dot(r).item()
    denom = 1.0 - b*b                      # = |ui|^2|uj|^2 - (ui.uj)^2
    if denom > 1e-9:
        sa = max(-hi, min(hi, (b*e - d)/denom))
    else:                                  # odcinki prawie rownolegle
        sa = 0.0
    sb = b*sa + e
    if sb < -hj:
        sb = -hj
        sa = max(-hi, min(hi, b*sb - d))
    elif sb > hj:
        sb = hj
        sa = max(-hi, min(hi, b*sb - d))
    return Ci + sa*ui, Cj + sb*uj


class SilaUderzenia:
    """Kontakt bryla-bryla (penalty) miedzy dwiema KAPSULAMI (odcinek wzdluz
    lokalnej osi z + promien). Przypadek punkt-kapsula to polowa_wys_i = 0.

    Gdy najblizsze punkty kapsul zblizaja sie ponizej sumy promieni: sila
    normalna N = max(0, k*wnikanie - c*predkosc_zblizania) wzdluz normalnej
    laczacej te punkty. Rowna i przeciwna na oba ciala, zaczepiona w punktach
    kontaktu, wiec przekazuje ped i generuje obrot/wahniecie. Uzycia:
    piesc/stopa (kapsula i lub punkt) vs worek; tulow vs worek (by worek nie
    przenikal ciala); konczyna vs przeciwnik (dwaj zawodnicy).

    czlon_i, s_i: cialo uderzajace i srodek jego odcinka w ukladzie ciala
    (dla punktu: sam punkt, polowa_wys_i=0); czlon_j: cialo trafiane (odcinek
    wzdluz lokalnej z, srodek w COM); promien = suma promieni kapsul;
    polowa_wys_j, polowa_wys_i: polowy dlugosci rdzeni; k, c: sztywnosc i
    tlumienie kontaktu."""

    def __init__(self, czlon_i, s_i, czlon_j, promien, polowa_wys_j,
                 polowa_wys_i=0.0, k=3.0e4, c=150.0):
        self.i = czlon_i
        self.j = czlon_j
        self.s_i = np.asarray(s_i, float).reshape(3, 1)
        self.promien = promien
        self.polowa_wys_j = polowa_wys_j
        self.polowa_wys_i = polowa_wys_i
        self.k = k
        self.c = c
        self.zeruj_metryki()

    def zeruj_metryki(self):
        """Zeruje liczniki metryk uderzenia (przed symulacja)."""
        self.F_szczyt = 0.0        # szczytowa sila normalna kontaktu [N]
        self.suma_Fn = 0.0         # suma sil po wywolaniach (impuls = suma*dt)
        self.n_kontaktu = 0        # liczba wywolan z kontaktem

    def impuls(self, dt):
        """Impuls uderzenia [N s]: calka sily po czasie. Wazne dla sym2
        (jedno wywolanie sila() na krok dt); dla sym3 przeskaluj liczba
        wywolan na krok."""
        return self.suma_Fn*dt

    def _geometria(self, q, N):
        """Najblizsze punkty kapsul i wnikanie."""
        ri, rj = r_i(self.i, q), r_i(self.j, q)
        pi_, pj = p_i(self.i, q, N), p_i(self.j, q, N)
        Ri, Rj = R(pi_), R(pj)
        Ci = ri + Ri.dot(self.s_i)             # srodek odcinka i
        ui, uj = Ri[:, 2:3], Rj[:, 2:3]        # lokalne osie z (swiat)
        Pi, Pj = _najblizsze_odcinki(Ci, ui, self.polowa_wys_i,
                                     rj, uj, self.polowa_wys_j)
        d = Pi - Pj
        dist = float(np.linalg.norm(d))
        wnikanie = self.promien - dist
        return Pi, Pj, d, dist, wnikanie, Ri, Rj, pi_, pj

    def energia_potencjalna(self, q, N):
        _, _, _, _, wnikanie, *_ = self._geometria(q, N)
        if wnikanie <= 0:
            return 0.0
        return 0.5*self.k*wnikanie**2

    def sila(self, q, dq, N):
        Qr_i = np.zeros((3, 1))
        Qp_i = np.zeros((4, 1))
        Qr_j = np.zeros((3, 1))
        Qp_j = np.zeros((4, 1))

        Pi, Pj, d, dist, wnikanie, Ri, Rj, pi_, pj = self._geometria(q, N)
        if wnikanie <= 0 or dist < 1e-9:
            return Qr_i, Qp_i, Qr_j, Qp_j
        n = d/dist                              # od punktu j do punktu i (3,1)

        # ramiona kontaktu w ukladach cial i predkosci punktow kontaktu
        s_i = Ri.transpose().dot(Pi - r_i(self.i, q))
        s_j = Rj.transpose().dot(Pj - r_i(self.j, q))
        Gi, Gj = G(pi_), G(pj)
        om_i = 2*Gi.dot(dp_i(self.i, dq, N))
        om_j = 2*Gj.dot(dp_i(self.j, dq, N))
        v_i = dr_i(self.i, dq) + Ri.dot(skew(om_i)).dot(s_i)
        v_j = dr_i(self.j, dq) + Rj.dot(skew(om_j)).dot(s_j)
        v_zbl = (v_i - v_j).T.dot(n).item()     # predkosc zblizania wzdluz n

        Fn = max(0.0, self.k*wnikanie - self.c*v_zbl)
        F = Fn*n                                # sila na cialo i (odpychajaca)

        # metryki uderzenia (szczyt i suma do impulsu)
        if Fn > 0.0:
            self.F_szczyt = max(self.F_szczyt, Fn)
            self.suma_Fn += Fn
            self.n_kontaktu += 1

        Qr_i = F
        Qp_i = 2*Gi.transpose().dot(skew(s_i).dot(Ri.transpose().dot(F)))
        Qr_j = -F                               # reakcja na cialo j
        Qp_j = 2*Gj.transpose().dot(skew(s_j).dot(Rj.transpose().dot(-F)))
        return Qr_i, Qp_i, Qr_j, Qp_j


class MomentWzgledny:
    """Aktuator obrotowy (sprezyna-tlumik) w przegubie miedzy czlonami i, j
    wokol osi obrotu. Moment tau = k*(theta_cel - theta) - c*dtheta jest
    przykladany +na czlon j oraz -na czlon i (para wewnetrzna).

    axis, ref: os obrotu i wektor odniesienia (prostopadly do osi), oba
    w ukladzie ciala; przy zbudowaniu ukladu czlony i, j maja byc ustawione
    tak, ze przy theta=0 wektory ref obu czlonow sie pokrywaja. theta_cel
    mozna zmieniac miedzy segmentami symulacji (sterowanie).

    moment_max: opcjonalne ograniczenie wartosci bezwzglednej momentu [N*m]
    (saturacja jak w realnym napedzie/miesniu); None = bez limitu. UWAGA: przy
    saturacji czlon k*blad przewyzsza limit, wiec tlumienie (rowniez objete
    limitem) traci wplyw i przegub moze oscylowac; realny staw ma pasywne
    tlumienie tkanek, ktore mozna dodac osobnym MomentWzgledny (k=0, c>0,
    bez moment_max).

    ki: wzmocnienie czlonu calkujacego (regulator PID zamiast PD). Znosi blad
    ustalony pod stalym obciazeniem (grawitacja) - staw dochodzi do celu bez
    sagu. Calke aktualizuje sym2 raz na krok (metoda aktualizuj_calke); dziala
    tylko w sym2, nie w sym (RK45). calka_max: anti-windup - ogranicza wartosc
    czlonu calkujacego |ki*calka| [N*m].

    Dla i=0 (podstawa) os i ref sa w ukladzie globalnym."""

    def __init__(self, i, j, axis, ref, k, theta_cel, c, moment_max=None,
                 ki=0.0, calka_max=None):
        self.i = i
        self.j = j
        self.a = np.asarray(axis, dtype=float).reshape(3, 1)
        self.a /= np.linalg.norm(self.a)
        self.u = np.asarray(ref, dtype=float).reshape(3, 1)
        self.u /= np.linalg.norm(self.u)
        self.k = k
        self.theta_cel = theta_cel
        self.c = c
        self.moment_max = moment_max
        self.ki = ki
        self.calka_max = calka_max
        self._calka = 0.0

    def _blad_kata(self, theta):
        # blad zawiniety do (-pi, pi] (najkrotsza droga)
        return (self.theta_cel - theta + np.pi) % (2*np.pi) - np.pi

    def aktualizuj_calke(self, q, N, dt):
        """Krok czlonu calkujacego (wywolywane raz na krok przez sym2)."""
        if self.ki == 0:
            return
        self._calka += self._blad_kata(self.kat(q, N))*dt
        if self.calka_max is not None:
            lim = self.calka_max/self.ki
            self._calka = max(-lim, min(lim, self._calka))

    def energia_potencjalna(self, q, N):
        if self.k == 0:
            return 0.0
        return 0.5*self.k*(self.kat(q, N) - self.theta_cel)**2

    def _ramki(self, q, dq, N):
        pj = p_i(self.j, q, N)
        Rj = R(pj)
        Gj = G(pj)
        om_j = Rj.dot(2*Gj.dot(dp_i(self.j, dq, N)))  # predkosc katowa globalna
        u_j = Rj.dot(self.u)
        if self.i == 0:
            a_g = self.a
            u_i = self.u
            om_i = np.zeros((3, 1))
            Ri = Gi = None
        else:
            pi = p_i(self.i, q, N)
            Ri = R(pi)
            Gi = G(pi)
            a_g = Ri.dot(self.a)
            u_i = Ri.dot(self.u)
            om_i = Ri.dot(2*Gi.dot(dp_i(self.i, dq, N)))
        return Rj, Gj, om_j, u_j, a_g, u_i, om_i, Ri, Gi

    def kat(self, q, N):
        """Kat wzgledny theta czlonu j wobec i wokol osi obrotu."""
        pj = p_i(self.j, q, N)
        u_j = R(pj).dot(self.u)
        if self.i == 0:
            a_g, u_i = self.a, self.u
        else:
            Ri = R(p_i(self.i, q, N))
            a_g, u_i = Ri.dot(self.a), Ri.dot(self.u)
        sin_cz = float(np.cross(u_i.ravel(), u_j.ravel()).dot(a_g.ravel()))
        cos_cz = float(u_i.ravel().dot(u_j.ravel()))
        return np.arctan2(sin_cz, cos_cz)

    def sila(self, q, dq, N):
        Rj, Gj, om_j, u_j, a_g, u_i, om_i, Ri, Gi = self._ramki(q, dq, N)

        sin_cz = float(np.cross(u_i.ravel(), u_j.ravel()).dot(a_g.ravel()))
        cos_cz = float(u_i.ravel().dot(u_j.ravel()))
        theta = np.arctan2(sin_cz, cos_cz)
        dtheta = float(a_g.ravel().dot((om_j - om_i).ravel()))

        # blad kata zawiniety do (-pi, pi] -> najkrotsza droga; bez tego kat
        # przechodzacy przez +-pi (np. biodro ~pi) prowadzilby staw "dookola"
        blad = self._blad_kata(theta)
        tau = self.k*blad - self.c*dtheta + self.ki*self._calka
        if self.moment_max is not None:
            tau = max(-self.moment_max, min(self.moment_max, tau))
        M = tau*a_g  # moment globalny na czlon j

        Qr_i = np.zeros((3, 1))
        Qr_j = np.zeros((3, 1))
        Qp_j = 2*Gj.transpose().dot(Rj.transpose().dot(M))
        if self.i == 0:
            Qp_i = np.zeros((4, 1))
        else:
            Qp_i = -2*Gi.transpose().dot(Ri.transpose().dot(M))
        return Qr_i, Qp_i, Qr_j, Qp_j


class MomentSferyczny:
    """Napedzany staw kulisty (3 stopnie swobody): aktuator momentowy
    sprowadzajacy orientacje czlonu j do zadanej orientacji wzgledem czlonu i
    (przestrzenna sprezyna-tlumik na obrocie 3D). Uzywac razem z
    Para_Sferyczna, ktora wiaze polozenie (kula w panewce); ta klasa steruje
    tylko orientacja.

    p_cel: docelowy obrot wzgledny (parametry Eulera, z ukladu i do j);
    domyslnie [1,0,0,0] = czlon j ustawiony jak i. p_cel mozna podmieniac
    miedzy segmentami symulacji (sterowanie 3D, np. bark, biodro). Moment
    liczony jest w ukladzie globalnym: M = -k*phi - c*(om_j - om_i), gdzie
    phi to wektor obrotu bledu (od orientacji docelowej do biezacej), i
    przykladany jako para wewnetrzna (+M na j, -M na i).

    moment_max: opcjonalne ograniczenie wartosci bezwzglednej (normy) momentu
    [N*m] (saturacja jak w realnym napedzie/miesniu); None = bez limitu.

    ki: wzmocnienie czlonu calkujacego (PID) - calkuje wektor bledu obrotu,
    znosi blad ustalony pod grawitacja. Calke aktualizuje sym2 raz na krok
    (aktualizuj_calke). calka_max: anti-windup na normie |ki*calka| [N*m]."""

    def __init__(self, i, j, k, c, p_cel=None, moment_max=None,
                 ki=0.0, calka_max=None):
        self.i = i
        self.j = j
        self.k = k
        self.c = c
        self.p_cel = (np.array([1.0, 0.0, 0.0, 0.0]) if p_cel is None
                      else np.asarray(p_cel, dtype=float).ravel())
        self.moment_max = moment_max
        self.ki = ki
        self.calka_max = calka_max
        self._calka = np.zeros((3, 1))

    def aktualizuj_calke(self, q, N, dt):
        """Krok czlonu calkujacego (wywolywane raz na krok przez sym2)."""
        if self.ki == 0:
            return
        self._calka = self._calka + self._blad(q, N).reshape(3, 1)*dt
        if self.calka_max is not None:
            n = self.ki*float(np.linalg.norm(self._calka))
            if n > self.calka_max:
                self._calka = self._calka*(self.calka_max/n)

    def _blad(self, q, N):
        """Wektor obrotu bledu (uklad globalny) orientacji czlonu j."""
        pj = p_i(self.j, q, N).ravel()
        pi = (np.array([1.0, 0.0, 0.0, 0.0]) if self.i == 0
              else p_i(self.i, q, N).ravel())
        p_docelowy = mnoz_kwaterniony(pi, self.p_cel)   # zadana orientacja j
        q_err = mnoz_kwaterniony(pj, sprzezenie_kwaternionu(p_docelowy))
        return kwaternion_na_wektor_obrotu(q_err)

    def kat(self, q, N):
        """Kat bledu orientacji (rad)."""
        return float(np.linalg.norm(self._blad(q, N)))

    def energia_potencjalna(self, q, N):
        if self.k == 0:
            return 0.0
        return 0.5*self.k*self.kat(q, N)**2

    def sila(self, q, dq, N):
        pj = p_i(self.j, q, N)
        Rj = R(pj)
        Gj = G(pj)
        om_j = Rj.dot(2*Gj.dot(dp_i(self.j, dq, N)))    # predkosc katowa globalna

        if self.i == 0:
            Ri = Gi = None
            om_i = np.zeros((3, 1))
        else:
            pi = p_i(self.i, q, N)
            Ri = R(pi)
            Gi = G(pi)
            om_i = Ri.dot(2*Gi.dot(dp_i(self.i, dq, N)))

        phi = self._blad(q, N).reshape(3, 1)            # blad orientacji (global)
        om_rel = om_j - om_i
        M = -self.k*phi - self.c*om_rel - self.ki*self._calka   # PID, moment na j
        if self.moment_max is not None:
            nrm = float(np.linalg.norm(M))
            if nrm > self.moment_max:
                M = M*(self.moment_max/nrm)

        Qr_i = np.zeros((3, 1))
        Qr_j = np.zeros((3, 1))
        Qp_j = 2*Gj.transpose().dot(Rj.transpose().dot(M))
        if self.i == 0:
            Qp_i = np.zeros((4, 1))
        else:
            Qp_i = 2*Gi.transpose().dot(Ri.transpose().dot(-M))
        return Qr_i, Qp_i, Qr_j, Qp_j


class OgranicznikKata:
    """Miekki ogranicznik zakresu ruchu przegubu obrotowego 1 DOF
    (np. lokiec, kolano). Jednostronna sprezyna-tlumik aktywna tylko, gdy kat
    wzgledny wyjdzie poza [kat_min, kat_max]: moment przywraca do zakresu, a
    tlumienie dziala tylko przy wchodzeniu glebiej (nie wysysa energii przy
    powrocie). Uzywac obok przegubu (Polaczenie_Obr) i jego napedu
    (MomentWzgledny); geometria (i, j, axis, ref) taka sama jak w MomentWzgledny.

    Dla i=0 (podstawa) os i ref sa w ukladzie globalnym."""

    def __init__(self, i, j, axis, ref, kat_min, kat_max, k, c):
        self.i = i
        self.j = j
        self.a = np.asarray(axis, dtype=float).reshape(3, 1)
        self.a /= np.linalg.norm(self.a)
        self.u = np.asarray(ref, dtype=float).reshape(3, 1)
        self.u /= np.linalg.norm(self.u)
        self.kat_min = kat_min
        self.kat_max = kat_max
        self.k = k
        self.c = c

    def kat(self, q, N):
        pj = p_i(self.j, q, N)
        u_j = R(pj).dot(self.u)
        if self.i == 0:
            a_g, u_i = self.a, self.u
        else:
            Ri = R(p_i(self.i, q, N))
            a_g, u_i = Ri.dot(self.a), Ri.dot(self.u)
        sin_cz = float(np.cross(u_i.ravel(), u_j.ravel()).dot(a_g.ravel()))
        cos_cz = float(u_i.ravel().dot(u_j.ravel()))
        return np.arctan2(sin_cz, cos_cz)

    def _przekroczenie(self, theta):
        if theta > self.kat_max:
            return theta - self.kat_max      # dodatnie: powyzej gornej granicy
        if theta < self.kat_min:
            return theta - self.kat_min      # ujemne: ponizej dolnej granicy
        return 0.0

    def energia_potencjalna(self, q, N):
        if self.k == 0:
            return 0.0
        p = self._przekroczenie(self.kat(q, N))
        return 0.5*self.k*p*p

    def sila(self, q, dq, N):
        pj = p_i(self.j, q, N)
        Rj = R(pj)
        Gj = G(pj)
        om_j = Rj.dot(2*Gj.dot(dp_i(self.j, dq, N)))
        u_j = Rj.dot(self.u)
        if self.i == 0:
            a_g, u_i = self.a, self.u
            om_i = np.zeros((3, 1))
            Ri = Gi = None
        else:
            pi = p_i(self.i, q, N)
            Ri = R(pi)
            Gi = G(pi)
            a_g = Ri.dot(self.a)
            u_i = Ri.dot(self.u)
            om_i = Ri.dot(2*Gi.dot(dp_i(self.i, dq, N)))

        sin_cz = float(np.cross(u_i.ravel(), u_j.ravel()).dot(a_g.ravel()))
        cos_cz = float(u_i.ravel().dot(u_j.ravel()))
        theta = np.arctan2(sin_cz, cos_cz)
        przekr = self._przekroczenie(theta)

        Qr_i = np.zeros((3, 1))
        Qr_j = np.zeros((3, 1))
        if przekr == 0.0:
            return Qr_i, np.zeros((4, 1)), Qr_j, np.zeros((4, 1))

        dtheta = float(a_g.ravel().dot((om_j - om_i).ravel()))
        tau = -self.k*przekr
        # tlumienie tylko, gdy przegub wchodzi glebiej w ogranicznik
        if przekr > 0 and dtheta > 0:
            tau -= self.c*dtheta
        elif przekr < 0 and dtheta < 0:
            tau -= self.c*dtheta
        M = tau*a_g

        Qp_j = 2*Gj.transpose().dot(Rj.transpose().dot(M))
        if self.i == 0:
            Qp_i = np.zeros((4, 1))
        else:
            Qp_i = -2*Gi.transpose().dot(Ri.transpose().dot(M))
        return Qr_i, Qp_i, Qr_j, Qp_j


class OgranicznikStozka:
    """Miekki ogranicznik stozka zakresu ruchu stawu kulistego (np. bark,
    biodro). Ogranicza kat odchylenia osi `os` czlonu j od tej samej osi
    czlonu i (kierunku neutralnego) do `kat_max`. Powyzej: moment przywraca
    os j do stozka, z tlumieniem tylko przy dalszym odchylaniu. Uzywac ze
    stawem kulistym (Para_Sferyczna + MomentSferyczny).

    Dla i=0 kierunek neutralny jest staly w ukladzie globalnym."""

    def __init__(self, i, j, os, kat_max, k, c):
        self.i = i
        self.j = j
        self.a = np.asarray(os, dtype=float).reshape(3, 1)
        self.a /= np.linalg.norm(self.a)
        self.kat_max = kat_max
        self.k = k
        self.c = c

    def _osie(self, q, N):
        a_j = R(p_i(self.j, q, N)).dot(self.a).ravel()
        a_i = (self.a.ravel() if self.i == 0
               else R(p_i(self.i, q, N)).dot(self.a).ravel())
        return a_i, a_j

    def kat(self, q, N):
        a_i, a_j = self._osie(q, N)
        return float(np.arccos(np.clip(a_i.dot(a_j), -1.0, 1.0)))

    def energia_potencjalna(self, q, N):
        if self.k == 0:
            return 0.0
        nad = max(0.0, self.kat(q, N) - self.kat_max)
        return 0.5*self.k*nad*nad

    def sila(self, q, dq, N):
        Qr_i = np.zeros((3, 1))
        Qr_j = np.zeros((3, 1))
        a_i, a_j = self._osie(q, N)
        kat = float(np.arccos(np.clip(a_i.dot(a_j), -1.0, 1.0)))
        nad = kat - self.kat_max
        if nad <= 0.0:
            return Qr_i, np.zeros((4, 1)), Qr_j, np.zeros((4, 1))

        # os obrotu przywracajaca a_j do a_i (moment na j w te strone)
        n = np.cross(a_j, a_i)
        nn = np.linalg.norm(n)
        if nn < 1e-9:
            return Qr_i, np.zeros((4, 1)), Qr_j, np.zeros((4, 1))
        n = n/nn

        pj = p_i(self.j, q, N)
        Rj = R(pj)
        Gj = G(pj)
        om_j = Rj.dot(2*Gj.dot(dp_i(self.j, dq, N))).ravel()
        if self.i == 0:
            om_i = np.zeros(3)
            Ri = Gi = None
        else:
            pi = p_i(self.i, q, N)
            Ri = R(pi)
            Gi = G(pi)
            om_i = Ri.dot(2*Gi.dot(dp_i(self.i, dq, N))).ravel()

        # tempo dalszego odchylania (skladowa predkosci wzglednej wzdluz -n)
        d_nad = -float((om_j - om_i).dot(n))
        tau = self.k*nad
        if d_nad > 0:
            tau += self.c*d_nad
        M = (tau*n).reshape(3, 1)

        Qp_j = 2*Gj.transpose().dot(Rj.transpose().dot(M))
        if self.i == 0:
            Qp_i = np.zeros((4, 1))
        else:
            Qp_i = -2*Gi.transpose().dot(Ri.transpose().dot(M))
        return Qr_i, Qp_i, Qr_j, Qp_j


class SilaZewn:
    """Sila (Fx/Fy/Fz) lub moment (nx/ny/nz) o stalej wielkosci dzialajacy na czlon."""
    def __init__(self,czlon, rodzaj, wielkosc):
        self.czlon = czlon
        self.rodzaj = rodzaj
        self.wielkosc = wielkosc


