# -*- coding: utf-8 -*-
# uw_dyn: dynamika 3D ukladow wieloczlonowych
# autor: Marcin Szalajski (praca magisterska, 2016; pakiet od 2026)

"""Uklad wieloczlonowy: skladanie rownan ruchu i symulacja."""

from scipy.integrate import solve_ivp
import scipy.linalg
import numpy as np
import functools

from uw_dyn.algebra import (p_i, dp_i, G, dG, wektor_p,
                            jakobian_p_kolumny, jakobian_r_kolumny)
from uw_dyn.czlony import Czlon
from uw_dyn.wiezy import Polaczenie

# przyspieszenie ziemskie [m/s^2]
GRAWITACJA = 9.80665


class Uklad:
    """Uklad wieloczlonowy: czlony, wiezy i sily oraz procedury symulacji (sym, sym2), rozwiazywania warunkow poczatkowych (newraph) i zapisu wynikow."""

    def __init__(self):
        self.czlony = []
        self.wiezy_k = [] # wiezy kinematyczne
        self.wiezy_d = [] # wiezy kierujace
        self.silyWewn = []
        self.silyZewn = []
        self.N = 0 #ilosc czlonow
        self.M = 0 #ilosc wiezow kinematycznych
        self.Mi = 0 #ilosc wiezow kierujacych (war.poczatkowe)
        self.Y = [] #wyniki symulacji
        self.grawitacja = True
        # pamiec podreczna (przyspieszenie obliczen)
        self._zbM = None      # macierz masowa (stala)
        self._zbJ = None      # tensor bezwladnosci (staly)
        self._jakK_klucz = None  # q, dla ktorego policzono jakobianK
        self._jakK = None


    def dodajCzlon(self, czlon):
        #sprawdzenie czy 'czlon' jest typu Czlon
        if isinstance(czlon, Czlon):
            self.czlony.append(czlon)
            self.N += 1
            self._zbM = None
            self._zbJ = None
            self._jakK_klucz = None
        else:
            raise Exception("Obiekt nie jest typu Czlon")

    def dodajWiez(self, wiez):
        if isinstance(wiez, Polaczenie):
            self.wiezy_k.append(wiez)
            self.M += wiez.m
            self._jakK_klucz = None
        else:
            raise Exception('Obiekt nie jest typu Polaczenie')
            
    def dodajWiezD(self, wiez):
        if isinstance(wiez, Polaczenie):
            self.wiezy_d.append(wiez)
            self.Mi += wiez.m
        else:
            raise Exception('Obiekt nie jest typu Polaczenie')        
            
    def dodajSileWewn(self, sila):
        self.silyWewn.append(sila)
        
    def dodajSileZewn(self, sila):
        self.silyZewn.append(sila)
            

    # wiezy kinematyczne zbiorczo
    def wiezyK(self, q):
    
        #
        F= np.zeros([self.M,1])
        k=0
        for w in self.wiezy_k:
            F[k:k+w.m,:] = w.wiezyK(q,self.N)
            k=k+w.m    

        return F
    
    # wiezy kierujace zbiorczo    
    def wiezyD(self, q):
    
        #
        F= np.zeros([self.Mi,1])
        k=0
        for w in self.wiezy_d:
            F[k:k+w.m,:] = w.wiezyD(q,self.N)
            k=k+w.m    

        return F
        
    
    # wiezy parametrow eulera zbiorczo
    def wiezyP(self, q):
        
        #liczba czlonow w ukladzie
        N = self.N
        
        Fp = np.zeros([N,1])
        
        for k in range(0,N):
            pi = p_i(k+1,q,N)
            Fp[k,:] = pi.transpose().dot(pi) - 1
    
        return Fp
        
    # wiezy kinematyczne i parametrow Eulera zbiorczo   
    def wiezyKP(self, q):
        
        wiezy_K = self.wiezyK(q)
        wiezy_P = self.wiezyP(q)

        F = np.concatenate((wiezy_K, wiezy_P), axis=0)
        return F
    
    # wiezy kinematyczne, kierujace i parametrow Eulera zbiorczo       
    def wiezyKPD(self, q):
        
        wiezy_K = self.wiezyK(q)
        wiezy_P = self.wiezyP(q)
        wiezy_D = self.wiezyD(q)
        
        F = np.concatenate((wiezy_K, wiezy_P, wiezy_D), axis=0)
        return F
        
    # jakobian wiezow kinematycznych zbiorczo
    # (memoizacja: Lstrona i Pstrona pytaja o ten sam q w jednym kroku)
    def jakobianK(self, q):

        M=self.M
        N=self.N

        klucz = q.tobytes() if isinstance(q, np.ndarray) else None
        if klucz is not None and klucz == self._jakK_klucz:
            return self._jakK

        Fq=np.zeros([M,7*N])

        k=0
        for w in self.wiezy_k:
            Fqi, Fqj = w.jakobianK(q,N)
            if w.i == 0:
                Fq[k:k+w.m, 7*(w.j-1):7*(w.j-1)+7]=Fqj
            else:
                Fq[k:k+w.m, 7*(w.i-1):7*(w.i-1)+7]=Fqi
                Fq[k:k+w.m, 7*(w.j-1):7*(w.j-1)+7]=Fqj
            k+=w.m

        self._jakK_klucz = klucz
        self._jakK = Fq

        return Fq
     
    # jakobian wiezow kierujacych zbiorczo
    def jakobianD(self, q):
        
        Mi=self.Mi
        N=self.N
        
        Fq=np.zeros([Mi,7*N])
        
        k=0
        for w in self.wiezy_d:
            Fqi, Fqj = w.jakobianD(q,N)
            if w.i == 0:
                Fq[k:k+w.m, 7*(w.j-1):7*(w.j-1)+7]=Fqj
            else:
                Fq[k:k+w.m, 7*(w.i-1):7*(w.i-1)+7]=Fqi
                Fq[k:k+w.m, 7*(w.j-1):7*(w.j-1)+7]=Fqj
            k+=w.m
        
    
        return Fq
        
    # jakobian wiezow parametrow eulera zbiorczo
    def jakobianP(self, q):
        N=self.N
        Fq_p=np.zeros([N,7*N])
        for k in range(0,N):
            pi=q[3*N+4*k:3*N+4*k+4]
            #pi=p_i(k,q)
            Fq_p[k,k*7+3:k*7+7]=2*pi
            
        return Fq_p
    
    # jakobian wiezow kinematycznych i parametrow eulera zbiorczo
    def jakobianKP(self, q):
        
        jak_K = self.jakobianK(q)
        jak_P = self.jakobianP(q)
        
        jak = np.concatenate((jak_K, jak_P), axis=0)
    
        return jak
        
    # jakobian wiezow kin., kier. i parametrow eulera zbiorczo    
    def jakobianKPD(self, q):
        
        jak_K = self.jakobianK(q)
        jak_P = self.jakobianP(q)
        jak_D = self.jakobianD(q)
        
        jak = np.concatenate((jak_K, jak_P, jak_D), axis=0)
    
        return jak
       
    # wektor gamma wiezow kin. zbiorczo
    def gammaK(self, q,dq):
        
        M = self.M
    
        gamK= np.zeros([M,1])
        k=0
        for w in self.wiezy_k:
            gamK[k:k+w.m,:] = w.gammaK(q,dq,self.N)
            k=k+w.m
            
        return gamK
        
    # wektor gamma wiezow par. eulera zbiorczo
    def gammaP(self, dq):
        
        N=self.N
        
        gammP = np.zeros([N,1])
        
        for k in range(0,N):
            dpi=dq[3*N+4*k:3*N+4*k+4]
            gammP[k,:]=-2*dpi.transpose().dot(dpi)
    
        return gammP
    
    # macierz masowa zbiorczo (stala; liczona raz)
    def zbiorczeM(self):

        if self._zbM is None:
            czM =[]
            for cz in self.czlony:
                czM.append(cz.M())
            self._zbM = functools.reduce(scipy.linalg.block_diag, czM)

        return self._zbM

    # tensor bezwladnosci zbiorczo (staly; liczony raz)
    def zbiorczeJ(self):

        if self._zbJ is None:
            czJ =[]
            for cz in self.czlony:
                czJ.append(cz.J)
            self._zbJ = functools.reduce(scipy.linalg.block_diag, czJ)

        return self._zbJ

    # macierz G zbiorczo
    def zbiorczeG(self,q):

        N = self.N

        zbG = np.zeros([3*N, 4*N])
        for k in range(0,N):
            pi=q[3*N+4*k:3*N+4*k+4]
            pi = wektor_p(pi[0],pi[1],pi[2],pi[3])
            zbG[3*k:3*k+3, 4*k:4*k+4] = G(pi)

        return zbG

    # macierz dG zbiorczo
    def zbiorcze_dG(self,dq):
        N = self.N

        zbdG = np.zeros([3*N, 4*N])
        for k in range(0,N):
            dpi=dq[3*N+4*k:3*N+4*k+4]
            dpi = wektor_p(dpi[0],dpi[1],dpi[2],dpi[3])
            zbdG[3*k:3*k+3, 4*k:4*k+4] = dG(dpi)

        return zbdG
        

    # kolumny jakobianu odpowiadajace r (do wybierania podmacierzy)
    def _kolumny_r(self):
        return np.fromiter(jakobian_r_kolumny(self.N), dtype=int)

    # kolumny jakobianu odpowiadajace p
    def _kolumny_p(self):
        return np.fromiter(jakobian_p_kolumny(self.N), dtype=int)

    # jakobianK polozen zbiorczo
    def zbiorczeF_r(self,q):

        F_q = self.jakobianK(q)
        #wybranie kolumn odpowiadajacych r
        F_r = F_q[:, self._kolumny_r()]
        return F_r

    # jakobianK par. eulera zbiorczo
    def zbiorczeF_p(self,q):

        F_q = self.jakobianK(q)
        #wybranie kolumn odpowiadajacych p
        F_p = F_q[:, self._kolumny_p()]
        return F_p

    # jakobianP par. eulera zbiorczo
    def zbiorczeFp_p(self,q):

        Fp = self.jakobianP(q)
        #wybranie kolumn odpowiadajacych p
        Fp_p = Fp[:, self._kolumny_p()]
        return Fp_p
        
    # parametry eulera zbiorczo    
    def zbiorcze_p(self,q):
        N=self.N
        p=[]
        for k in range(0,N):
            pi=q[3*N+4*k:3*N+4*k+4]
            #pi = p_i(k,dq)
            p.append(pi[0])
            p.append(pi[1])
            p.append(pi[2])
            p.append(pi[3])
            
        p = np.array(p)
        #zamiana z array 1D na 2D 
        return p.reshape((-1,1))
    
    # lewa strona rownania dynamiki
    def Lstrona(self,q,dq):
        N = self.N
        M = self.M
        #M - masa , może sie mylic z liczba wiezow kinematycznych M
        zbM = self.zbiorczeM()
        J = self.zbiorczeJ()
        G = self.zbiorczeG(q)
        F_r = self.zbiorczeF_r(q)
        F_p = self.zbiorczeF_p(q)
        Fp_p = self.zbiorczeFp_p(q)
        
        GT = G.transpose()
        Iloczyn = 4*GT.dot(J).dot(G)

        # zlozenie duzej macierzy lewej strony (prealokacja blokow)
        n = 7*N + M + N
        Lstr = np.zeros([n, n])
        Lstr[0:3*N, 0:3*N] = zbM
        Lstr[0:3*N, 7*N:7*N+M] = F_r.transpose()
        Lstr[3*N:7*N, 3*N:7*N] = Iloczyn
        Lstr[3*N:7*N, 7*N:7*N+M] = F_p.transpose()
        Lstr[3*N:7*N, 7*N+M:n] = Fp_p.transpose()
        Lstr[7*N:7*N+M, 0:3*N] = F_r
        Lstr[7*N:7*N+M, 3*N:7*N] = F_p
        Lstr[7*N+M:n, 3*N:7*N] = Fp_p

        return Lstr
        
        
    # prawa strona rownania dynamiki
    def Pstrona(self,q,dq):
        N = self.N
        
        G = self.zbiorczeG(q)
        dG = self.zbiorcze_dG(dq)
        dGT = dG.transpose()
        
        J = self.zbiorczeJ()
        
        
        FA = np.zeros([3*N,1])
        nA = np.zeros([3*N,1]) 

        # dodanie sil zewnetrznych jesli istnieja
        if len(self.silyZewn):
            for s in self.silyZewn:
                i = s.czlon-1
                if s.rodzaj == 'Fx':
                    FA[3*i,:] = s.wielkosc
                elif s.rodzaj == 'Fy':
                    FA[3*i+1,:] = s.wielkosc
                elif s.rodzaj == 'Fz':
                    FA[3*i+2,:] = s.wielkosc
                elif s.rodzaj == 'nx':
                    nA[3*i,:] = s.wielkosc
                elif s.rodzaj == 'ny':
                    nA[3*i+1,:] = s.wielkosc
                elif s.rodzaj == 'nz':
                    nA[3*i+2,:] = s.wielkosc
                else:
                    raise Exception('zla sila')
                    
        #dodanie sily grawitacyjnej
        if self.grawitacja:
            for cz in self.czlony:
                FA[3*(cz.i-1)+2,:] += -GRAWITACJA*cz.m
        
        p = self.zbiorcze_p(q)
                                        
        Ilo = 2*G.transpose().dot(nA) + 8*dGT.dot(J).dot(dG).dot(p)
        
        #jesli sa dodane jakies sily wewn
        if len(self.silyWewn):
            for s in self.silyWewn:
                if s.i == 0:
                    Qr_i, Qp_i, Qr_j, Qp_j = s.sila(q,dq,self.N)
                    #print(Qr_j)
                    FA[3*(s.j-1):3*(s.j-1)+3] += Qr_j
                    Ilo[4*(s.j-1):4*(s.j-1)+4] += Qp_j
                else:
                   Qr_i, Qp_i, Qr_j, Qp_j = s.sila(q,dq,self.N)
                   FA[3*(s.i-1):3*(s.i-1)+3] += Qr_i
                   Ilo[4*(s.i-1):4*(s.i-1)+4] += Qp_i
                   FA[3*(s.j-1):3*(s.j-1)+3] += Qr_j
                   Ilo[4*(s.j-1):4*(s.j-1)+4] += Qp_j
                   

        
        # zamiana kolejnosci wspolrzednych dopasowanych do jakobianu
        # q-> r1 r2 rn p1 p2 pn. bedzie r1 p1 r2 p2 rn pn
        dq_jak = np.concatenate((np.asarray(dq[0:3*N]).reshape(N,3),
                                 np.asarray(dq[3*N:7*N]).reshape(N,4)),
                                axis=1).reshape(-1,1)


        Ps = np.concatenate((self.gammaK(q,dq), self.gammaP(dq)), axis=0)
        
        #obliczenie wspolczynnikow metody Baugmarte'a
        Baug = -2*self.alfa*self.jakobianKP(q).dot(dq_jak) - self.beta*self.beta*self.wiezyKP(q)

        Ps += Baug
    
        Pstr = np.concatenate((FA, Ilo, Ps), axis=0)
    
        return Pstr
            
    # jakobian wiezow KP z kolumnami w porzadku wektora q (r..., p...)
    def _jakobianKP_q(self, q):
        Fq = self.jakobianKP(q)
        return np.hstack((Fq[:, self._kolumny_r()], Fq[:, self._kolumny_p()]))

    # normalizacja kwaternionow wszystkich czlonow w wektorze q
    def _normalizuj_kwaterniony(self, q):
        N = self.N
        for k in range(0, N):
            p = q[3*N+4*k:3*N+4*k+4]
            p /= np.linalg.norm(p)

    # rzutowanie polozen na rozmaitosc wiezow
    def projekcja_polozen(self, q, tol=1e-10, maks_iter=20):
        """Poprawia q tak, by wiezy kinematyczne i normy kwaternionow byly
        spelnione (Newton z poprawka minimalnej normy). Zwraca nowe q.

        Kwaterniony sa normalizowane przed i po kazdej iteracji: jakobiany
        wiezow sa dokladne tylko dla kwaternionow jednostkowych, wiec
        poprawka nie moze zawierac skladowej skalujacej kwaternion."""
        q = np.array(q, dtype=float)
        self._normalizuj_kwaterniony(q)
        F = self.wiezyKP(q)
        norma = np.linalg.norm(F)
        for _ in range(maks_iter):
            if norma < tol:
                break
            Jq = self._jakobianKP_q(q)
            popr = Jq.T.dot(np.linalg.solve(Jq.dot(Jq.T), F)).ravel()

            # kontrola kroku (tlumienie w razie rozbiegania)
            wsp = 1.0
            for _ in range(8):
                q_nowe = q - wsp*popr
                self._normalizuj_kwaterniony(q_nowe)
                F_nowe = self.wiezyKP(q_nowe)
                norma_nowa = np.linalg.norm(F_nowe)
                if norma_nowa < norma:
                    break
                wsp *= 0.5
            else:
                break

            q, F, norma = q_nowe, F_nowe, norma_nowa
        return q

    # rzutowanie predkosci na wiezy (zgodne z macierza mas)
    def projekcja_predkosci(self, q, dq):
        """Najblizsze (w metryce macierzy mas) predkosci spelniajace J dq = 0.

        Fizycznie odpowiada uderzeniu plastycznemu: przy naglej zmianie
        wiezow (np. postawienie stopy robota) daje predkosci po uderzeniu."""
        N = self.N
        M = self.M
        q = np.asarray(q, dtype=float)
        dq = np.asarray(dq, dtype=float)

        Gz = self.zbiorczeG(q)
        Mq = np.zeros([7*N, 7*N])
        Mq[0:3*N, 0:3*N] = self.zbiorczeM()
        Mq[3*N:, 3*N:] = 4*Gz.T.dot(self.zbiorczeJ()).dot(Gz)

        Jq = self._jakobianKP_q(q)
        nw = 7*N + M + N
        A = np.zeros([nw, nw])
        b = np.zeros(nw)
        A[0:7*N, 0:7*N] = Mq
        A[0:7*N, 7*N:] = Jq.T
        A[7*N:, 0:7*N] = Jq
        b[0:7*N] = Mq.dot(dq)

        return np.linalg.solve(A, b)[0:7*N]

    # integracja poljawnym schematem Eulera (symplektycznym) ze stalym krokiem
    def sym2(self, y0,t0,tK,dt, alfa=0, beta=0, stabilizacja='rzutowanie'):
        """Symulacja ukladu poljawnym schematem Eulera ze stalym krokiem dt.

        stabilizacja='rzutowanie' (domyslna): po kazdym kroku polozenia sa
        rzutowane na rozmaitosc wiezow, a predkosci na wiezy predkosciowe;
        alfa i beta sa wtedy pomijane. stabilizacja='baumgarte': klasyczna
        stabilizacja Baumgarte'a z parametrami alfa, beta i normalizacja
        kwaternionow. y0 nie jest modyfikowane. Wyniki w self.Y."""
        N = self.N
        rzutuj = (stabilizacja == 'rzutowanie')
        if stabilizacja not in ('rzutowanie', 'baumgarte'):
            raise ValueError("stabilizacja: 'rzutowanie' albo 'baumgarte'")
        self.alfa = 0 if rzutuj else alfa
        self.beta = 0 if rzutuj else beta
        y = np.array(y0, dtype=float).copy()
        q = y[0:7*N]
        dq = y[7*N:14*N]
        wyniki = [y.copy()]

        for t in np.arange(t0,tK,dt):
            LS = self.Lstrona(q,dq)
            PS = self.Pstrona(q,dq)
            ddq = np.linalg.solve(LS, PS).ravel()

            dq += ddq[0:7*N]*dt
            q += dq*dt

            if not np.isfinite(q).all() or np.abs(q).max() > 1e6:
                raise RuntimeError(
                    f'symulacja rozbiegla sie w t={t:.4f}: '
                    'zmniejsz krok dt albo sprawdz model')

            if rzutuj:
                q[:] = self.projekcja_polozen(q)
                dq[:] = self.projekcja_predkosci(q, dq)
            else:
                # normalizacja kwaternionow po kroku calkowania
                for k in range(0,N):
                    p = q[3*N+4*k:3*N+4*k+4]
                    p /= np.linalg.norm(p)

            wyniki.append(y.copy())

        self.Y = np.array(wyniki)


    # integracja adaptacyjna (scipy solve_ivp, RK45)
    def sym(self,y0,t0,tK,dt,alfa,beta, rtol=1e-8, atol=1e-8):
        """Symulacja ukladu metoda adaptacyjna RK45 (scipy.solve_ivp).

        Dokladniejsza (ale wolniejsza) niz sym2; dt okresla tylko gestosc
        zapisu wynikow. y0 nie jest modyfikowane. Wyniki w self.Y."""
        N = self.N
        self.alfa = alfa
        self.beta = beta

        def prawa(t,y):
            q = y[0:7*N]
            dq = y[7*N:14*N]

            LS = self.Lstrona(q,dq)
            PS = self.Pstrona(q,dq)

            ddq = np.linalg.solve(LS, PS).ravel()[0:7*N]
            return np.concatenate((dq, ddq))

        t_eval = np.arange(t0, tK+dt/2, dt)
        roz = solve_ivp(prawa, (t0, tK), np.array(y0, dtype=float),
                        t_eval=t_eval, method='RK45', rtol=rtol, atol=atol)
        if not roz.success:
            raise RuntimeError('sym: calkowanie nieudane: ' + roz.message)

        self.Y = roz.y.T.copy()

    # rozwiazanie nieliniowego ukladu rownan
    # metoda Newtona-Raphsona
    # dla obliczenia warunkow pocz.
    def newraph(self,q0):
        q=q0
        F=self.wiezyKPD(q)
        it=1
        norma =[]
        norma.append(np.linalg.norm(F))
        while (np.linalg.norm(F)>0.001 and it<25):
            F=self.wiezyKPD(q)
            J=self.jakobianKPD(q)
            
            dq = np.linalg.solve(J,F).ravel()
            
            r = []
            p = []
            for i in range(0, self.N):
            
                ri = dq[7*i:7*i+3]
                pi = dq[7*i+3:7*i+7]
                r.append(ri[0])
                r.append(ri[1])
                r.append(ri[2])
                
                p.append(pi[0])
                p.append(pi[1])
                p.append(pi[2])
                p.append(pi[3])
        
            dq2 = np.concatenate((r,p))
            
            q=q-dq2
            norma.append(np.linalg.norm(F))
            it=it+1  
            if(it==24):
                print('cos nie tak')  
                
        print('warunki początkowe obliczone')
        return q        


    def zapiszWyniki(self, nazwaPliku):
        if len(self.Y) == 0:
            raise Exception('brak wynikow symulacji, najpierw zasymuluj uklad')
        else:
            np.savetxt(nazwaPliku, self.Y, delimiter=";")

    # --- energia mechaniczna ---

    def energia_kinetyczna(self, y):
        """Energia kinetyczna ukladu dla wektora stanu y = [q, dq]."""
        N = self.N
        y = np.asarray(y, dtype=float)
        q = y[0:7*N]
        dq = y[7*N:14*N]
        E = 0.0
        for cz in self.czlony:
            i = cz.i
            dr = dq[3*(i-1):3*(i-1)+3]
            p = p_i(i, q, N)
            dp = dp_i(i, dq, N)
            om = 2*G(p).dot(dp)  # predkosc katowa w ukladzie ciala
            E += 0.5*cz.m*float(dr.dot(dr))
            E += 0.5*om.T.dot(cz.J).dot(om).item()
        return E

    def energia_potencjalna(self, y):
        """Energia potencjalna: grawitacja (gdy wlaczona) + sprezyny
        elementow SilaWewnProst (czlon tlumika i sily stalej pominiety)."""
        N = self.N
        y = np.asarray(y, dtype=float)
        q = y[0:7*N]
        E = 0.0
        if self.grawitacja:
            for cz in self.czlony:
                E += cz.m*GRAWITACJA*float(q[3*(cz.i-1)+2])
        for s in self.silyWewn:
            if s.k != 0:
                l = s.dlugosc(q, N)
                E += 0.5*s.k*(l - s.l0)**2
        return E

    def energia(self, y):
        """Energia mechaniczna (kinetyczna + potencjalna) dla stanu y."""
        return self.energia_kinetyczna(y) + self.energia_potencjalna(y)     
    

        





        