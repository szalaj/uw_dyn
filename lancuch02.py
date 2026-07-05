
from uw_dyn import *
import numpy as np

J1 = np.array([[10, 0, 0],[0, 10, 0],[0, 0, 10]])
J2 = np.array([[803.7, 0, 0],[0, 803.7, 0],[0, 0, 368.1]])
J3 = np.array([[136.3, 0, 0],[0, 21140.4, 0],[0, 0, 21097.8]])
J4 = np.array([[39.7, 0, 0],[0, 3545.3, 0],[0, 0, 3537.2]])

ukl = Uklad()
ukl.dodajCzlon(Czlon(1, 1, J1))
ukl.dodajCzlon(Czlon(2, 1, J1))
ukl.dodajCzlon(Czlon(3, 1, J1))
ukl.dodajCzlon(Czlon(4, 1, J1))

ukl.dodajWiez( Polaczenie_Obr(0,1, wektor(0,0,0), wektor(0,0, 2), wektor(1,0,0), wektor(0,0,1), wektor(0,1,0)))
ukl.dodajWiez( Polaczenie_Obr(1,2, wektor(0,0,-2), wektor(0,0, 2), wektor(1,0,0), wektor(0,0,1), wektor(0,1,0)))
ukl.dodajWiez( Polaczenie_Obr(2,3, wektor(0,0,-2), wektor(0,0, 2), wektor(1,0,0), wektor(0,0,1), wektor(0,1,0)))
ukl.dodajWiez( Polaczenie_Obr(3,4, wektor(0,0,-2), wektor(0,0, 2), wektor(1,0,0), wektor(0,0,1), wektor(0,1,0)))


c = 0
f1= 40
f2 =0
k=0

t=1

#ukl.dodajSileWewn(SilaWewnProst(0,4, wektor(0, 0, 120), wektor(0, 0, -2), k, 0, c, 0))

ukl.dodajSileWewn(SilaWewnProst(0,1, wektor(-t, 0, 0), wektor(-t, 0, -2), k, 4, c, f1))
ukl.dodajSileWewn(SilaWewnProst(1,2, wektor(-t, 0, 2), wektor(-t, 0, -2), k, 4, c, 0))
ukl.dodajSileWewn(SilaWewnProst(2,3, wektor(-t, 0, 2), wektor(-t, 0, -2), k, 4, c, 0))
ukl.dodajSileWewn(SilaWewnProst(3,4, wektor(-t, 0, 2), wektor(-t, 0, -2), k, 4, c, 0))

ukl.dodajSileWewn(SilaWewnProst(0,1, wektor(t, 0, 0), wektor(t, 0, -2), k, 4, c, f2))
ukl.dodajSileWewn(SilaWewnProst(1,2, wektor(t, 0, 2), wektor(t, 0, -2), k, 4, c, f2))
ukl.dodajSileWewn(SilaWewnProst(2,3, wektor(t, 0, 2), wektor(t, 0, -2), k, 4, c, f2))
ukl.dodajSileWewn(SilaWewnProst(3,4, wektor(t, 0, 2), wektor(t, 0, -2), k, 4, c, f2))

#ukl.dodajWiezD(Odleglosc(0,1,wektor(0,0,0), wektor(0,0,0),wektor(1,0,0), 0))
#ukl.dodajWiezD(Kat(1,2,wektor(1,0,0),wektor(0,1,0),wektor(1,0,0), -np.pi/2))

#ukl.dodajWiezD(Kat(2,3,wektor(1,0,0),wektor(0,0,1),wektor(1,0,0), -np.pi/4))
#ukl.dodajWiezD(Kat(3,4,wektor(1,0,0),wektor(0,0,1),wektor(1,0,0), np.pi/2))


#ukl.dodajSileZewn(SilaZewn(1,'ny',3000))
# ukl.dodajSileZewn(SilaZewn(3,'Fz', 50))
# ukl.dodajSileZewn(SilaZewn(4,'Fz',-30))
#ukl.dodajSileZewn(SilaZewn(4,'Fx', 2))

ukl.grawitacja = True


q0 = np.zeros(7*ukl.N)

q0[2] = -2
q0[5] = -6
q0[8] = -10
q0[11] = -14

q0[12] = 1
q0[16] = 1
q0[20] = 1
q0[24] = 1


dq0 = np.zeros(7*ukl.N)

#y0 = np.concatenate((q0, dq0),axis =1 )  
y0 = np.concatenate((q0, dq0) )  

t0 = 0
tK = 50
dt = 0.1

alfa = 1
beta = 1


ukl.sym2(y0,t0,tK,dt, alfa, beta)
ukl.zapiszWyniki('lancuch.csv')