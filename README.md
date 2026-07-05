# uw_dyn: dynamika 3D układów wieloczłonowych

Program do obliczeń dynamiki przestrzennej układów wieloczłonowych (multibody dynamics), napisany w Pythonie w ramach pracy magisterskiej (2016).

Autor: Marcin Szalajski

## Opis metody

Ruch każdego członu (bryły sztywnej) opisany jest 7 współrzędnymi: 3 współrzędnymi położenia środka masy oraz 4 parametrami Eulera (kwaternionem) określającymi orientację. Równania ruchu formułowane są jako układ równań różniczkowo-algebraicznych z mnożnikami Lagrange'a, a więzy stabilizowane są metodą Baumgarte'a (parametry `alfa` i `beta`). Całkowanie po czasie realizuje procedura `sym2`, a wyniki zapisywane są do pliku CSV.

## Zawartość repozytorium

| Plik | Opis |
|------|------|
| `uw_dyn.py` | Główny moduł obliczeniowy: klasy członów, więzów, sił oraz procedury symulacji |
| `lancuch02.py` | Przykład użycia: symulacja łańcucha czterech członów połączonych przegubami obrotowymi |
| `lancuch.blend` | Scena Blendera do wizualizacji ruchu łańcucha |
| `MSzalajski_mgr4.pdf` | Praca magisterska dokumentująca metodę i obliczenia |

## Główne elementy modułu `uw_dyn.py`

- `Uklad`: klasa zbiorcza układu; metody `dodajCzlon`, `dodajWiez`, `dodajWiezD`, `dodajSileWewn`, `dodajSileZewn`, symulacja (`sym`, `sym2`), rozwiązywanie więzów metodą Newtona-Raphsona (`newraph`) oraz zapis wyników (`zapiszWyniki`).
- `Czlon`: człon układu zdefiniowany masą i tensorem bezwładności.
- Więzy kinematyczne (pary kinematyczne):
  - `Para_Sferyczna`: przegub kulisty,
  - `Polaczenie_Obr`: przegub obrotowy,
  - `Polaczenie_Cyl`: para cylindryczna,
  - `Polaczenie_Przes`: para przesuwna,
  - `Para_Prostopadla`, `Para_Prostopadla_D`: więzy prostopadłości wektorów.
- Więzy kierujące: `Odleglosc` (zadana odległość punktów) oraz `Kat` (zadany kąt między wektorami).
- Siły:
  - `SilaWewnProst`: element sprężysto-tłumiący z siłą stałą (sprężyna, tłumik, siłownik) między punktami dwóch członów,
  - `SilaZewn`: siła lub moment zewnętrzny działający na wybrany człon,
  - grawitacja włączana flagą `ukl.grawitacja = True`.

## Przykład użycia

Przykład `lancuch02.py` buduje łańcuch czterech członów zawieszony w punkcie stałym, połączony przegubami obrotowymi, z elementami sprężysto-tłumiącymi i wymuszeniem siłowym:

```python
from uw_dyn import *
import numpy as np

ukl = Uklad()
ukl.dodajCzlon(Czlon(1, 1, J1))
ukl.dodajWiez(Polaczenie_Obr(0, 1, wektor(0,0,0), wektor(0,0,2),
                             wektor(1,0,0), wektor(0,0,1), wektor(0,1,0)))
ukl.grawitacja = True

ukl.sym2(y0, t0, tK, dt, alfa, beta)
ukl.zapiszWyniki('lancuch.csv')
```

Uruchomienie:

```bash
python lancuch02.py
```

Wynikowy plik `lancuch.csv` zawiera przebiegi współrzędnych w czasie i może posłużyć do wizualizacji, na przykład w Blenderze (`lancuch.blend`).

## Wymagania

- Python 3
- NumPy, SciPy, Matplotlib
