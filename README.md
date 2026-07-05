# uw_dyn: dynamika 3D układów wieloczłonowych

Biblioteka Pythona do obliczeń dynamiki przestrzennej układów wieloczłonowych
(multibody dynamics). Powstała w ramach pracy magisterskiej (2016), obecnie
rozwijana jako pakiet wielokrotnego użytku. Autor: Marcin Szalajski.

## Opis metody

Ruch każdego członu (bryły sztywnej) opisany jest 7 współrzędnymi: 3 współrzędnymi położenia środka masy oraz 4 parametrami Eulera (kwaternionem) określającymi orientację. Równania ruchu formułowane są jako układ równań różniczkowo-algebraicznych z mnożnikami Lagrange'a, a więzy stabilizowane są metodą Baumgarte'a (parametry `alfa` i `beta`). Całkowanie po czasie realizuje procedura `sym2`, a wyniki zapisywane są do pliku CSV.

Teoria i oznaczenia: `docs/MSzalajski_mgr4.pdf`.

## Struktura repozytorium

| Ścieżka | Opis |
|------|------|
| `src/uw_dyn/` | pakiet: moduły `algebra`, `czlony`, `wiezy`, `sily`, `uklad` |
| `tests/` | testy pytest: algebra, walidacja fizyczna, rzutowanie, energia |
| `przyklady/lancuch02.py` | łańcuch czterech członów (wynik do CSV) |
| `przyklady/przysiad.py` | staw kolanowy podczas przysiadu (mięśnie sprężysto-tłumiące) |
| `przyklady/robot_kroczacy.py` | najprostszy robot kroczący (chód cyrklowy) |
| `przyklady/lancuch.blend` | scena Blendera do wizualizacji ruchu łańcucha |
| `web/` | wizualizacje Three.js: `przysiad.html`, `robot.html` |
| `docs/MSzalajski_mgr4.pdf` | praca magisterska dokumentująca metodę i obliczenia |
| `PLAN.md` | mapa drogowa rozwoju i stan prac |

## Instalacja i uruchomienie

Projekt używa [uv](https://docs.astral.sh/uv/):

```bash
uv sync                                # instalacja pakietu i zależności
uv run pytest                          # testy
uv run python przyklady/lancuch02.py   # przykładowa symulacja -> lancuch.csv

# przykłady z wizualizacją web (Three.js, wymagany internet dla CDN):
uv run python przyklady/przysiad.py         # generuje web/dane_przysiad.js
uv run python przyklady/robot_kroczacy.py   # generuje web/dane_robot.js
cd web && python3 -m http.server 8000       # potem otworzyć np. localhost:8000/przysiad.html
```

Jako zależność w innym projekcie:

```bash
uv add uw-dyn --path ../uw_dyn         # albo: pip install -e ../uw_dyn
```

## Główne elementy pakietu

- `Uklad`: klasa zbiorcza układu; metody `dodajCzlon`, `dodajWiez`, `dodajWiezD`, `dodajSileWewn`, `dodajSileZewn`, symulacja (`sym`, `sym2`), rozwiązywanie więzów metodą Newtona-Raphsona (`newraph`), rzutowanie na więzy (`projekcja_polozen`, `projekcja_predkosci`, przydatne też jako uderzenie plastyczne przy zmianie więzów), energia mechaniczna (`energia`, `energia_kinetyczna`, `energia_potencjalna`) oraz zapis wyników (`zapiszWyniki`).
- `sym2` domyślnie stabilizuje więzy rzutowaniem (dokładne więzy, bez strojenia); klasyczne Baumgarte przez `stabilizacja='baumgarte'` z parametrami `alfa`, `beta`.
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

```python
from uw_dyn import *
import numpy as np

ukl = Uklad()
ukl.dodajCzlon(Czlon(1, 1, np.diag([10., 10., 10.])))
ukl.dodajWiez(Polaczenie_Obr(0, 1, wektor(0,0,0), wektor(0,0,2),
                             wektor(1,0,0), wektor(0,0,1), wektor(0,1,0)))
ukl.grawitacja = True

q0 = np.zeros(7); q0[2] = -2; q0[3] = 1
y0 = np.concatenate((q0, np.zeros(7)))

ukl.sym2(y0, t0=0, tK=50, dt=0.01)
ukl.zapiszWyniki('wyniki.csv')
```

Wektor stanu ma długość `14*N`: najpierw `3*N` współrzędnych położeń, potem `4*N` parametrów Eulera, następnie prędkości w tym samym porządku. Wyniki (CSV, wiersz na krok czasowy) mogą posłużyć do wizualizacji, na przykład w Blenderze.

## Testy

Testy walidacyjne sprawdzają fizykę względem rozwiązań analitycznych:

- okres małych drgań wahadła fizycznego zgodny ze wzorem `T = 2π·sqrt(I/(mgL))` z dokładnością 2%,
- zachowanie energii mechanicznej bez tłumienia,
- spełnienie więzów kinematycznych i normy kwaternionów w trakcie ruchu,
- własności macierzy obrotu i funkcji kwaternionowych,
- regresja: stabilność symulacji łańcucha 4 członów.

## Kierunek rozwoju

Zobacz `PLAN.md`. W skrócie: całość rozwijana w Pythonie (fizyka, testy, optymalizacja, wizualizacja web z Three.js). Port rdzenia do Rust (PyO3, WASM) jest przewidziany warunkowo, tylko gdyby wydajność wersji Python przestała wystarczać w realnych zastosowaniach.
