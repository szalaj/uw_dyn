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
| `przyklady/piesek.py` | mini-piesek: czworonóg w chodzie pełzającym |
| `przyklady/dron.py` | kwadrokopter z ładunkiem podwieszonym na linie |
| `przyklady/bokser.py` | kickboxing na worku: sierpowy i front kick; worek to wahadło z masą, ciosy trafiają w nie przez kontakt bryła-bryła |
| `przyklady/balans.py` | pełna sylwetka stojąca na stopach z regulatorem balansu PID |
| `przyklady/krok.py` | praca nóg: krok w bok z balansem PID (mechanika + balans; kierunkowy krok w toku) |
| `przyklady/pompka.py` | pompka (push-up): ciało jako dźwignia na palcach stóp, ramiona-zawiasy, dłonie na kontakcie, łokcie napędzane PID |
| `przyklady/lancuch.blend` | scena Blendera do wizualizacji ruchu łańcucha |
| `web/` | wizualizacje Three.js: `przysiad.html`, `robot.html`, `piesek.html`, `dron.html`, `bokser.html`, `balans.html` |
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
uv run python przyklady/piesek.py           # generuje web/dane_piesek.js
uv run python przyklady/dron.py             # generuje web/dane_dron.js
uv run python przyklady/bokser.py           # generuje web/dane_bokser.js
cd web && python3 -m http.server 8000       # potem otworzyć np. localhost:8000/bokser.html
```

Jako zależność w innym projekcie:

```bash
uv add uw-dyn --path ../uw_dyn         # albo: pip install -e ../uw_dyn
```

## Główne elementy pakietu

- `Uklad`: klasa zbiorcza układu; metody `dodajCzlon`, `dodajWiez`, `dodajWiezD`, `dodajSileWewn`, `dodajSileZewn`, symulacja (`sym`, `sym2`), rozwiązywanie więzów metodą Newtona-Raphsona (`newraph`), rzutowanie na więzy (`projekcja_polozen`, `projekcja_predkosci`, przydatne też jako uderzenie plastyczne przy zmianie więzów), energia mechaniczna (`energia`, `energia_kinetyczna`, `energia_potencjalna`) oraz zapis wyników (`zapiszWyniki`).
- `sym2` domyślnie stabilizuje więzy rzutowaniem (dokładne więzy, bez strojenia); klasyczne Baumgarte przez `stabilizacja='baumgarte'` z parametrami `alfa`, `beta`.
- `sym3` to integrator typu RATTLE (Verlet prędkościowy z więzami rozwiązywanymi wewnątrz kroku): rząd 2, więzy na poziomie maszynowym, energia o rzędy wielkości dokładniejsza niż `sym2` przy tym samym koszcie, więc można użyć znacznie większego kroku. Przy sztywnym tłumieniu (`c/m·dt/2 > 1`) lub sztywnych sprężynach (`ω·dt > 2`) jawny `sym3` traci stabilność; tryb `sym3(..., polniejawne=True)` liczy krok półniejawnie, wcielając rzadkie jakobiany sił do układu KKT (tłumienie `C=∂Q/∂v` w kroku prędkości, sztywność `K=∂Q/∂q` w kroku położeń). Jakobiany liczone są lokalnie (per siła, po stopniach swobody jej członów), więc koszt to ~3 rozwiązania układu na krok, niezależnie od `N`, a rząd 2 zostaje zachowany. Dzięki temu pełna sylwetka na stopach liczy się przy `dt=5e-3` (5× większy krok niż w `sym2`) szybciej i dokładniej.
- `Czlon`: człon układu zdefiniowany masą i tensorem bezwładności.
- Więzy kinematyczne (pary kinematyczne):
  - `Para_Sferyczna`: przegub kulisty,
  - `Polaczenie_Obr`: przegub obrotowy,
  - `Polaczenie_Cyl`: para cylindryczna,
  - `Polaczenie_Przes`: para przesuwna,
  - `Para_Prostopadla`, `Para_Prostopadla_D`: więzy prostopadłości wektorów.
- Więzy kierujące: `Odleglosc` (zadana odległość punktów) oraz `Kat` (zadany kąt między wektorami).
- Siły:
  - `SilaWewnProst`: element sprężysto-tłumiący z siłą stałą (sprężyna, tłumik, siłownik) między punktami dwóch członów; z `tylko_rozciaganie=True` działa jak lina,
  - `SilaWPunkcie`: siła zadana w układzie ciała, zaczepiona w punkcie (np. ciąg wirnika drona, podąża za orientacją),
  - `SilaKontaktu`: jednostronny kontakt punktu ciała z podłożem z=0 (model penalty ze sprężyną, tłumieniem i tarciem),
  - `SilaUderzenia`: kontakt penalty bryła-bryła (punkt jednego ciała vs kapsuła drugiego, np. pięść/stopa w worek treningowy); siła równa i przeciwna na oba ciała, zaczepiona w punkcie kontaktu, więc przekazuje pęd i wprawia trafione ciało w ruch,
  - `MomentWzgledny`: aktuator obrotowy 1 DOF w przegubie (regulator PD, opcjonalnie PID przez `ki`) dążący do zadanego kąta, do stawów zawiasowych (kolano, łokieć); opcjonalny limit momentu `moment_max`,
  - `MomentSferyczny`: napędzany staw kulisty 3 DOF (regulator PD/PID na orientacji 3D) dążący do zadanej orientacji względnej; łączyć z `Para_Sferyczna`; opcjonalny limit `moment_max`,
  - `OgranicznikKata`: miękki ogranicznik zakresu ruchu przegubu zawiasowego (łokieć, kolano); jednostronna sprężyna-tłumik aktywna poza `[kat_min, kat_max]`,
  - `OgranicznikStozka`: miękki ogranicznik stożka dla stawu kulistego (bark, biodro); ogranicza odchylenie osi od kierunku neutralnego do `kat_max`,
  - `SilaZewn`: siła lub moment zewnętrzny w osiach globalnych działający na wybrany człon,
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
