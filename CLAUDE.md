# CLAUDE.md

Wskazówki dla Claude Code przy pracy z tym repozytorium.

## Czym jest ten projekt

Solver dynamiki przestrzennej układów wieloczłonowych (multibody dynamics) w Pythonie, stworzony w ramach pracy magisterskiej (2016, Marcin Szalajski). Orientacja członów opisana parametrami Eulera (kwaternionami), równania ruchu z mnożnikami Lagrange'a, stabilizacja więzów metodą Baumgarte'a.

## Struktura

- `uw_dyn.py`: cały silnik obliczeniowy w jednym pliku. Klasy: `Uklad` (układ zbiorczy, symulacja `sym`/`sym2`, Newton-Raphson `newraph`, zapis `zapiszWyniki`), `Czlon`, pary kinematyczne (`Para_Sferyczna`, `Polaczenie_Obr`, `Polaczenie_Cyl`, `Polaczenie_Przes`, `Para_Prostopadla`, `Para_Prostopadla_D`), więzy kierujące (`Odleglosc`, `Kat`), siły (`SilaWewnProst`, `SilaZewn`).
- `lancuch02.py`: przykładowa symulacja łańcucha czterech członów; wynik zapisywany do `lancuch.csv`.
- `lancuch.blend`: scena Blendera do wizualizacji wyników.
- `MSzalajski_mgr4.pdf`: praca magisterska z opisem metody; tu należy szukać teorii i oznaczeń.

## Uruchamianie

```bash
python lancuch02.py
```

Wymaga: NumPy, SciPy, Matplotlib. Brak testów automatycznych i brak systemu budowania.

## Konwencje

- Kod, nazwy klas, metod i komentarze są po polsku (bez polskich znaków w identyfikatorach). Nowy kod pisz w tym samym stylu.
- Wektory kolumnowe tworzą funkcje `wektor(ax, ay, az)` oraz `wektor_p(e0, e1, e2, e3)`; używaj ich zamiast ręcznego budowania tablic NumPy.
- Wektor stanu układu ma długość `14*N`: najpierw `3*N` współrzędnych położeń, potem `4*N` parametrów Eulera, następnie prędkości w tym samym porządku.
- Numeracja członów zaczyna się od 1; indeks 0 oznacza podstawę (układ odniesienia).
- Każda para kinematyczna implementuje metody `wiezyK`, `jakobianK`, `gammaK`; więzy kierujące dodatkowo `wiezyD` i `jakobianD`. Nowe typy więzów muszą zachować ten interfejs.

## Uwagi

- To kod naukowy o wartości historycznej (praca magisterska): przy zmianach zachowuj zgodność wyników z opisem w PDF, nie refaktoryzuj bez wyraźnej prośby.
- Komunikaty commitów po polsku, bez stopki `Co-Authored-By`.
