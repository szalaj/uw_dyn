# CLAUDE.md

Wskazówki dla Claude Code przy pracy z tym repozytorium.

## Czym jest ten projekt

Biblioteka Pythona do dynamiki przestrzennej układów wieloczłonowych (multibody dynamics), wywodząca się z pracy magisterskiej (2016, Marcin Szalajski). Orientacja członów opisana parametrami Eulera (kwaternionami), równania ruchu z mnożnikami Lagrange'a, stabilizacja więzów metodą Baumgarte'a.

**Najpierw przeczytaj `PLAN.md`**: tam jest mapa drogowa, podjęte decyzje architektoniczne i aktualny stan prac. Najważniejsza decyzja: rozwijamy wszystko w Pythonie; port do Rust jest warunkowy i wchodzi w grę dopiero, gdyby wydajności nie dało się uratować optymalizacją.

## Struktura

- `src/uw_dyn/dynamika.py`: cały silnik obliczeniowy. Klasy: `Uklad` (układ zbiorczy, symulacja `sym`/`sym2`, Newton-Raphson `newraph`, zapis `zapiszWyniki`), `Czlon`, pary kinematyczne (`Para_Sferyczna`, `Polaczenie_Obr`, `Polaczenie_Cyl`, `Polaczenie_Przes`, `Para_Prostopadla`, `Para_Prostopadla_D`), więzy kierujące (`Odleglosc`, `Kat`), siły (`SilaWewnProst`, `SilaZewn`).
- `src/uw_dyn/__init__.py`: publiczne API pakietu (jawna lista `__all__`).
- `tests/`: pytest; `conftest.py` zawiera budowę wahadła testowego i obliczanie energii mechanicznej.
- `przyklady/lancuch02.py`: przykładowa symulacja łańcucha czterech członów; wynik zapisywany do `lancuch.csv` (ignorowany przez git).
- `przyklady/lancuch.blend`: scena Blendera do wizualizacji wyników.
- `docs/MSzalajski_mgr4.pdf`: praca magisterska; tu należy szukać teorii i oznaczeń.

## Uruchamianie

```bash
uv sync                                # instalacja
uv run pytest                          # testy (~6 s, wszystkie muszą przechodzić)
uv run python przyklady/lancuch02.py   # przykład
```

## Konwencje

- Kod, nazwy klas, metod, testów i komentarze są po polsku (bez polskich znaków w identyfikatorach). Nowy kod pisz w tym samym stylu.
- Wektory kolumnowe tworzą funkcje `wektor(ax, ay, az)` oraz `wektor_p(e0, e1, e2, e3)`; używaj ich zamiast ręcznego budowania tablic NumPy.
- Wektor stanu układu ma długość `14*N`: najpierw `3*N` współrzędnych położeń, potem `4*N` parametrów Eulera, następnie prędkości w tym samym porządku.
- Numeracja członów zaczyna się od 1; indeks 0 oznacza podstawę (układ odniesienia).
- Każda para kinematyczna implementuje metody `wiezyK`, `jakobianK`, `gammaK`; więzy kierujące dodatkowo `wiezyD` i `jakobianD`. Nowe typy więzów muszą zachować ten interfejs.

## Pułapki

- Nowy NumPy nie pozwala na `float()` z macierzy (1,1): używaj `.item()`.
- Więzy kierujące (`Kat`, `Odleglosc`, dodawane przez `dodajWiezD`) działają tylko w `newraph` (obliczanie warunków początkowych); w samej dynamice nie są egzekwowane.
- `newraph` wymaga, żeby liczba więzów (kinematyczne + normy kwaternionów + kierujące) była równa `7*N`.
- `Uklad.jakobianK` jest memoizowany po `q` (i unieważniany w `dodajCzlon`/`dodajWiez`); jeśli dodajesz nowe metody mutujące układ, unieważnij `self._jakK_klucz`.

## Uwagi

- Zmiany w fizyce muszą zachowywać zgodność z opisem w PDF i przechodzić testy walidacyjne (okres wahadła, energia, więzy); nie refaktoryzuj silnika bez wyraźnej prośby.
- Komunikaty commitów po polsku, bez stopki `Co-Authored-By`.
- Po każdym istotnym kroku zaktualizuj status w `PLAN.md`.
