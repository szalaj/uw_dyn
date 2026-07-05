# CLAUDE.md

Wskazówki dla Claude Code przy pracy z tym repozytorium.

## Czym jest ten projekt

Biblioteka Pythona do dynamiki przestrzennej układów wieloczłonowych (multibody dynamics), wywodząca się z pracy magisterskiej (2016, Marcin Szalajski). Orientacja członów opisana parametrami Eulera (kwaternionami), równania ruchu z mnożnikami Lagrange'a, stabilizacja więzów metodą Baumgarte'a.

**Najpierw przeczytaj `PLAN.md`**: tam jest mapa drogowa, podjęte decyzje architektoniczne i aktualny stan prac. Najważniejsze decyzje: rozwijamy wszystko w Pythonie (port do Rust tylko warunkowo, gdyby optymalizacja nie wystarczyła); bibliotekę rozbudowujemy i weryfikujemy na trzech przykładach kanonicznych: przysiad, robot kroczący (docelowo mini-piesek) i kwadrokopter.

## Struktura

- `src/uw_dyn/algebra.py`: wektory, kwaterniony (parametry Eulera), macierze R/G, skew.
- `src/uw_dyn/czlony.py`: `Czlon`.
- `src/uw_dyn/wiezy.py`: pary kinematyczne (`Para_Sferyczna`, `Polaczenie_Obr`, `Polaczenie_Cyl`, `Polaczenie_Przes`, `Para_Prostopadla`, `Para_Prostopadla_D`) i więzy kierujące (`Odleglosc`, `Kat`).
- `src/uw_dyn/sily.py`: `SilaWewnProst` (sprężyna/tłumik/siła stała), `SilaZewn`.
- `src/uw_dyn/uklad.py`: `Uklad` (składanie równań ruchu, symulacja `sym`/`sym2`, `newraph`, rzutowania `projekcja_polozen`/`projekcja_predkosci`, metody energii, `zapiszWyniki`).
- `src/uw_dyn/dynamika.py`: alias zgodności wstecznej (re-eksport wszystkiego).
- `src/uw_dyn/__init__.py`: publiczne API pakietu (jawna lista `__all__`).
- `tests/`: pytest; `conftest.py` zawiera budowę wahadła testowego i obliczanie energii mechanicznej.
- `przyklady/`: `lancuch02.py` (CSV), `przysiad.py`, `robot_kroczacy.py`; `lancuch.blend` to scena Blendera.
- `web/`: wizualizacje Three.js (`przysiad.html`, `robot.html`); pliki `dane_*.js` generują skrypty z `przyklady/`.
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
- Jakobiany więzów są dokładne tylko dla kwaternionów jednostkowych (macierz `R(p)` jest kwadratowa w `p`). Każda iteracja Newtona po `q` musi normalizować kwaterniony, inaczej poprawki „odbijają" i mogą dywergować; `projekcja_polozen` już to robi.
- `sym2` z domyślnym rzutowaniem wymaga `dt` w granicach stabilności półjawnego Eulera; za duży krok kończy się `RuntimeError` z komunikatem (Baumgarte przy grubym kroku przeżywa, ale kosztem naruszenia więzów).
- `SilaZewn` w `Pstrona` nadpisuje (nie sumuje) wartość dla danego członu i kierunku: dwie siły `ny` na ten sam człon nie zadziałają łącznie.

## Uwagi

- Zmiany w fizyce muszą zachowywać zgodność z opisem w PDF i przechodzić testy walidacyjne (okres wahadła, energia, więzy); nie refaktoryzuj silnika bez wyraźnej prośby.
- Komunikaty commitów po polsku, bez stopki `Co-Authored-By`.
- Po każdym istotnym kroku zaktualizuj status w `PLAN.md`.
