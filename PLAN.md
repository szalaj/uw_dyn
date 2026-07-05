# Plan rozwoju projektu uw_dyn

Dokument dla przyszłych sesji: co jest zrobione, co dalej i jakie decyzje już zapadły.
Aktualizuj status po każdym istotnym kroku.

## Cel docelowy

Szybka biblioteka do dynamiki 3D układów wieloczłonowych, używalna z Pythona
(w projektach takich jak `~/repos/logistyka`, wzorzec: biblioteka `procenty`)
oraz z wizualizacją w przeglądarce.

## Decyzje architektoniczne (ustalone z Marcinem)

1. **Póki co robimy wszystko w Pythonie** (decyzja z 2026-07-05): rozwój
   fizyki, API, testów i optymalizacji zostaje w wersji Python. To ona jest
   produktem i biblioteką używaną w innych projektach.
2. **Rust tylko warunkowo**: gdyby wydajność wersji Python przestała
   wystarczać w realnych zastosowaniach (mimo optymalizacji z kroku 3),
   przenosimy rdzeń do Rust (`nalgebra`), wystawiony:
   - dla Pythona przez PyO3 + maturin (pakiet instalowalny przez uv/pip),
   - dla przeglądarki przez WebAssembly (wasm-bindgen).
   Wtedy wersja Python zostaje implementacją wzorcową do testów zgodności.
3. **Wizualizacja: web + Three.js** jako podstawa (interaktywna, na żywo);
   Blender opcjonalnie do renderów (scena `przyklady/lancuch.blend`).

## Mapa drogowa

### Krok 1: pakiet Python z testami [ZROBIONE, 2026-07-05]

- [x] Struktura pakietu `src/uw_dyn/` + `pyproject.toml` (uv, hatchling)
- [x] Poprawki zgodności z aktualnym NumPy/SciPy (skew, np.delete)
- [x] Usunięcie martwych importów (matplotlib, csv, sys, scipy.optimize)
- [x] Testy (pytest, 15 testów): algebra (skew, R, G, u2p, EA_to_EP),
      wahadło fizyczne (równowaga, okres drgań kontra wzór analityczny,
      zachowanie energii, spełnienie więzów), regresja łańcucha 4 członów
- [x] Przykład przeniesiony do `przyklady/`, praca magisterska do `docs/`

### Krok 2: dopracowanie fizyki (wersja wzorcowa) [ZROBIONE, 2026-07-05]

- [x] Normalizacja kwaternionów po każdym kroku całkowania w `sym2`
      (normy jednostkowe do 1e-12, wcześniej dryf do ~0.996)
- [x] Naprawa aliasowania: `sym2` już nie modyfikuje `y0`
- [x] `sym` przepisane z ode/dopri5 na `solve_ivp` (RK45, adaptacyjny);
      test: okres wahadła z dokładnością 0.5%, dryf energii < 1e-4
- [x] Naprawione błędy fizyki: prędkości kątowe w `gammaK` i `SilaWewnProst`
      liczone jako `2*G(p)*p` (tożsamościowo zero!) zamiast `2*G(p)*dp`;
      brakujące `ri`/`rj` w `Para_Prostopadla_D.jakobianK` (i != 0);
      niezdefiniowane `Fqi` w `Kat.jakobianD` (i == 0)
- [x] Testy pozostałych więzów i sił: `Polaczenie_Przes` (spadek swobodny
      kontra wzór analityczny), `SilaZewn` (równowaga momentu), sprężyna
      i tłumik (równowaga statyczna, częstość drgań), `Kat` i `Odleglosc`
      (przez `newraph`), energia wahadła podwójnego (razem 25 testów)
- [x] Docstringi klas
- [ ] (przeniesione na później) API wysokopoziomowe, np. helper budowy łańcucha

### Krok 3: wydajność wersji Python [ZROBIONE, 2026-07-05]

- [x] Profilowanie (`cProfile`): wąskie gardło to potrójne liczenie
      `jakobianK` na krok i podwójne wywołania jakobianów w parach złożonych
- [x] Optymalizacje (trajektoria bitowo identyczna z wersją przed zmianą):
      memoizacja `jakobianK` per q, cache stałych macierzy masowych,
      prealokacja bloków zamiast `reduce(block_diag)`/`np.concatenate`,
      wybieranie kolumn zamiast `np.delete`, wektoryzacja `dq_jak`
- [x] Wynik: łańcuch 4 członów z 4.85 do 1.65 ms/krok (3x szybciej)
- [ ] (w razie potrzeby) numba lub głębsza wektoryzacja `gammaK`;
      profil jest już rozproszony po faktycznej fizyce

### Krok 4: wizualizacja [DO ZROBIENIA]

- [ ] Strona z Three.js: animacja z wyników symulacji (CSV/JSON),
      suwaki parametrów; symulacja liczona w Pythonie
- [ ] Opcjonalnie: skrypt importu CSV do Blendera (`przyklady/lancuch.blend`)

### Krok 5: port rdzenia do Rust [WARUNKOWY: tylko gdy performance siada]

Uruchamiamy dopiero, gdy krok 3 nie wystarczy w realnym zastosowaniu.

- [ ] Crate `uw_dyn` (nalgebra), struktura kodu odwzorowująca wersję Python
- [ ] Bindingi PyO3 + maturin, to samo API co pakiet Python
- [ ] Testy zgodności: te same wejścia, porównanie trajektorii z wersją Python
- [ ] Benchmarki (kryterium sukcesu: co najmniej 100x szybciej niż Python)
- [ ] Kompilacja do WASM (wasm-bindgen), podpięcie pod wizualizację z kroku 4

## Stan techniczny (na 2026-07-05, po krokach 2 i 3)

- Pakiet: `uv sync` + `uv run pytest` (25 testów przechodzi w ~8 s)
- Przykład: `uv run python przyklady/lancuch02.py` tworzy `lancuch.csv`
  (501 wierszy x 56 kolumn dla 4 członów; format: 3N pozycji, 4N parametrów
  Eulera, potem prędkości w tym samym porządku)
- Dwa integratory: `sym2` (półjawny Euler, stały krok, normalizacja
  kwaternionów, szybki) i `sym` (solve_ivp RK45, adaptacyjny, dokładny;
  dt określa tylko gęstość zapisu wyników)
- Wydajność: łańcuch 4 członów ~1.65 ms/krok (po optymalizacji 3x)
- Znane ograniczenia: stabilizacja Baumgarte'a z ręcznie dobieranymi
  alfa/beta, macierz układu gęsta (koszt rośnie sześciennie z N),
  więzy kierujące (`Kat`, `Odleglosc`) działają tylko w `newraph`
  (warunki początkowe), nie w dynamice
- Teoria i oznaczenia: `docs/MSzalajski_mgr4.pdf` (praca magisterska, 2016)
