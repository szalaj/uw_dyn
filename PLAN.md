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

### Krok 2: dopracowanie fizyki (wersja wzorcowa) [DO ZROBIENIA]

- [ ] Normalizacja kwaternionów po każdym kroku całkowania w `sym2`
      (dziś norma dryfuje do ~0.996 przy dłuższych symulacjach)
- [ ] Naprawa aliasowania: `sym2` modyfikuje `y0` w miejscu przez widoki
      `q`/`dq` (testy obchodzą to przez `y0.copy()`)
- [ ] Użycie/naprawa ścieżki `sym` (dopri5) lub przejście na `solve_ivp`;
      porównanie dokładności z `sym2`
- [ ] Testy dla pozostałych więzów: `Polaczenie_Cyl`, `Polaczenie_Przes`,
      `Odleglosc`, `Kat`, `SilaZewn`, sprężyna/tłumik w `SilaWewnProst`
- [ ] Sensowne API wysokopoziomowe (np. budowa łańcucha helperem) i docstringi

### Krok 3: wydajność wersji Python [DO ZROBIENIA, gdy zajdzie potrzeba]

- [ ] Profilowanie (`cProfile`); wąskim gardłem jest narzut Pythona przy
      budowie małych macierzy, nie algebra (macierz układu dla N=4 to ~52x52)
- [ ] Wektoryzacja składania macierzy (prealokacja zamiast `np.concatenate`)
- [ ] Ewentualnie numba dla gorących ścieżek

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

## Stan techniczny (na 2026-07-05)

- Pakiet: `uv sync` + `uv run pytest` (15 testów przechodzi w ~6 s)
- Przykład: `uv run python przyklady/lancuch02.py` tworzy `lancuch.csv`
  (501 wierszy x 56 kolumn dla 4 członów; format: 3N pozycji, 4N parametrów
  Eulera, potem prędkości w tym samym porządku)
- Znane słabości fizyki: półjawny Euler ze stałym krokiem w `sym2`,
  dryf norm kwaternionów, stabilizacja Baumgarte'a z ręcznie dobieranymi
  alfa/beta, macierz układu gęsta (koszt rośnie sześciennie z N)
- Teoria i oznaczenia: `docs/MSzalajski_mgr4.pdf` (praca magisterska, 2016)
