# Plan rozwoju projektu uw_dyn

Dokument dla przyszłych sesji: co jest zrobione, co dalej i jakie decyzje już zapadły.
Aktualizuj status po każdym istotnym kroku.

## Cel docelowy

Szybka biblioteka do dynamiki 3D układów wieloczłonowych, używalna z Pythona
jako zwykły pakiet (wzorzec: biblioteka `procenty`) oraz z wizualizacją
w przeglądarce.

Decyzja z 2026-07-05: NIE łączymy tego repo z projektem `~/repos/logistyka`;
logistyka to ujęcie makro, a uw_dyn to mikro-fizyka. Przykład transportowy
(`transport_teren.py`) został usunięty; gdyby kiedyś był potrzebny, jest
w historii gita (commit 68dc92d).

## Trzy przykłady kanoniczne (decyzja z 2026-07-05)

Bibliotekę rozbudowujemy, testujemy i weryfikujemy na trzech przykładach:

1. **Przysiad** (`przyklady/przysiad.py`): biomechanika, mięśnie
   sprężysto-tłumiące, równowaga statyczna, stateczność.
2. **Robot kroczący** (`przyklady/robot_kroczacy.py`, rozwijany w stronę
   mini-pieska czworonożnego): zmiany więzów w locie, uderzenia,
   sterowanie chodem, kontakt z podłożem.
3. **Kwadrokopter** (`przyklady/dron.py`): siły ciągu w układzie ciała,
   regulator dyskretny (PD po segmentach), ładunek podwieszony.

Każde nowe rozszerzenie biblioteki powinno być umotywowane potrzebą któregoś
z tych przykładów i pokryte testami; przykłady służą też jako testy
weryfikacyjne całości (fizyka + stabilność numeryczna + wizualizacja).

### Stan trzech filarów (2026-07-05)

- **Przysiad**: gotowy (mięśnie sprężysto-tłumiące, równowaga statyczna).
- **Robot kroczący → mini-piesek**: `robot_kroczacy.py` (compass gait, 2 człony)
  zostaje jako najprostszy przykład; `piesek.py` to czworonóg (tułów + 4 nogi
  dwuczłonowe, 9 członów, kontakt stóp `SilaKontaktu`, stawy `MomentWzgledny`,
  chód pełzający). Piesek idzie do przodu i nie przewraca się (~7.5 m / 5 s);
  chód jest otwarty (bez balansu), więc czuły na parametry: NIE stroić dalej
  bez potrzeby, obecne wartości w `piesek.py` są sprawdzone i spójne z danymi.
- **Kwadrokopter**: gotowy (`dron.py`, regulator kaskadowy PD 100 Hz, ładunek
  na linie, misja 3 punktów).

### Nowe prymitywy biblioteki wprowadzone dla filarów

- `SilaWPunkcie` (follower force) — ciąg wirnika drona.
- `SilaKontaktu` (kontakt penalty z tarciem) — stopy pieska, ładunek drona.
- `MomentWzgledny` (aktuator obrotowy w przegubie) — stawy pieska.
- `SilaWewnProst(tylko_rozciaganie=True)` — lina drona.
Wszystkie pokryte testami (40 testów). Dalsze kierunki: `docs/ULEPSZENIA.md`.

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

### Krok 3.5: restrukturyzacja i nowe metody numeryczne [ZROBIONE, 2026-07-05]

- [x] Podział monolitu na moduły: `algebra`, `czlony`, `wiezy`, `sily`,
      `uklad` (`uw_dyn.dynamika` został jako alias zgodności)
- [x] Stabilizacja więzów rzutowaniem (nowa domyślna w `sym2`):
      `projekcja_polozen` (Newton z kontrolą kroku; kwaterniony
      normalizowane przed iteracją, bo jakobiany są dokładne tylko dla
      kwaternionów jednostkowych) oraz `projekcja_predkosci` (rzutowanie
      w metryce macierzy mas, czyli uderzenie plastyczne). Baumgarte
      dostępne przez `stabilizacja='baumgarte'`.
- [x] Wykrywanie rozbiegania symulacji z czytelnym komunikatem
      (rzutowanie wymaga dt w granicach stabilności integratora;
      przykład łańcucha przeszedł z dt=0.1 na dt=0.02)
- [x] Metody energii: `energia_kinetyczna`, `energia_potencjalna`
      (grawitacja + sprężyny), `energia`; nowa metoda `SilaWewnProst.dlugosc`
- [x] Testy rzutowania i energii (razem 30 testów)

### Krok 4: wizualizacja i przykłady [ZROBIONE w wersji podstawowej, 2026-07-05]

- [x] `przyklady/przysiad.py` + `web/przysiad.html`: staw kolanowy podczas
      przysiadu (3 człony, mięśnie sprężysto-tłumiące, długości swobodne
      z równowagi statycznej, kontrola stateczności hesjanem energii)
- [x] `przyklady/robot_kroczacy.py` + `web/robot.html`: najprostszy robot
      kroczący (chód cyrklowy): noga podporowa przypinana przegubem,
      siłownik biodra ze stałymi momentami, tłumik między nogami; zmiana
      podpory = projekcja położeń + uderzenie plastyczne; 8 kroków, ~1.9 m.
      Parametry z przeszukiwania siatki (COM=0.8 przy biodrze kluczowe).
- [x] Strony web: Three.js z CDN (wymagają internetu); uruchamianie:
      `cd web && python3 -m http.server 8000` i otworzyć stronę
- [ ] Opcjonalnie: skrypt importu CSV do Blendera (`przyklady/lancuch.blend`)

### Krok 5: port rdzenia do Rust [WARUNKOWY: tylko gdy performance siada]

Uruchamiamy dopiero, gdy krok 3 nie wystarczy w realnym zastosowaniu.

- [ ] Crate `uw_dyn` (nalgebra), struktura kodu odwzorowująca wersję Python
- [ ] Bindingi PyO3 + maturin, to samo API co pakiet Python
- [ ] Testy zgodności: te same wejścia, porównanie trajektorii z wersją Python
- [ ] Benchmarki (kryterium sukcesu: co najmniej 100x szybciej niż Python)
- [ ] Kompilacja do WASM (wasm-bindgen), podpięcie pod wizualizację z kroku 4

## Stan techniczny (na 2026-07-05, po krokach 1-4)

- Pakiet: `uv sync` + `uv run pytest` (30 testów przechodzi w ~12 s)
- Moduły: `uw_dyn.algebra` (wektory, kwaterniony, macierze), `uw_dyn.czlony`,
  `uw_dyn.wiezy`, `uw_dyn.sily`, `uw_dyn.uklad`; `uw_dyn.dynamika` to alias
- Dwa integratory: `sym2` (półjawny Euler, stały krok; domyślnie
  stabilizacja rzutowaniem, opcjonalnie Baumgarte) i `sym` (solve_ivp
  RK45, adaptacyjny; dt określa tylko gęstość zapisu wyników)
- Nowe metody `Uklad`: `projekcja_polozen`, `projekcja_predkosci`
  (uderzenie plastyczne), `energia_kinetyczna/potencjalna/energia`
- Przykłady: `lancuch02.py` (CSV), `przysiad.py`, `robot_kroczacy.py`;
  wizualizacje w `web/*.html`
  (Three.js z CDN; `cd web && python3 -m http.server 8000`)
- Wydajność: łańcuch 4 członów ~1.65 ms/krok z Baumgarte,
  ~3.4 ms/krok z rzutowaniem (dokładne więzy)
- Znane ograniczenia: macierz układu gęsta (koszt rośnie sześciennie
  z N), więzy kierujące (`Kat`, `Odleglosc`) działają tylko w `newraph`
  (warunki początkowe), nie w dynamice; rzutowanie wymaga dt w granicach
  stabilności półjawnego Eulera (za duży krok kończy się czytelnym
  RuntimeError, wtedy zmniejszyć dt)
- Ważna subtelność: jakobiany więzów są dokładne tylko dla kwaternionów
  jednostkowych; wszelkie poprawki Newtona muszą normalizować kwaterniony
  (zrobione w `projekcja_polozen`)
- Teoria i oznaczenia: `docs/MSzalajski_mgr4.pdf` (praca magisterska, 2016)
