# Plan rozwoju projektu uw_dyn

Dokument dla przyszłych sesji: co jest zrobione, co dalej i jakie decyzje już zapadły.
Aktualizuj status po każdym istotnym kroku.

## Cel docelowy

Szybka biblioteka do dynamiki 3D układów wieloczłonowych, używalna z Pythona
jako zwykły pakiet (wzorzec: biblioteka `procenty`) oraz z wizualizacją
w przeglądarce.

**Cel długoterminowy (decyzja z 2026-07-06): symulacja dynamiki sportów walki
człowieka, w szczególności kickboxingu, za pomocą uw_dyn.** Przykład `bokser.py`
(prawy sierpowy na cień) to pierwszy krok w tym kierunku. Docelowo: wiarygodny
model biomechaniczny zawodnika, ciosy i kopnięcia, praca nóg, garda, a w dalszej
perspektywie interakcja dwóch zawodników (trafienia, bloki). Dedykowana mapa
drogowa tego kierunku jest w sekcji „Kierunek: dynamika sportów walki" poniżej.

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
  zostaje jako najprostszy przykład chodu; `piesek.py` to czworonóg (tułów +
  4 nogi dwuczłonowe, 9 członów, kontakt stóp `SilaKontaktu`, stawy
  `MomentWzgledny`). Zaczynamy od PRZYSIADÓW, nie od chodu: piesek stoi na
  czterech nogach i symetrycznie ugina wszystkie stawy (kinematyka kolano =
  −2·biodro trzyma stopy dokładnie pod biodrami), tułów opada z 0.274 do
  0.169 m i wstaje. To stabilne i niewrażliwe na parametry. Chód (otwarty,
  bez balansu) był niestabilny i został odłożony; wersja z chodem jest
  w historii gita. Uwaga: głębsze zgięcie kolan wymaga `dt=0.0005` (przy
  `dt=0.001` całkowanie się rozbiega).
- **Kwadrokopter**: gotowy (`dron.py`, regulator kaskadowy PD 100 Hz, ładunek
  na linie, misja 3 punktów).

### Nowe prymitywy biblioteki wprowadzone dla filarów

- `SilaWPunkcie` (follower force) — ciąg wirnika drona.
- `SilaKontaktu` (kontakt penalty z tarciem) — stopy pieska, ładunek drona.
- `MomentWzgledny` (aktuator obrotowy w przegubie) — stawy pieska.
- `SilaWewnProst(tylko_rozciaganie=True)` — lina drona.
- `MomentWzgledny` (aktuator obrotowy w przegubie) — stawy pieska i boksera.
Wszystkie pokryte testami (40 testów). Dalsze kierunki: `docs/ULEPSZENIA.md`.

### Dodatkowy przyklad pokazowy

- **Bokser** (`przyklady/bokser.py` + `web/bokser.html`): kickboxing, walka
  z cieniem, prawy sierpowy. Gorna czesc ciala (tulow + 2 ramiona po 2 czlony,
  5 czlonow) sterowana aktuatorami `MomentWzgledny`; sekwencja garda ->
  wyprowadzenie ciosu -> powrot. Sierpowy to ruch poziomy, wiec WSZYSTKIE osie
  obrotu sa pionowe (z), a czlony ramion rozciagaja sie wzdluz lokalnej osi x
  (grawitacja nie daje momentu wzgledem pionu -> wiez trzyma ramie bez
  oklapniecia; jeden przegub obrotowy nie zlapie 3D barku, stad plaszczyzna
  pozioma). Zweryfikowane: garda piesc przy brodzie (x=0.16), szczyt ciosu
  x=0.40 z obrotem tulowia o 46 st. `dt=0.0005` (sztywne aktuatory).

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

## Kierunek: dynamika sportów walki (kickboxing) [PLAN, 2026-07-06]

Docelowe zastosowanie biblioteki: symulacja biomechaniki zawodnika sportów
walki. `bokser.py` (prawy sierpowy na cień) to punkt wyjścia. Poniżej mapa
drogowa; kolejność od najtańszego i najbardziej fundamentalnego.

### Etap A: solidny model człowieka (biomechanika)

- [ ] Pełna sylwetka jako łańcuch członów: miednica, tułów (1–2 segmenty),
      głowa, ramiona (bark+łokieć, docelowo +nadgarstek), nogi (biodro+kolano
      +kostka). Masy i tensory bezwładności z tablic antropometrycznych
      (np. de Leva / Winter), skalowane wzrostem i masą zawodnika.
- [ ] Realistyczne stawy 3D. Kluczowe ograniczenie obecnej biblioteki:
      `Polaczenie_Obr` to 1 stopień swobody. Bark i biodro są kuliste (3 DOF),
      a łokieć/kolano zawiasowe (1 DOF). ROZWIĄZANIE: albo złożyć staw kulisty
      z 2–3 przegubów obrotowych w serii z drobnym członem pośrednim, albo
      dodać do biblioteki napędzany staw kulisty (`Para_Sferyczna` + aktuator
      3-osiowy). To warunkuje naturalne pozy (garda z łokciami w dół), których
      `bokser.py` nie ma (model płaszczyzny poziomej — patrz jego opis).
- [ ] Limity zakresu ruchu w stawach (miękkie ograniczniki: jednostronna
      sprężyna-tłumik przy przekroczeniu kąta granicznego).
- [ ] Model „mięśni/napędów" stawów: aktuator momentu z ograniczeniem
      maksymalnego momentu (nie liniowa sprężyna bez limitu), żeby siły ciosów
      były fizycznie sensowne.

### Etap B: sterowanie ruchem zawodnika

- [ ] Biblioteka gotowych ruchów jako trajektorie kątów stawów w czasie
      (jab, cross, hak, podbródkowy, kopnięcia: front/round/side), parametryzowane
      tempem i siłą. Wzorzec z `bokser.py`: sekwencja garda → cios → powrót,
      cele kątowe modulowane fazą (dyskretny sterownik co segment).
- [ ] Utrzymanie równowagi na nogach z kontaktem stóp (`SilaKontaktu`) —
      to najtrudniejsze; open-loop się przewraca (jak chód pieska). Potrzebny
      regulator postawy (środek masy nad wielobokiem podparcia) i/lub praca nóg.
- [ ] Praca nóg: krok, unik, zejście z linii — na bazie kontaktu i przenoszenia
      ciężaru.

### Etap C: interakcja i pomiary

- [ ] Trafienie: kontakt pięści/stopy z celem (workiem, tarczą, przeciwnikiem)
      przez `SilaKontaktu`/`SilaWPunkcie`; pomiar siły i pędu uderzenia.
- [ ] Metryki ciosu: prędkość i energia pięści w chwili kontaktu, impuls,
      moment obrotowy przenoszony na cel. Częściowo już liczone w `bokser.py`
      (prędkość pięści).
- [ ] Dwaj zawodnicy: dwa modele w jednym `Uklad`, wzajemne trafienia i bloki.

### Potrzebne rozszerzenia biblioteki (wynikają z etapów A–C)

- napędzany staw kulisty (3 DOF) — Etap A, najważniejsze;
- ograniczniki zakresu ruchu — Etap A;
- aktuator momentu z limitem — Etap A/B;
- kontakt bryła–bryła (nie tylko punkt–podłoże) do trafień — Etap C;
- regulator postawy/balansu — Etap B (najtrudniejsze numerycznie).

Powiązane: te rozszerzenia częściowo pokrywają się z listą w `docs/ULEPSZENIA.md`
(np. integrator dla sztywnych sprężyn, wykrywanie zdarzeń dla momentu kontaktu).

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
