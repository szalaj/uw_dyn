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
  5 czlonow); sekwencja garda -> wyprowadzenie ciosu -> powrot.
  **Przebudowany na model 3D (2026-07-06):** barki to napedzane stawy kuliste
  (`Para_Sferyczna` + `MomentSferyczny`), lokcie to zawiasy (`Polaczenie_Obr`
  + `MomentWzgledny`), pas to przegub obrotowy wokol pionu (skret korpusu).
  Czlony ramion to prety wzdluz lokalnej z; orientacje ramion buduje
  `orientacja(dz, dx)` z kierunku, cele sterowania interpoluje slerp.
  Dzieki stawom kulistym garda jest NATURALNA (lokcie w dol, piesci przy
  twarzy), a cios idzie lukiem w przod i w poprzek na cien.
  Zweryfikowane: garda piesc przy policzku (0.17, -0.08, 1.41), szczyt ciosu
  (0.34, 0.14, 1.35) z uniesionym lokciem i obrotem korpusu (yaw 0.28 -> 0.60,
  CCW). `dt=0.0005`. Poprzednia wersja (plaszczyzna pozioma) w historii gita.

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
   Blender opcjonalnie do renderów.

## Mapa drogowa

### Krok 1: pakiet Python z testami [ZROBIONE, 2026-07-05]

- [x] Struktura pakietu `src/uw_dyn/` + `pyproject.toml` (uv, hatchling)
- [x] Poprawki zgodności z aktualnym NumPy/SciPy (skew, np.delete)
- [x] Usunięcie martwych importów (matplotlib, csv, sys, scipy.optimize)
- [x] Testy (pytest, 15 testów): algebra (skew, R, G, u2p, EA_to_EP),
      wahadło fizyczne (równowaga, okres drgań kontra wzór analityczny,
      zachowanie energii, spełnienie więzów), regresja łańcucha 4 członów
- [x] Przykład przeniesiony do `przyklady/`

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
- [ ] Opcjonalnie: skrypt importu CSV do Blendera

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

- [x] Pełna sylwetka z antropometrii — **ZROBIONE 2026-07-06**: moduł
      `antropometria.py` z tablicą Wintera (masy jako ułamek masy ciała,
      długości jako ułamek wzrostu, środki mas, promienie bezwładności) i
      funkcją `segmenty(masa, wzrost)` (12 segmentów, suma mas = masa ciała;
      dłoń scalona z przedramieniem). Tensory: pręt/walec wzdłuż lokalnej z.
      Builder `zbuduj_postac(masa, wzrost)` składa stojącą postać (tułów jako
      podstawa przypięta w miednicy, głowa, ramiona bark+łokieć, nogi
      biodro+kolano+kostka; barki/biodra/szyja kuliste, łokcie/kolana/kostki
      zawiasowe; aktuatory trzymają pozę neutralną, cele = zmierzona poza →
      punkt stały bez zgadywania znaków). Postać stoi stabilnie przy
      `dt=1e-4` (układ sztywny; większy krok wpada w rezonans). 8 testów
      (`test_antropometria.py`). Uwaga na przyszłość: rzutowanie położeń przy
      stawach kulistych bywa czułe (pozycja barku zależy od orientacji
      ramienia) — mały krok to obchodzi; docelowo warto rozważyć rzutowanie
      ważone macierzą mas.
- [x] Realistyczne stawy 3D — **ZROBIONE 2026-07-06**: dodano `MomentSferyczny`
      (napędzany staw kulisty 3 DOF = `Para_Sferyczna` + moment 3D w ramce
      globalnej: M = -k·φ - c·(ω_j-ω_i), φ z wektora obrotu błędu orientacji).
      Sprowadza człon do dowolnej zadanej orientacji względnej (`p_cel`,
      podmienialny w sterowaniu), z zerowym błędem, i trzyma pozę pod
      grawitacją. Bark i biodro modelujemy teraz jako staw kulisty, łokieć/
      kolano jako `Polaczenie_Obr` (1 DOF). 8 testów (`test_staw_kulisty.py`).
      To odblokowuje naturalne pozy (garda z łokciami w dół), których
      `bokser.py` nie miał (płaszczyzna pozioma). **ZROBIONE 2026-07-06:**
      `bokser.py` przebudowany na barki kuliste — garda jest teraz naturalna.
- [x] Limity zakresu ruchu w stawach — **ZROBIONE 2026-07-06**: `OgranicznikKata`
      (przegub 1 DOF: łokieć, kolano; jednostronna sprężyna-tłumik poza
      `[kat_min, kat_max]`, tłumienie tylko przy wchodzeniu głębiej) oraz
      `OgranicznikStozka` (staw kulisty: ogranicza odchylenie osi od kierunku
      neutralnego do `kat_max`). Model penalty (miękki: dopuszcza małe
      przejściowe wniknięcie, malejące ze sztywnością). 5 testów
      (`test_ograniczniki.py`). Można je dołożyć do stawów boksera/modelu.
- [x] Model „mięśni/napędów" stawów — **ZROBIONE 2026-07-06**: `MomentWzgledny`
      i `MomentSferyczny` mają opcjonalny `moment_max` (saturacja wartości/normy
      momentu). Słaby napęd nie dźwignie ciężaru do celu (siły ciosów fizyczne).
      UWAGA: w saturacji człon `k*błąd` przewyższa limit, więc tłumienie
      (też objęte limitem) traci wpływ i staw może oscylować; realny staw ma
      pasywne tłumienie tkanek, które dodaje się osobnym `MomentWzgledny`
      (k=0, c>0, bez limitu). 4 testy (`test_moment_max.py`).

**Pułapka symetrii kończyn (lekcja z boksera, 2026-07-06):** budując parzyste
kończyny (dwie ręce, dwie nogi) przez lustrzane odbicie względem płaszczyzny
strzałkowej (y → −y), trzeba ODWRÓCIĆ TAKŻE znak kątów zawiasów (łokieć,
kolano). Odbicie zmienia skrętność, więc to samo dodatnie zgięcie o tej samej
osi wygina lewy staw w drugą stronę (nieanatomicznie, np. przedramię do tyłu).
W `bokser.py` lewa ręka ma `flex = −2.60` (prawa `+2.60`). To samo dotyczy
całej pozy budowanej z kierunków: kierunki odbijamy w y, a kąty flexji negujemy.

### Etap B: sterowanie ruchem zawodnika

- [ ] Biblioteka gotowych ruchów jako trajektorie kątów stawów w czasie
      (jab, cross, hak, podbródkowy, kopnięcia: front/round/side), parametryzowane
      tempem i siłą. Wzorzec z `bokser.py`: sekwencja garda → cios → powrót,
      cele kątowe modulowane fazą (dyskretny sterownik co segment).
- [~] Postawienie postaci na nogach z kontaktem — **CZĘŚCIOWO 2026-07-06**:
      `zbuduj_postac(..., podparcie='stopy')` daje swobodną podstawę stojącą na
      stopach (bez pinu miednicy, kontakt `SilaKontaktu` w 4 punktach każdej
      stopy = wielobok podparcia, stawy trzymane sztywno). Kontakt podpiera
      ciężar (zweryfikowane statycznie). ALE: to wysoki odwrócony wahadeł na
      podatnym kontakcie, więc bez aktywnego balansu przewraca się bokiem
      (~0.25 m dryfu w 0.3 s). Wymaga też małego kroku (dt~5e-5, model sztywny)
      — sim jest wolny (0.05 s ≈ 35 s ściany), więc test kontaktu jest
      statyczny, nie symulacyjny.
- [x] Rzutowanie ważone macierzą mas — **ZROBIONE 2026-07-06**:
      `projekcja_polozen` liczy poprawkę w metryce mas (układ KKT, jak
      `projekcja_predkosci`). Dryf figury na stopach spadł z 0.247 do 0.115 m
      (mniej pasożytniczego obrotu), zestaw testów 118 → 47 s. Sufit `dt`
      figury wyznacza teraz sztywność aktuatorów/kontaktu, nie projekcja
      (dalej: integrator półniejawny/RATTLE, `docs/ULEPSZENIA.md`).
- [x] Regulator równowagi PID — **ZROBIONE 2026-07-06**: `przyklady/balans.py`.
      Dwie warstwy PID: (1) stawy trzymają pozę regulatorem PID (`ki` znosi
      sag pod grawitacją → wystarcza niższa sztywność → większy krok
      `dt=2e-4`); (2) regulator balansu PID sprzęga poziome położenie CoM
      z kostkami (przód-tył) i biodrami (bok). Kluczowe odkrycie: sama
      zmiana stawów na PID + rzutowanie ważone macierzą mas daje stabilne
      stanie (dryf CoM ~5 mm/0.4 s, wcześniej przewracał się o 0.1–0.25 m).
      Regulator balansu dokłada odporność na pchnięcie (dryf 0.23 → 0.13 m
      pod pchnięciem 0.5 m/s). Znak sprzężenia kostki ujemny (wyznaczony
      eksperymentalnie). Wizualizacja: `web/balans.html`. Ograniczenie:
      model kosztowny (dt=2e-4, ~340 s symulacji na 1 s ruchu); mocnego
      pchnięcia bez kroku w bok nie da się w pełni odrzucić.
- [~] Praca nóg / krok w bok — **CZĘŚCIOWO 2026-07-06**: `przyklady/krok.py`.
      Mechanika kroku złożona i działa: przeniesienie ciężaru (roll bioder),
      uniesienie nogi (kontakt penalty sam zwalnia stopę), wymach (odwiedzenie),
      postawienie; do tego lateralny PID (biodra) + PID przód-tył (kostki)
      utrzymują figurę na nogach przez fazę jednonożną (bez tego przewracała
      się). NIE osiągnięto czystego stabilnego kroku KIERUNKOWEGO: silne
      odwiedzenie stawia stopę wyraźnie w bok (~0.2 m), ale figura osiada
      (kolana się uginają pod dynamiką); umiarkowane trzyma balans, lecz
      pochył równoważący ściąga nogę do środka. Sprzężenie balans↔krok wymaga
      regulatora ustawienia stopy (capture point) i drobniejszego strojenia,
      które blokuje koszt symulacji (dt=2e-4). Wizualizacja: `web/krok.html`.
- [ ] Regulator ustawienia stopy (capture point) + szybszy integrator (RATTLE)
      dla stabilnego kroku i lokomocji.
- [ ] Praca nóg: krok, unik, zejście z linii — na bazie kontaktu i przenoszenia
      ciężaru.

### Etap C: interakcja i pomiary

- [x] Trafienie: kontakt pięści/stopy z celem — **ZROBIONE 2026-07-06**:
      nowa siła `SilaUderzenia` (kontakt penalty bryła-bryła: punkt vs kapsuła),
      worek bokserski w `bokser.py` to prawdziwe ciało (wahadło z masą na stawie
      kulistym), sierpowy i front kick trafiają w niego i wprawiają w wahnięcie
      (CoM worka do ~0.19 m). Wizualizacja: worek + lina w `web/bokser.html`.
- [x] Metryki uderzenia w worek — **ZROBIONE 2026-07-06**: `SilaUderzenia`
      liczy szczytową siłę kontaktu (`F_szczyt`) i impuls (`impuls(dt)`);
      w `bokser.py` raportowane osobno dla sierpowego i front kicku
      (np. kick ~875 N, impuls ~23 N·s), wyświetlane też w `web/bokser.html`.
- [ ] Metryki ciosu: prędkość i energia pięści w chwili kontaktu, impuls,
      moment obrotowy przenoszony na cel. Częściowo już liczone w `bokser.py`
      (prędkość pięści).
- [ ] Dwaj zawodnicy: dwa modele w jednym `Uklad`, wzajemne trafienia i bloki.

### Potrzebne rozszerzenia biblioteki (wynikają z etapów A–C)

- napędzany staw kulisty (3 DOF) — Etap A, najważniejsze;
- ograniczniki zakresu ruchu — Etap A;
- aktuator momentu z limitem — Etap A/B;
- kontakt bryła–bryła (nie tylko punkt–podłoże) do trafień — Etap C, ZROBIONE
  (`SilaUderzenia`: punkt vs kapsuła);
- regulator postawy/balansu — Etap B (najtrudniejsze numerycznie).

Powiązane: te rozszerzenia częściowo pokrywają się z listą w `docs/ULEPSZENIA.md`
(np. integrator dla sztywnych sprężyn, wykrywanie zdarzeń dla momentu kontaktu).

## Stan techniczny (na 2026-07-05, po krokach 1-4)

- Pakiet: `uv sync` + `uv run pytest` (30 testów przechodzi w ~12 s)
- Moduły: `uw_dyn.algebra` (wektory, kwaterniony, macierze), `uw_dyn.czlony`,
  `uw_dyn.wiezy`, `uw_dyn.sily`, `uw_dyn.uklad`; `uw_dyn.dynamika` to alias
- Trzy integratory: `sym2` (półjawny Euler, rząd 1; najodporniejszy na
  sztywne tłumienie — używać dla figury z kontaktem), `sym3` (RATTLE/Verlet,
  rząd 2, więzy maszynowo dokładne, energia ~500× lepsza; jawny wybucha przy
  sztywnym tłumieniu c/m·dt/2>1 lub sztywnych sprężynach ω·dt>2, tryb
  `polniejawne=True` znosi oba przez rzadkie jakobiany sił wcielane do KKT
  (tłumienie w kroku prędkości, sztywność w kroku położeń) — ~3 solve/krok,
  rząd 2 zachowany, figura biegnie na dt=5e-3 szybciej i dokładniej niż sym2)
  i `sym` (solve_ivp RK45, adaptacyjny; bez rzutowania i bez PID)
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
