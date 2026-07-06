# Research: co poprawić w module uw_dyn

Stan na 2026-07-05, po restrukturyzacji pakietu i wprowadzeniu stabilizacji
rzutowaniem. Analiza pod kątem trzech przykładów kanonicznych (przysiad,
robot kroczący/mini-piesek, kwadrokopter), które wyznaczają priorytety.

## 0. Audyt matematyczny (2026-07-06): sformułowanie ZWERYFIKOWANE

Numeryczna weryfikacja poprawności (utrwalona w `tests/test_matematyka.py`):

1. **Jakobiany więzów** = pochodne kierunkowe Φ wzdłuż wszystkich kierunków
   stycznych (dp ⊥ p), dla wszystkich 9 typów par (gałęzie i=0 oraz i≠0),
   zgodność ~1e-9. Wzdłuż kierunku radialnego kwaternionu są (udokumentowanie)
   niedokładne, stąd normalizacja w iteracjach Newtona.
2. **Człony gamma** (prawa strona przyspieszeń): γ == −(dJ/dt)·v dla losowych
   stanów, zgodność ~1e-9 dla wszystkich par.
3. **Siły potencjalne**: Q == −∂V/∂q wzdłuż kierunków stycznych (sprężyna,
   staw kuliste PD, ograniczniki, kontakt); zawiasy egzaktne na rozmaitości
   przegubu (poza nią wzory kątowe z definicji się rozjeżdżają — nieistotne,
   bo więz trzyma układ na rozmaitości).
4. **Bryła swobodna** (asymetryczny tensor, oś pośrednia, test Dzhanibekova):
   kręt w układzie świata zachowany, a błąd maleje DOKŁADNIE liniowo z dt
   (stosunki 2.00 przy połowieniu kroku) → jedyne źródło błędu to rząd 1
   półjawnego Eulera, nie sformułowanie równań (człony 4GᵀJG, 8ĠᵀJĠp są
   egzaktne).

Wniosek: matematyka biblioteki jest poprawna; głównym kierunkiem poprawy
pozostaje RZĄD i koszt integracji (punkt 3 poniżej: RATTLE / grupa Liego),
nie poprawność.

## 1. Siły: największa luka funkcjonalna

Obecnie biblioteka zna tylko trzy rodzaje sił: grawitację, stałą siłę lub
moment w osiach globalnych (`SilaZewn`) oraz element sprężysto-tłumiący
między punktami (`SilaWewnProst`). Wszystkie przykłady obchodzą ten brak,
tnąc symulację na segmenty i podmieniając parametry sił między segmentami
(dyskretny regulator). To uczciwa technika, ale są luki, których nie da się
nią załatać:

- **Siła w punkcie ciała, w układzie ciała** (follower force): ciąg wirnika
  drona działa wzdłuż osi kadłuba i jest zaczepiony poza środkiem masy.
  `SilaZewn` działa w osiach globalnych i w środku masy. To rozszerzenie
  jest tanie: uogólnione siły to `Qr = R f'` oraz `Qp = 2 G^T (s' x f')`.
  [wdrażane teraz jako `SilaWPunkcie`]
- **Kontakt z podłożem** (jednostronna sprężyna-tłumik pod punktem ciała,
  model penalty): bez tego czworonóg wymagałby żonglowania więzami jak
  chodziarz cyrklowy, co przy 4 nogach przestaje być praktyczne.
  [wdrażane teraz jako `SilaKontaktu`]
- **Siły zależne od prędkości** (opór aerodynamiczny ~ v|v|): potrzebne,
  gdyby dron miał latać szybko albo pojawił się przykład lotniczy.
- Docelowo: wspólny interfejs sił `sila(q, dq, t)` i możliwość podania
  własnej funkcji (callback), co domknie wszystkie powyższe przypadki.
- **Człon całkujący w aktuatorach (PID)** — **ZROBIONE 2026-07-06**:
  `MomentWzgledny` i `MomentSferyczny` maja opcjonalny `ki` (regulator PID
  zamiast PD) z anti-windup (`calka_max`). Całkę aktualizuje `sym2` raz na
  krok (hook `aktualizuj_calke`), więc działa tylko w `sym2` (nie w `sym`
  RK45). Znosi błąd ustalony pod grawitacją (wahadło: sag PD 0.25 rad →
  PID 0.02 rad). Domyślnie `ki=0` (pełna zgodność wsteczna, PD). Uwaga:
  za duże `ki` destabilizuje (klasyczny nadmiar całki) — stroić ostrożnie.
  Dla szybko zmiennych celów (ciosy) integral mniej pomaga (PD dominuje);
  główny zysk to trzymanie pozy (postawa, stanie).

## 2. Sterowanie: callback zamiast ręcznego cięcia na segmenty

Wzorzec z przykładów (pętla po segmentach + mutacja parametrów sił) powinien
trafić do biblioteki jako opcjonalny callback w `sym2`, na przykład
`sterowanie(t, y, ukl)` wywoływany co `co_ile` kroków. Zmniejszy to
powielanie kodu i pozwoli na zdarzenia (patrz punkt 4).

## 3. Integracja numeryczna

- Półjawny Euler + rzutowanie to dobry, przewidywalny domyślny schemat;
  jego realny limit to **sztywne sprężyny** (przysiad: k=1e5 wymusza
  dt=1e-3). Kierunki, od najtańszego:
  1. linearyzacja sprężyn w kroku (półniejawny w siłach sprężystych),
  2. RATTLE/SHAKE — **ZROBIONE 2026-07-06 jako `sym3`**: Verlet
     prędkościowy z więzami wewnątrz kroku (SHAKE = `projekcja_polozen`
     w metryce mas na pozycjach, RATTLE = `projekcja_predkosci`), plus
     iteracja punktu stałego na prędkości końcowej (konieczna dla rzędu 2
     przy członie żyroskopowym). Zweryfikowane: rząd 2 (stosunek błędu
     dokładnie 4.00 przy dt/2), więzy ~1e-13, energia wahadła/łańcucha
     500–740× dokładniejsza niż sym2 przy tym samym koszcie kroku; przy
     dt 10× większym nadal dokładniejszy i ~9× szybszy. Dron przechodzi
     pełną misję na sym3. OGRANICZENIE: przy bardzo sztywnym tłumieniu
     na lekkich członach (figura: c/J·dt/2 > 1) iteracja punktu stałego
     nie zbiega i sym3 wybucha tam, gdzie sym2 (półjawny w prędkości)
     przeżywa — dla pełnej sylwetki z kontaktem zostaje sym2. Docelowe
     rozwiązanie: tłumienie półniejawne w kroku prędkości.
  3. integratory na grupie Liego (SO(3)): kwaternion aktualizowany
     mnożeniem przez kwaternion przyrostowy, norma zachowana z definicji,
     znika też problem jakobianów przy nieunormowanych kwaternionach.
- `sym` (RK45 z solve_ivp) jest dokładny, ale nie stabilizuje więzów
  (tylko Baumgarte przez alfa/beta) i nie przyjmuje callbacków; po
  wprowadzeniu wspólnego interfejsu sił warto ujednolicić.
- **Rzutowanie ważone macierzą mas** (`projekcja_polozen`) — **ZROBIONE
  2026-07-06**: zamiast poprawki min-normowej `projekcja_polozen` rozwiazuje
  uklad KKT `[[M, J^T],[J,0]][dq;l]=[0;F]` (minimalizuje `dq^T M dq` przy
  `J dq = F`), spójnie z `projekcja_predkosci`. Efekt: mniej pasożytniczego
  obrotu lekkich członów (dryf figury na stopach z 0.247 do 0.115 m w 0.3 s),
  lepsze uwarunkowanie i szybszy zestaw testów (118 → 47 s). Uwaga: sufit
  kroku `dt` pełnej sylwetki wyznacza teraz SZTYWNOŚĆ aktuatorów i kontaktu
  (stabilność jawnego Eulera), nie projekcja — to osobny kierunek (integrator
  półniejawny/RATTLE, punkt 3 powyżej).

## 4. Zdarzenia (event detection)

Chodziarz wykrywa lądowanie skanem trajektorii po fakcie, z krokiem dt.
Bisekcja zdarzenia (funkcja zdarzenia + zawężanie) dałaby dokładny moment
przełączenia i mniejsze artefakty uderzeń. Naturalne API: `sym2(...,
zdarzenie=funkcja)` przerywające całkowanie w miejscu zera funkcji.

## 5. Wydajność

Profil po optymalizacjach: czas rozkłada się na `gammaK`/`jakobianK`
(faktyczna fizyka) i drobne operacje NumPy. Rezerwy, od najtańszej:

- cache macierzy `R(p)`, `G(p)` per człon w obrębie kroku (liczone
  wielokrotnie przez różne więzy tego samego ciała),
- numba/cython na `gammaK` par kinematycznych,
- macierz układu jest blokowo-rzadka: scipy.sparse albo eliminacja Schura
  (najpierw zmierzyć, dziś N jest małe i gęsta LU wygrywa),
- sformułowania O(N) (Featherstone, articulated body) dopiero przy długich
  łańcuchach; duża przebudowa, nieopłacalna przy N < ~20.

Mini-piesek (9 ciał, ~53 więzy) będzie pierwszym realnym testem skalowania.

## 6. Poprawność i weryfikacja

- Jakobiany więzów są dokładne tylko dla kwaternionów jednostkowych
  (macierz R jest kwadratowa w p). Obejście (normalizacja w iteracjach
  Newtona) działa, ale integrator Liego (punkt 3.3) usunąłby problem
  u źródła.
- Warto dodać zadania benchmarkowe IFToMM (np. mechanizm Bricarda,
  wahadło N-członowe z literatury) jako testy referencyjne oraz test
  rzędu zbieżności integratora (błąd vs dt).
- Więzy kierujące (`Kat`, `Odleglosc`) działają tylko w `newraph`;
  włączenie ich do dynamiki dałoby napędy kinematyczne (zadany przebieg
  kąta w czasie), przydatne do prowadzenia chodu i manipulatorów.

## 7. API i ergonomia

- Helper budowy łańcuchów/mechanizmów (powtarzalny kod w przykładach:
  budowa q0 z kątów, przeguby seryjne).
- Sprężyna i tłumik obrotowy w przegubie (moment ~ kąt względny i prędkość
  względna): dziś emulowane linkami liniowymi (mięśnie, tłumik biodra),
  co działa, ale ma nieliniową geometrię ramion.
- Walidacja wejścia (wymiary wektorów, spójność indeksów członów)
  z czytelnymi błędami.
- Zapis/odczyt modelu (np. YAML) w dalszej perspektywie.

## Priorytety (pod trzy przykłady kanoniczne)

| Priorytet | Co | Dla kogo |
|---|---|---|
| P1 (teraz) | `SilaWPunkcie` (ciąg w układzie ciała) | kwadrokopter |
| P1 (teraz) | `SilaKontaktu` (penalty, tarcie wiskotyczne) | mini-piesek |
| P2 | callback sterowania w `sym2` | wszystkie |
| P2 | sprężyna/tłumik obrotowy przegubu | piesek, przysiad |
| P2 | wykrywanie zdarzeń | piesek, robot |
| P3 | integrator dla sztywnych sprężyn (RATTLE / Lie) | przysiad |
| P3 | cache R/G per krok, numba na gammaK | wszystkie (N rośnie) |
| P4 | sparse / Featherstone O(N) | dopiero przy dużych N |
| P4 | benchmarki IFToMM, testy rzędu zbieżności | wiarygodność |
