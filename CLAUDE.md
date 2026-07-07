# CLAUDE.md

Wskazówki dla Claude Code przy pracy z tym repozytorium.

## Czym jest ten projekt

Biblioteka Pythona do dynamiki przestrzennej układów wieloczłonowych (multibody dynamics). Orientacja członów opisana parametrami Eulera (kwaternionami), równania ruchu z mnożnikami Lagrange'a, stabilizacja więzów metodą Baumgarte'a.

**Najpierw przeczytaj `PLAN.md`**: tam jest mapa drogowa, podjęte decyzje architektoniczne i aktualny stan prac. Najważniejsze decyzje: rozwijamy wszystko w Pythonie (port do Rust tylko warunkowo, gdyby optymalizacja nie wystarczyła); bibliotekę rozbudowujemy i weryfikujemy na trzech przykładach kanonicznych: przysiad, robot kroczący (docelowo mini-piesek) i kwadrokopter.

## Struktura

- `src/uw_dyn/algebra.py`: wektory, kwaterniony (parametry Eulera), macierze R/G, skew.
- `src/uw_dyn/czlony.py`: `Czlon`.
- `src/uw_dyn/wiezy.py`: pary kinematyczne (`Para_Sferyczna`, `Polaczenie_Obr`, `Polaczenie_Cyl`, `Polaczenie_Przes`, `Para_Prostopadla`, `Para_Prostopadla_D`) i więzy kierujące (`Odleglosc`, `Kat`).
- `src/uw_dyn/sily.py`: `SilaWewnProst` (sprężyna/tłumik/siła stała, opcja `tylko_rozciaganie` dla lin), `SilaWPunkcie` (siła w układzie ciała zaczepiona w punkcie, np. ciąg wirnika), `SilaKontaktu` (kontakt penalty z podłożem z=0, tarcie regularyzowane), `SilaUderzenia` (kontakt penalty bryła-bryła kapsuła-kapsuła: najbliższe punkty dwóch odcinków wzdłuż lokalnych osi z + promienie; `polowa_wys_i=0` daje przypadek punkt-kapsuła, np. pięść/stopa w worek; kapsuła-kapsuła np. tułów-worek albo kończyna-przeciwnik; siła równa-przeciwna zaczepiona w punktach kontaktu → wahnięcie/odbicie; metryki uderzenia: `F_szczyt` [N] i `impuls(dt)` [N·s], zerowane `zeruj_metryki()`), `MomentWzgledny` (aktuator obrotowy 1 DOF: regulator PD, opcjonalnie PID przez `ki`+`calka_max`; błąd kąta zawijany do (−π,π]; opcjonalny `moment_max`), `MomentSferyczny` (napędzany staw kulisty 3 DOF: PD/PID na orientacji, łączyć z `Para_Sferyczna`), `OgranicznikKata` (miękki limit zakresu przegubu 1 DOF: łokieć, kolano), `OgranicznikStozka` (miękki limit stożka stawu kulistego: bark, biodro), `SilaZewn`.
- `src/uw_dyn/uklad.py`: `Uklad` (składanie równań ruchu, symulacja `sym`/`sym2`/`sym3`, `newraph`, rzutowania `projekcja_polozen`/`projekcja_predkosci`, metody energii, `zapiszWyniki`). Trzy integratory: `sym2` (półjawny Euler, rząd 1, najodporniejszy na sztywne tłumienie — domyślny dla figury z kontaktem), `sym3` (RATTLE/Verlet, rząd 2, więzy maszynowo dokładne, energia ~500× lepsza — do układów zachowawczych i umiarkowanie tłumionych; jawny wybucha przy sztywnym tłumieniu c/m·dt/2>1 lub sztywnych sprężynach ω·dt>2, a `sym3(..., polniejawne=True)` znosi oba: rzadkie jakobiany sił wcielane do KKT — tłumienie C=∂Q/∂v w kroku prędkości i sztywność K=∂Q/∂q w kroku położeń; koszt ~3 solve/krok, rząd 2 zachowany. Figura na stopach biegnie wtedy na dt=5e-3, 5× sufit sym2, szybciej i dokładniej), `sym` (RK45 adaptacyjny, bez rzutowania i PID).
- `src/uw_dyn/antropometria.py`: parametry segmentów ciała (tablica Wintera: masy, długości, środki mas, tensory bezwładności skalowane wzrostem i masą), `segmenty(masa, wzrost)`, oraz builder `zbuduj_postac(masa, wzrost, podparcie=...)` składający stojącą postać (12 członów, stawy kuliste w barkach/biodrach/szyi, zawiasy w łokciach/kolanach/kostkach). `podparcie='miednica'` (domyślne): tułów przypięty w miednicy, stabilne, `dt=1e-4`. `podparcie='stopy'`: swobodna podstawa na kontakcie stóp (4 punkty/stopa), quasi-sztywna; wymaga `dt~5e-5` i aktywnego balansu (bez niego przewraca się bokiem, patrz `PLAN.md` Etap B).
- `src/uw_dyn/dynamika.py`: alias zgodności wstecznej (re-eksport wszystkiego).
- `src/uw_dyn/__init__.py`: publiczne API pakietu (jawna lista `__all__`).
- `tests/`: pytest; `conftest.py` zawiera budowę wahadła testowego i obliczanie energii mechanicznej.
- `przyklady/`: `lancuch02.py` (CSV), trzy przykłady kanoniczne (`przysiad.py`, `robot_kroczacy.py` compass gait → `piesek.py` czworonóg robiący przysiady `dt=0.0005`, `dron.py` kwadrokopter) oraz `bokser.py` (kickboxing: sierpowy + front kick, nogi ze stawem kolanowym; worek to prawdziwe ciało/ciężkie wahadło 22 kg z tłumionym przegubem (`MomentSferyczny` k=0 jako tłumik, worek wraca do pionu), kontaktem `SilaUderzenia` i metrykami uderzenia), `pompka.py` (push-up: deska ciała na przegubie w palcach stóp, ramiona-zawiasy, dłonie na kontakcie, łokcie PID; `sym3` półniejawny) i `balans.py` (Etap B: pełna sylwetka na stopach z regulatorem balansu PID, `dt=2e-4`, kosztowny) i `czworaka.py` (człowiek na czworaka pełznący do przodu: proporcje Wintera, stawy kuliste w barkach/biodrach + zawiasy w łokciach/kolanach, kontakt dłoni i kolan, chód pełzający z przenoszeniem ciężaru; sterowanie postawą względem tułowia ODNIESIENIA — cele liczone dla idealnej wysokości/poziomu, więc kończyny odpychają zapadający się tułów; `moment_max` chroni przed kopnięciami, `dt=3e-4`).
- `web/`: wizualizacje Three.js (`przysiad.html`, `robot.html`, `piesek.html`, `dron.html`, `bokser.html`, `balans.html`, `krok.html`, `pompka.html`, `czworaka.html`); pliki `dane_*.js` generują skrypty z `przyklady/`.
- `docs/ULEPSZENIA.md`: research kierunków rozwoju biblioteki z priorytetami.

## Nowe siły (interfejs `sila(q, dq, N) -> (Qr_i, Qp_i, Qr_j, Qp_j)`)

Każda siła wewnętrzna implementuje `sila(...)` oraz `energia_potencjalna(q, N)`. `SilaWPunkcie`, `SilaKontaktu` i `MomentWzgledny` ustawiają `i=0` i działają na człon `j` (jak siły „do podstawy"). Ich wektory (ciąg, cel kąta, `p_cel` orientacji) można podmieniać między segmentami symulacji, co daje dyskretny sterownik (regulator drona, chód pieska) bez zmian w silniku. `MomentSferyczny` nie wiąże położenia — to napędzany staw kulisty składany z `Para_Sferyczna` (więz kuli w panewce) + moment 3D; cel `p_cel` to orientacja względna (parametry Eulera, z układu i do j). To fundament Etapu A (biomechaniczny model człowieka: bark i biodro jako stawy kuliste).

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
- Człon całkujący aktuatorów (`ki`, PID) jest aktualizowany tylko w `sym2` (hook `aktualizuj_calke` raz na krok), nie w `sym` (RK45). Za duże `ki` destabilizuje; stroić ostrożnie, używać `calka_max` (anti-windup).

## Uwagi

- Zmiany w fizyce muszą zachowywać zgodność z opisem w PDF i przechodzić testy walidacyjne (okres wahadła, energia, więzy); nie refaktoryzuj silnika bez wyraźnej prośby.
- Komunikaty commitów po polsku, bez stopki `Co-Authored-By`.
- Po każdym istotnym kroku zaktualizuj status w `PLAN.md`.
