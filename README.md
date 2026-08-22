# Atlas punktów akupunktury 3D

Przeglądarka 3D prawdziwego modelu anatomicznego z **361 punktami** klasycznych 14 kanałów.
Warstwy ciała (skóra / mięśnie / kości) przełączasz w locie, kanały włączasz pojedynczo,
a kliknięcie punktu pokazuje jego opis i **mięsień, który leży pod igłą**.

Mięśnie są osobnymi obiektami: klikasz konkretny mięsień, dostajesz jego nazwę, długość,
objętość i listę punktów akupunkturowych leżących nad nim — a suwakami zmieniasz jego
napięcie, tonus spoczynkowy, przerost i sztywność.

Python + PyVista (VTK).

![Wszystkie kanały, widok z przodu](docs/przod-wszystkie-kanaly.png)

## Skąd biorą się położenia punktów

Punkty **nie mają wpisanych współrzędnych XYZ**. Każdy opisany jest kotwicą proporcjonalną —
dokładnie tak, jak opisuje go podręcznik:

```json
{ "code": "ST-36", "anchor": { "seg": "shank", "cun": 3, "from": "p0", "az": 22 } }
```

„Na podudziu, 3 cun poniżej ST-35, szerokość palca bocznie od grzebienia piszczeli”.
Dopiero rig zamienia to na metry. Linijka *cun* nie jest zgadnięta — liczy się ją
z **prawdziwych kości** modelu (`tools/calibrate_rig.py`): wyrostków kolczystych kręgów,
wcięcia szyjnego mostka, wyrostka mieczykowatego, górnego brzegu spojenia łonowego,
krętarza większego, szpary stawu kolanowego i kostek.

Dla tego modelu wychodzi 2,48 cm na 1 cun w pionie tułowia i 2,07 cm w poziomie —
i stąd BL-23 ląduje 1,5 cun od linii pośrodkowej, a BL-52 3 cun, obydwa na wysokości L2.

Konsekwencja praktyczna: **wymiana modelu ciała nie unieważnia żadnej z 361 kotwic.**

## Warstwy

![Skóra, mięśnie, kości](docs/warstwy.png)

| warstwa | źródło | liczby |
|---|---|---|
| skóra | BodyParts3D (`FMA7163`) | 90 tys. trójkątów |
| mięśnie | Z-Anatomy, 477 osobnych obiektów | 1,08 mln trójkątów |
| kości | Z-Anatomy, 1244 obiekty scalone | 220 tys. trójkątów |

Oba zbiory są w innych układach odniesienia. Przesunięcia nie zgadujemy — liczymy je,
porównując położenie **tej samej kości** (krzyżowej i żuchwy) w obu zbiorach.
Skala wychodzi 1.0, co potwierdza, że to fizycznie ten sam model.

## Mięśnie: wybór i modulacja

![Wybrany mięsień z parametrami](docs/miesien.png)

Kliknięcie w mięsień zaznacza go na żółto i pokazuje panel. Mięsień kurczy się wzdłuż
swojej osi głównej, a przyczepy zostają na kościach:

```
s    = s_max · (napięcie · sztywność + tonus)     skrócenie, s_max = 0,16
w(t) = sin(π·t)                                   0 na przyczepach, 1 w brzuścu
u'   = u · (1 − s·w)                              wzdłuż osi
r'   = r / √(1 − s·w)                             w poprzek — objętość bez zmian
```

Pogrubienie brzuśca nie jest osobnym suwakiem: wynika z zachowania objętości.
Przerost i zanik to druga, niezależna gałka, bo to inne zjawisko.

Powiązanie działa w obie strony:

- wybierz mięsień → dostajesz punkty leżące nad nim wraz z odległością
  (mięsień brzuchaty łydki wskazuje BL-57 Chengshan w 4 mm, naramienny — LI-15 w 6 mm),
- kliknij punkt → dostajesz mięsień pod nim i głębokość, na jakiej się zaczyna.

## Instalacja

```bash
python -m venv .venv
.venv\Scripts\pip install -r requirements.txt
```

Siatek anatomicznych **nie ma w repozytorium** (setki MB, licencja CC BY-SA).
Pobiera je i przetwarza jeden skrypt — potrzebny jest zainstalowany
[Blender](https://www.blender.org/) (Z-Anatomy trzyma całą anatomię w pliku `.blend`):

```bash
.venv\Scripts\python tools\fetch_assets.py
```

Skrypt pobiera Z-Anatomy (~87 MB) i BodyParts3D (~140 MB), eksportuje siatki, obraca je
do naszego układu, decymuje i kalibruje rig. Trwa to kilkanaście minut, raz.

Bez pobranych zasobów aplikacja i tak działa — przechodzi na wbudowany model proceduralny
(`python -m acu3d --procedural`), z tymi samymi punktami na uproszczonym manekinie.

## Uruchomienie

```bash
.venv\Scripts\python -m acu3d
```

```bash
.venv\Scripts\python -m acu3d --channel ST --labels --view front
```

```bash
.venv\Scripts\python -m acu3d --muscle gastrocnemius --activation 0.9 --layer muscle
```

## Sterowanie

| | |
|---|---|
| lewy przycisk myszy | obrót; kółko — zoom; prawy — przesunięcie |
| klik w punkt | opis w oknie + pełny opis na konsoli + **automatyczny wybór mięśnia pod punktem** |
| klik w mięsień | zaznaczenie i panel parametrów (gdy warstwa mięśni jest włączona) |
| **`f`** | **fokus**: reszta ciała staje się prawie przezroczysta, a wybrany punkt i mięsień trafiają w środek obrotu; `f` ponownie — powrót do poprzedniego widoku |
| `1` `2` `3` `4` | widok przód / tył / bok / góra, `0` — widok domyślny |
| `a` / `n` / `i` | wszystkie kanały / żaden / odwróć zaznaczenie |
| `c` `v` | solo kolejnego kanału (z etykietami), `x` — wszystkie z powrotem |
| `m` `l` | linie kanałów / etykiety kodów |
| `k` `j` `b` | skóra / mięśnie / kości |
| `r` | wyzeruj parametry wybranego mięśnia |

Po lewej: przełączniki warstw (dół), lista kanałów z przyciskami **wszystkie / żaden / odwróć**.
Po prawej: suwaki parametrów mięśnia, przyciski **fokus** i **zeruj mięsień**, krycie skóry.

### Wydajność

Ruch suwakiem przelicza jedną klatkę w **ok. 5 ms**. Osiąga się to dwoma sztuczkami:
wybrany mięsień jest rysowany jako osobny aktor (a w wielkiej siatce wszystkich mięśni
jego komórki stają się przezroczyste), więc na kartę idzie kilka tysięcy wierzchołków
zamiast pół miliona; mapowanie punkt↔mięsień liczy się raz na mięsień i zostaje w pamięci,
bo nie zależy od parametrów.

![Kanał żołądka](docs/kanal-st.png)

## Struktura

```
acu3d/
  anatomy.py    cialo z siatek + linijka cun z kosci; powierzchnia szukana promieniem
  rig.py        model proceduralny (zapasowy) - proporcje z tabeli
  anchors.py    kotwica -> wspolrzedne 3D
  muscles.py    477 miesni w jednej siatce, deformacja aktywacyjna
  model.py      spina cialo, punkty i meridiany
  viewer.py     okno, warstwy, wybor punktu i miesnia
data/
  points.json      361 kotwic + polskie notatki o lokalizacji
  muscles.json     katalog miesni: os dzialania, objetosc, strona
  rig_fitted.json  punkty kostne i linijka cun zmierzone na modelu
  meridians.json   przebiegi kanalow
tools/
  fetch_assets.py      pobranie i przygotowanie calej anatomii (jeden skrypt)
  blender_export.py    eksport siatek z pliku .blend
  prepare_meshes.py    obrot, decymacja, os glowna miesnia
  prepare_skin.py      skora BodyParts3D dopasowana do anatomii
  calibrate_rig.py     linijka cun z prawdziwych kosci
  build_points.py      tabele kotwic -> data/points.json
  validate.py          kontrola spojnosci atlasu
```

## Dokładność

`tools/validate.py` sprawdza, czy wszystkie punkty się rozwiązują, leżą na skórze,
czy strony są symetryczne i czy odstępy między sąsiadami są prawdopodobne.

Każdy punkt ma znacznik `accuracy`:

- **anchored** (272) — kotwica odwzorowuje opis podręcznikowy wprost,
- **approx** (89) — model nie ma tej struktury (palce, małżowina, oczodół),
  więc położenie jest poglądowe i wymaga sprawdzenia wzrokiem.

To nie jest wyrób medyczny ani pomoc do samodzielnego nakłuwania. Narzędzie do nauki
topografii — punkt na ekranie zawsze weryfikuj z podręcznikiem i z żywym pacjentem.

## Licencje

Kod — MIT. Dane wyprowadzone z modeli anatomicznych — CC BY-SA.
Opisów punktów z podręcznika Deadmana **nie ma w repozytorium**.
Szczegóły i wymagana atrybucja: [ATTRIBUTION.md](ATTRIBUTION.md).
