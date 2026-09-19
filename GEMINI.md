# Instrukcja dla Gemini — atlas punktów akupunktury 3D

Pracujesz w repozytorium `acupuncture-3d` razem z Pawłem (lekarz, praktyka YAMA w Gdyni:
akupunktura i fizjoterapia). Pisz po polsku.

---

## 1. Czym to ma być

Narzędzie kliniczno-edukacyjne: **rysowanie map napięć** i **nauka akupunktury przez
pokazywanie adekwatnie umiejscowionych punktów**, z potencjałem weryfikacji diagnostycznej
i terapeutycznej. Dokładność położenia nie jest tu kwestią estetyki — punkt w złym miejscu
uczy złego miejsca.

## 2. Model, na którym stoi cały projekt

**Kanał to gałąź.** Ma znany przebieg i znane interpunkcje — punkty o określonej lokacji,
zawsze w tym samym miejscu względem siebie. **Wokół gałęzi rośnie tkanka**: mięśnie, ścięgna,
kości, powięzie. Tkanka jest tym, w czym gałąź jest osadzona, a nie tym, co ją definiuje.

Z tego wynika sposób pracy:

- **Punkty nie mają współrzędnych XYZ.** Każdy ma kotwicę proporcjonalną: segment ciała,
  odległość w cun, azymut albo odniesienie do krawędzi kości. Wymiana modelu ciała nie
  unieważnia niczego.
- **Punkt jest poprawny, gdy zgadzają się jego RELACJE**, nie gdy współrzędna trafia
  w liczbę: odstęp wzdłuż kanału, strona ścięgna, odległość od wyczuwalnej krawędzi kości,
  poziom kręgu.
- Sprawdzanie jednego punktu wobec własnego opisu to za mało. Odstęp między sąsiadami jest
  osobnym twierdzeniem i wiąże punkty, których opis nie zawiera żadnej liczby.

## 3. Zasady twarde — nie łam ich

1. **`data/point_meta.json` — NIE CZYTAJ i NIE CYTUJ.** To OCR książki Deadmana, materiał
   chroniony prawem autorskim. Wysłanie tych opisów do zewnętrznego modelu to wypuszczenie
   cudzego tekstu. Plik nigdy nie był i nie będzie commitowany.
   *Zamiast tego:* `data/cun_intervals.json` — kanoniczne odstępy w cun to fakty tradycji,
   nie czyjaś proza. Tego wolno używać i pokazywać.
2. **`assets/` nie trafia do repo** (siatki CC BY-SA, setki MB). Również `data/tension/`
   i `data/combos.json` zostają lokalnie.
3. **Przed publikacją czegokolwiek skanuj ZAWARTOŚĆ, nie nazwy plików:**
   `python tools/scan_leak.py`. Wyciek zdarzył się raz właśnie dlatego, że sprawdzano nazwy.
   Znane, świadomie odłożone trafienie: pole `categories` w `data/points.json` — **nie czyść
   go bez decyzji Pawła**.
4. **Nie commituj i nie pushuj bez wyraźnej prośby.** Przed wysyłką na GitHub pytaj wprost.
5. **Nie uruchamiaj `wsl --shutdown`** i nie wyłączaj Smart App Control.
6. Po każdej zmianie geometrii: `python tools/validate.py` musi dawać `WYNIK: OK`.

## 4. Jak uruchamiać

VTK jest blokowany na Windowsie przez Smart App Control, więc wszystko idzie przez WSL2:

```bash
wsl -d Ubuntu-24.04 -u root bash -c 'cd /mnt/c/Users/Paweł/repos/acupuncture-3d && /opt/acu3d-venv/bin/python tools/validate.py'
```

Narzędzia weryfikujące, od najszybszego:

| narzędzie | co sprawdza |
|---|---|
| `tools/validate.py` | czy 670 markerów leży na skórze i czy strony są symetryczne |
| `tools/verify_intervals.py` | odstępy w cun między sąsiednimi punktami kanału |
| `tools/verify_points.py` | pojedyncze punkty wobec twierdzeń z opisu |
| `tools/verify_edges.py` | czytelność krawędzi kości i wpasowanie kośćca |
| `tools/verify_fit.py` | ile tkanki leży nad kośćmi podskórnymi |
| `tools/verify_meridians.py` | linijka cun, przebiegi kanałów, przekazania (ok. 4 min) |
| `tools/measure_noise.py` | szum samych siatek — stąd biorą się progi |

## 5. Stan faktyczny — czytaj to zanim uwierzysz w procenty

Wskaźniki typu „54/54" opisują **tylko te punkty, których opis zawiera liczbę do wyłuskania**.
Realnie:

- **123 z 361 punktów (34 %)** ma choć jedno sprawdzone twierdzenie. HT 0/9, GB 4/44, DU 1/28.
- Liczba punktów i numeracja są za to w pełni poprawne: **361/361**, każdy kanał ciągły.
- Odstępy: **62/62** w tolerancji 0,5 cun, ale tabela pokrywa 70 z 341 par (20 %).
- **Punkty dystalne to najsłabsze miejsce:** 50 punktów na dłoni i stopie, wszystkie na
  kotwicach kostnych, tylko 6 z jakimkolwiek potwierdzeniem, **0 ze zmierzonym azymutem** —
  a tam kości Z-Anatomy wychodzą przez skórę BodyParts3D (30–100 % kości poza skórą, do 17 mm;
  na tułowiu 0 %).

Nierozstrzygnięte, nie próbuj tego „naprawić" bez nowych danych:

- **Brzeg przyśrodkowy piszczeli** — metoda „po kości" daje −100°, palpacyjna −25°.
  SP-6, SP-8, LR-5, LR-6 zostają nietknięte. Trzeba gęstszej siatki piszczeli z Blendera.
- **Miednica** jest węższa od skóry o ok. 7 mm równomiernie ze wszystkich stron — to różnica
  budowy dwóch dawców, nie przesunięcie. Żadne przesunięcie tego nie naprawi.
- **17 liczb w `data/cun_vectors.json`** ma `confirmed: false`. Cała linijka na nich stoi.
  Potwierdzić może tylko Paweł, z egzemplarzem normy WHO 2008 / GB/T 12346.
- **Konwencja −0,4 cun** (punkty Shu grzbietu pod czubkiem wyrostka) jest niezmierzona.
  Próba pomiaru z siatek kręgów dała wynik sprzeczny z anatomią — opis w `spinous_low()`.

## 6. Zadania, w kolejności wartości

1. **Relacje anatomiczne tam, gdzie nie ma cun.** Na twarzy, karku i tułowiu odstępy
   wyznaczają żebra, kręgi i linia sutkowa, nie odmierzone cun — dlatego 80 % par
   w `cun_intervals.json` ma `null`. Potrzebna jest analogiczna tabela **relacji**
   („na poziomie żebra N", „między ścięgnami X i Y", „w dołku przy kącie żuchwy"),
   którą da się sprawdzić maszynowo. To odblokowuje setki punktów naraz.
2. **Punkty dystalne.** 50 punktów, 6 sprawdzonych, kości źle wpasowane. Trzeba metody,
   która nie opiera się na kości wystającej ze skóry.
3. **Kolaterale, naczynia luo i przebiegi w głąb — przesiewowo.** To następny etap po
   punktach głównych. Nie punkt po punkcie: przesiew, który pokazuje, gdzie model
   w ogóle się nie trzyma.

## 7. Konwencje przy zmianach

- Zmieniasz kotwicę punktu → wpisz uzasadnienie w pole `az_source`, zaczynając od `pomiar:`,
  z liczbami przed i po. Narzędzie `tools/fit_azimuth.py` pomija takie punkty, żeby ich
  nie nadpisać.
- Nowe liczby z wiedzy ogólnej → `confirmed: false`, dopóki Paweł nie sprawdzi ich z normą.
- Nie masz pewności → wpisz `null` z notatką. **Dziura jest lepsza niż zmyślona liczba.**
- Nie dopasowuj kotwic do siatek mięśni „na oko". Kończyny mają 12–17 mm rozjazdu między
  skórą (BodyParts3D) a mięśniami i kośćmi (Z-Anatomy). Projekt raz już się o to potknął:
  `fit_azimuth.py` przyciągał punkty NA strukturę z opisu, zamiast obok niej.

## 8. Pułapka, która złapała mnie dwa razy — czytaj, zanim zgłosisz „atlas się myli"

Gdy nowy test masowo nie zgadza się z modelem, **najpierw podejrzewaj swoją tabelę
odniesienia i jednostkę pomiaru, dopiero potem atlas.**

Dwa autentyczne przykłady z budowy testu odstępów:

- Pierwsza wersja tabeli, pisana z pamięci, miała błędy 4–6 cun: LR-7 to punkt **przy
  kolanie**, nie 8 cun nad kostką; LU-6 → LU-7 to 5,5 cun, nie 1,5. Model miał rację.
- Drugi przebieg pokazał, że SP, KI i LR na goleni są „o 10–17 % za długie". To była
  **jednostka**: te punkty mają własną linijkę przyśrodkową 13 cun, a mierzone były
  standardową 16 na tej samej długości.

Dopiero po obu poprawkach wyszły dwie prawdziwe usterki (KI-10 bez własnej linijki,
LU-1 o 1,3 cun za nisko) — i te warto było zgłaszać.
