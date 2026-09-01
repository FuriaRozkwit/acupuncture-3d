# Źródła i licencje

Repozytorium łączy trzy rzeczy o **różnym statusie prawnym**. Poniżej dokładnie, co skąd
pochodzi i co z tego wynika.

## 1. Kod — MIT

Cały kod w `acu3d/` i `tools/` jest mój i objęty [licencją MIT](LICENSE).

## 2. Geometria anatomiczna — CC BY-SA, **nie ma jej w repozytorium**

Siatki 3D nie są tu trzymane. Pobiera i przetwarza je `tools/fetch_assets.py`
do katalogu `assets/`, który jest wyłączony z repozytorium (`.gitignore`).

| co | źródło | licencja |
|---|---|---|
| mięśnie, kości, kręgi, przyczepy | [Z-Anatomy](https://github.com/Z-Anatomy/Models-of-human-anatomy) — Gauthier Kervyn, Marcin Zieliński | CC BY-SA 4.0 |
| skóra (powłoka ciała, `FMA7163`) | [BodyParts3D](https://dbarchive.biosciencedbc.jp/en/bodyparts3d/) — Database Center for Life Science | CC BY-SA 2.1 JP |

Wymagana atrybucja przy każdym rozpowszechnianiu tych siatek lub ich przeróbek:

> BodyParts3D, © The Database Center for Life Science, licensed under
> CC Attribution-Share Alike 2.1 Japan.
>
> Z-Anatomy, © Gauthier Kervyn, licensed under CC Attribution-Share Alike 4.0 International.

**Share-alike działa dalej:** jeżeli udostępnisz przerobione siatki albo dane z nich
wyliczone, muszą iść na tej samej licencji. Dotyczy to również ośmiu plików, które
**są** w repozytorium, bo powstały z pomiarów tych siatek:

- `data/muscles.json` — katalog 477 mięśni: nazwa, strona, oś działania, objętość, środek
- `data/rig_fitted.json` — punkty kostne i linijka *cun* zmierzone na modelu
- `data/structures.json` — 61 drobnych struktur (paliczki, kości śródręcza, małżowina,
  gałka oczna): środek, oś, zakres, obrys
- `data/azimuth_fit.json` — kąty dopasowane przez pomiar odległości do mięśni; zawiera
  ich nazwy i wyniki pomiarów na siatkach Z-Anatomy
- `data/fascia.json` — katalog 198 struktur łącznotkankowych (powięzie, przegrody
  międzymięśniowe, troczki, rozcięgna, pochewki, kaletki): nazwa, klasa, strona,
  normalna płata, grubość, rozpiętość
- `data/septa_analysis.json` — wynik pomiaru odległości punktów i kanałów do śladu
  przegród międzymięśniowych na skórze (`tools/analyze_septa.py`)
- `data/wire.npz` — siatka cunowa: 9144 wierzchołki zmierzone na skórze, z grubością
  tkanki miękkiej, głębokością mięśnia i odległością do powięzi (`tools/bake_wire.py`)
- `data/wire.json` — opis tej siatki wraz z osiami działania mięśni

Wszystkie osiem jest objętych **CC BY-SA 4.0**, nie licencją MIT.

## 3. Opisy punktów — materiał chroniony, **nie ma go w repozytorium**

`data/point_meta.json` (lokalizacje, nakłucie, działanie, wskazania) to OCR książki
*Peter Deadman, „A Manual of Acupuncture"* — materiał chroniony prawem autorskim.
Pliku nie ma w repozytorium i nie wolno go tu dodawać.

Kto ma własny egzemplarz w postaci cyfrowej, generuje go u siebie:

```bash
python tools/import_meta.py <ścieżka/do/points.json>
```

Bez tego pliku aplikacja działa normalnie — panel pokazuje wtedy własne, polskie opisy
lokalizacji z `data/points.json` zamiast cytatów z podręcznika.

## 4. Kotwice punktów — moje

`data/points.json` to moja praca: 361 punktów opisanych kotwicą proporcjonalną
(segment ciała + odległość w *cun* + azymut) wraz z polskimi notatkami o lokalizacji.
Anatomiczne fakty pochodzą z literatury, ale zapis i układ danych są własne — MIT.
