# Навигация по проекту IE_finder

Краткая карта: **что где лежит**, что считать «рабочим», а что — черновиками и архивом.

## С чего начать

1. **Обзор модулей и запуск** — корневой [README.md](../README.md), установка — [SETUP.md](../SETUP.md).
2. **Сводные результаты и ноутбуки** — каталог [`results_organized/`](../results_organized/) (см. ниже).
3. **История планов и статусов** — [`docs/archive/`](archive/) (перенесено из корня, чтобы не мешать навигации).

## Каталоги верхнего уровня

| Путь | Назначение |
|------|------------|
| [`IE_finder/`](../IE_finder/) | Основной Snakemake-пайплайн поиска MGE; скрипты в `scripts/`. Подробности: [README.md](../IE_finder/README.md). |
| [`bakta_pipeline/`](../bakta_pipeline/) | Аннотация регионов Bakta + слияние с att-сайтами. |
| [`padloc_pipeline/`](../padloc_pipeline/) | PADLOC (системы защиты). |
| [`results_organized/`](../results_organized/) | **Точка входа для анализа:** объединённые таблицы, GTDB, кластеризация геномов, визуализации, ноутбуки `notebooks/`. Сводные TSV (MGE, PADLOC/defense, метаданные) лежат в [`results_organized/summary/`](../results_organized/summary/). |
| [`results_organized_rep/`](../results_organized_rep/) | Урезанная/репрезентативная выкладка (меньше полного `results_organized`). |
| [`classification_investigation/`](../classification_investigation/) | Кластеризация MGE (MMseqs), vContact, визуализации для классификации. |
| [`classification_investigation_bins/`](../classification_investigation_bins/) | Вариант того же исследования для **бинов** (отдельные скрипты и фигуры). |
| [`protein_clusterization/`](../protein_clusterization/) | MMseqs по белкам Bakta. |
| [`data/`](../data/) | Локальный NCBI-слой: **`ncbi/`** (FASTA), **`metadata/`** (TSV по Deinococcales), **`scripts/`**, **`docs/`**, **`logs/`**. Не путать с `IE_finder/data/`. См. [data/README.md](../data/README.md). |
| [`analysis/`](../analysis/) | Ранний/локальный анализ: `statistics/`, `plots/`, `selected_samples/`, ноутбуки в `notebooks/` (раньше каталог назывался `notebooks ` с пробелом — переименован). |
| [`additionals/`](../additionals/) | Вспомогательные утилиты вне Snakemake (сейчас в основном `collect_reps.py`); не ядро пайплайна. |
| [`genomes/`, `integrases/`, `pfam/`](../) | Небольшие вспомогательные данные и HMM-профили. |

## Большие данные (гигабайты)

Эти каталоги **в `.gitignore`** — это ожидаемые артефакты прогонов, а не «мусор в git».

| Путь | Тип содержимого |
|------|-----------------|
| `IE_finder/results/` | Побочные результаты по каждому геному |
| `IE_finder/data/` | Входные FASTA (в т.ч. `genomes_old/` — вероятная старая копия) |
| `bakta_pipeline/results_bakta_*/` | Выход Bakta (Deinococcales / Thermaceae и т.д.) |
| `padloc_pipeline/results_*/` | CSV/выходы PADLOC |
| `classification_investigation/clusterMGE/results_mge*` | MMseqs и промежуточные FASTA/БД |
| `classification_investigation/vcontact/results_vcontact/` | Выход vContact |
| `results_organized/genome_clustering/` | Тяжёлые профили MMseqs/HMM для ноутбуков кластеризации |

Логи после уборки можно снова накапливать в `*/logs/` и корневом `logs/` — при необходимости их снова можно чистить или ротировать.

## Документация (разбросанные `.md`)

| Где | О чём |
|-----|--------|
| Корень | [README.md](../README.md), [SETUP.md](../SETUP.md), [CONTRIBUTING.md](../CONTRIBUTING.md) |
| `data/` | Статусы загрузок и таксономия |
| `IE_finder/` | [README.md](../IE_finder/README.md) |
| `classification_investigation/` | roadmap, разборы tRNA, графов |
| `classification_investigation/visualizations/` | индекс рисунков |
| `results_organized/notebooks/` | handoff между ноутбуками, кластеризация |
| **`docs/archive/`** | старые планы из корня |

## Скрипты и дубли

### `results_organized/misc/`

Смесь **скриптов** (GTDB/таксономия, Thermaceae, общая статистика, графики) и **тяжёлых артефактов** (например матрицы кластеров белков MGE, Jaccard graphml, `rep_mge_proteins_all.faa`, выгрузки meta\*). Всё это рабочие входы/выходы анализа — не удалять без осознанного решения.

Скрипты под NCBI Thermaceae (в т.ч. `results_organized/misc/plot_*`, `export_*`, `padloc_pipeline/plot_thermaceae_*.py`, `padloc_pipeline/scripts/*thermus*`, `IE_finder/scripts/link_supplement_thermus_genomes.py`) берут корень датасета из **`THERMACEAE_GENOMES_DIR`** (по умолчанию `/home/lam34/Thermaceae_genomes`). Внутри ожидаются `ncbi_dataset/data/data_summary.tsv`, при необходимости `assembly_data_report.jsonl`, и `ncbi_thermus_supplement/ncbi_dataset/data`.

### `results_organized/scripts/`

| Файл | Назначение |
|------|------------|
| `build_results_organized_rep.py` | Сборка `results_organized_rep/` по `merged_deduplicated/rep_index.tsv` |
| `flatten_rep_annotations_only.py` | «Сплющивание» merged GBFF в одном файле на папку rep |
| `build_cut_sites_download_bundle.py` | Бандл для выгрузки cut sites |

### `analysis/` (gitignored)

Ранний слой: `statistics/mge_stats_300000_*`, `selected_samples/` (много однократно скопированных `.gbff`), деревья/алайнменты интеграз, `plots/`, ноутбуки Cytoscape. Не считается основным входом для текущих ноутбуков в `results_organized/notebooks/`.

- **`classification_investigation/visualizations/`** и **`classification_investigation_bins/visualizations_bins/`** — параллельные наборы фигур (геномы vs бины); при правках смотрите, какой набор актуален для публикации.

## Ноутбуки

| Расположение | Содержание |
|--------------|------------|
| [`results_organized/notebooks/`](../results_organized/notebooks/) | Основная цепочка анализа (в т.ч. supervised classification, синтения). |
| [`analysis/notebooks/`](../analysis/notebooks/) | Старые/локальные ноутбуки (Cytoscape, ранние графики). |

---

Если что-то из этого списка устарело после реорганизации — поправьте этот файл одним коммитом вместе с перемещением каталогов.
