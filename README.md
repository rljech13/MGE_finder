
# MGE Finder: Пайплайн для аннотации 

`MGE Finder` — модульный и расширяемый Snakemake-пайплайн для поиска мобильных элем

---

## Структура проекта

```
MGE_finder/
├── config.yaml                  # Главный конфигурационный файл
├── Snakefile                    # Snakemake-описание пайплайна
├── scripts/                     # Все вспомогательные скрипты
│   ├── prepare_fastas.py
│   ├── predict_orfs.py
│   ├── translate_orfs.py
│   └── hmm_search.py
├── logs/                        # Автоматически создаваемые логи
├── results/                     # Все выходные данные пайплайна
└── run.sh                       # Скрипт запуска пайплайна
```

---

## Установка

1. Убедитесь, что у вас установлен [`conda`](https://docs.conda.io/) или [`mamba`](https://mamba.readthedocs.io/).
2. Создайте окружение:

```bash
mamba env create -f envs/MGE_finder.yaml
```

3. Активируйте окружение:

```bash
conda activate MGE_finder
```

---

## Запуск пайплайна

```bash
bash run.sh
```

### Доступные флаги:

- `--dry-run` или `-n`: показать, что будет выполнено, без запуска
- `--unlock`: разблокировать пайплайн после сбоя
- `--force`: принудительно пересчитать цели
- `--rerun-incomplete`: повторить только неполные задания

---

## Что делает пайплайн?

| Шаг | Описание |
|-----|----------|
| **0. prepare_fasta**     | Копирует и конвертирует входные геномные FASTA |
| **1. predict_orfs**      | Предсказывает ORF с помощью Pyrodigal |
| **2. translate_orfs**    | Переводит ORF в белки и сохраняет координаты |
| **3. hmm_search**        | Использует `hmmscan` для поиска интеграз (PF00589 и PF22022) |

---

## Конфигурация (`config.yaml`)

```yaml
paths:
  genomes_dir: "data/genomes"
  results_dir: "results"

execution:
  conda_env: "MGE_finder"

input_sources:
  - "/path/to/ncbi_data"
  - "/path/to/hybracter_data"

pfam_profiles:
  - "pfam/PF00589.27.hmm"
  - "pfam/PF22022.2.hmm"
```

---

## Логи

Все логи автоматически сохраняются в папку `logs/` и дублируются в консоль. Формат логов:

```
[2025-03-24 12:34:56 - INFO]: [predict_orfs] 85 ORFs predicted from: ...
```

---

## Требования

- Python 3.8+
- Snakemake 7.x+
- Biopython
- Pyrodigal
- BCBio.GFF
- HMMER 3

---

## TODO

- [ ] Скачать и начать использовать HMM
- [ ] Вывести статистику для учета количества
- [ ] Внести логику для поиска самих MGE
- [ ] Сделать сливаемым с MCAPP
- [ ] Docker-образ для reproducibility



