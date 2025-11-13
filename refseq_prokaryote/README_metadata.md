# RefSeq Metadata Extraction

## Описание

Скрипт `extract_metadata.py` извлекает метаданные (ID файла, вид, штамм, тип сборки) из загруженных геномов RefSeq прокариот.

## Использование

### Базовое использование

```bash
# Извлечь метаданные для всех загруженных файлов
python3 extract_metadata.py

# Только бактерии
python3 extract_metadata.py --domain bacteria --output metadata_bacteria.tsv

# Только археи
python3 extract_metadata.py --domain archaea --output metadata_archaea.tsv
```

### Параметры

- `--input_dir` - Базовая директория с поддиректориями `bacteria/` и `archaea/` (по умолчанию: `/home/lam34/MGE_finder/refseq_prokaryote`)
- `--output` - Путь к выходному TSV файлу (по умолчанию: `metadata.tsv`)
- `--domain` - Домен для обработки: `bacteria`, `archaea`, или `both` (по умолчанию: `both`)

## Формат выходного файла

TSV файл со следующими колонками:

- `file_id` - ID файла (например, `bacteria.1.1`)
- `filename` - Имя файла
- `species` - Название вида (род и вид)
- `strain` - Штамм/изолят (если указан)
- `assembly_type` - Тип сборки (WGS, Complete, Chromosome)
- `num_contigs` - Количество контигов
- `path` - Полный путь к файлу

## Примеры

### Просмотр таблицы

```bash
# Просмотр первых 20 строк
head -20 metadata_bacteria.tsv

# Только ID и вид
cut -f1,3 metadata_bacteria.tsv | head -20

# Статистика по видам
tail -n +2 metadata_bacteria.tsv | cut -f3 | sort | uniq -c | sort -rn | head -20
```

### Поиск по виду

```bash
# Найти все файлы для определенного вида
grep "Escherichia coli" metadata_bacteria.tsv

# Найти все файлы для рода
grep "Pseudomonas" metadata_bacteria.tsv
```

## Текущие файлы

- `metadata_bacteria.tsv` - Полная таблица метаданных для бактерий
- `metadata_simple.tsv` - Упрощенная версия (только ID и вид)
- `metadata_final.tsv` - Финальная версия с основными полями

## Обновление метаданных

После загрузки новых файлов запустите скрипт снова для обновления таблицы:

```bash
python3 extract_metadata.py --domain both --output metadata.tsv
```

