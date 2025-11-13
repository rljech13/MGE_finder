import os
import shutil

# Путь к исходной директории, где лежат исходные папки
input_base_dir = "/home/lam34/MGE_finder/results"
# Путь к целевой директории
output_base_dir = "/home/lam34/MGE_finder/results_2"

# Имена нужных файлов
target_files = {
    "mge_annotated.gbk",
    "integrase_hits_summary.tsv",
    "integrase_trna.tsv",
    "attachment_sites.tsv",
    "mge_region.fa"
}

# Создание выходной директории, если не существует
os.makedirs(output_base_dir, exist_ok=True)

# Проход по всем подпапкам в input_base_dir
for root, dirs, files in os.walk(input_base_dir):
    current_files = set(files)
    if target_files.issubset(current_files):
        sample_name = os.path.basename(root)
        target_dir = os.path.join(output_base_dir, sample_name)
        os.makedirs(target_dir, exist_ok=True)

        # Копируем нужные файлы
        for file in target_files:
            src_path = os.path.join(root, file)
            dst_path = os.path.join(target_dir, file)
            shutil.copy2(src_path, dst_path)
            print(f"✔️ Копирован: {src_path} → {dst_path}")

print("✅ Все подходящие директории обработаны.")