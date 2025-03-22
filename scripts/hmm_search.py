import os

sample = snakemake.wildcards.sample
faa_path = snakemake.input.faa
out_path = snakemake.output.hits
hmm_list = snakemake.config["pfam_profiles"]
combined_hmm = f"results/{sample}/combined.hmm"

# Объединяем профили
with open(combined_hmm, "w") as out_hmm:
    for pfam in hmm_list:
        with open(pfam) as f:
            out_hmm.write(f.read())

# Прессуем и запускаем hmmscan
os.system(f"hmmpress {combined_hmm}")
os.system(f"hmmscan --tblout {out_path} {combined_hmm} {faa_path}")