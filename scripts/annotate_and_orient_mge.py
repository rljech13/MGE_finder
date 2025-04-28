import argparse
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from logger import Logger

logger = Logger(name="annotate_mge", level=Logger.Level.INFO).get_logger()

def reverse_coords(seq_len: int, start: int, end: int) -> tuple[int, int]:
    """Recalculate coordinates after reverse-complementing a sequence.

    Args:
        seq_len (int): Total sequence length.
        start (int): Original start coordinate (0-based).
        end (int): Original end coordinate (0-based, exclusive).

    Returns:
        tuple[int, int]: New start and end coordinates after reversal.
    """
    return seq_len - end, seq_len - start

def load_input_files(trna_fa_path: str, integrase_path: str, blast_path: str) -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    """Load tRNA FASTA, integrase hit summary, and BLAST hit tables.

    Args:
        trna_fa_path (str): Path to FASTA file of tRNA sequences.
        integrase_path (str): Path to integrase_hits_summary.tsv.
        blast_path (str): Path to mge_blast.tsv.

    Returns:
        tuple: (tRNA records dictionary, integrase DataFrame, BLAST hits DataFrame)
    """
    try:
        # В ключ словаря берем только первую часть идентификатора, чтобы обеспечить сопоставление по core_id.
        trna_records = {rec.id.split(":")[0]: rec for rec in SeqIO.parse(trna_fa_path, "fasta")}
        integrase_df = pd.read_csv(integrase_path, sep="\t")
        blast_df = pd.read_csv(blast_path, sep="\t")
        return trna_records, integrase_df, blast_df
    except pd.errors.EmptyDataError:
        logger.warning("One or more input tables are empty. Writing empty outputs.")
        return {}, pd.DataFrame(), pd.DataFrame()

def get_annotation_data(record_id: str, trna_records: dict, integrase_df: pd.DataFrame, blast_df: pd.DataFrame):
    """Extract relevant annotation data for a single MGE region.

    Args:
        record_id (str): Sequence record ID from the FASTA.
        trna_records (dict): Dictionary of tRNA SeqRecords keyed by ID.
        integrase_df (pd.DataFrame): Integrase hit summary table.
        blast_df (pd.DataFrame): BLAST hits with att site info.

    Returns:
        tuple or None: (core ID, tRNA record, integrase row, BLAST row), or None on failure.
    """
    core_id = record_id.split(":")[0]
    if core_id not in trna_records or core_id not in integrase_df["orf_id"].values:
        return None, None, None, None
    trna_rec = trna_records[core_id]
    int_row = integrase_df[integrase_df["orf_id"] == core_id].iloc[0]
    blast_row = blast_df[blast_df["integrase_id"] == core_id].sort_values("evalue").iloc[0]
    return core_id, trna_rec, int_row, blast_row

def annotate_record(record, trna_rec, int_row, blast_row):
    """Annotate a MGE region with att sites and recalculated integrase CDS.

    If the tRNA (from mge_query) is located on the "+" strand, the sequence region is reverse complemented,
    and the integrase coordinates are recalculated accordingly. In this case, the integrase CDS is set on the
    "-" strand after the reversal. If the tRNA is on the "-" strand, the region is left as is and the integrase CDS
    remains on the "+" strand.

    Args:
        record (SeqRecord): MGE region sequence record. Its id should include region start coordinates
            (e.g. "AE017221.1_1879:AE017221.1:1778425-1793358").
        trna_rec (SeqRecord): tRNA sequence record from mge_query, whose id includes the strand (e.g. "AE017221.1_1879:AE017221.1:1778425-1778501:+").
        int_row (pd.Series): A row from the integrase hits summary table containing the integrase genomic coordinates.
        blast_row (pd.Series): A row from the BLAST hit table used to determine the length of the att site.

    Returns:
        tuple: (updated record, list of SeqFeature annotations, att_data dictionary)
    """
    features = []
    # Исправлено: берем цепь tRNA из trna_rec, а не из record.
    trna_strand = trna_rec.id.split(":")[-1]
    # Извлекаем mge_start из диапазона, записанного в id (например, "1778425-1793358")
    mge_start = int(record.id.split(":")[2].split("-")[0])
    seq_len = len(record.seq)

    # Геномные координаты интегразы из таблицы
    int_start_genome = int(int_row["start"])
    int_end_genome = int(int_row["end"])
    attL_len = int(blast_row["length"])
    trna_len = len(trna_rec.seq)

    # Рассчитываем локальные координаты интегразы относительно начала MGE региона
    local_int_start = int_start_genome - mge_start
    local_int_end = int_end_genome - mge_start

    if trna_strand == "+":
        # Если tRNA на "+" цепи, проводим reverse complement региона
        record.seq = record.seq.reverse_complement()
        # Пересчет координат интегразы: новый интервал = [seq_len - local_int_end, seq_len - local_int_start)
        new_int_start, new_int_end = (seq_len - local_int_end, seq_len - local_int_start)
        # После переворота интегаза отображается на минусовой цепи
        int_feature_strand = -1
    else:
        # Если tRNA на "-" цепи, преобразований не требуется
        new_int_start, new_int_end = local_int_start, local_int_end
        int_feature_strand = 1

    # Определяем координаты att сайтов: attL в начале региона, attR в конце
    attL_start = 0
    attL_end = attL_len
    attR_start = seq_len - trna_len
    attR_end = seq_len

    # Формируем аннотационные features для GenBank
    features.append(
        SeqFeature(FeatureLocation(attL_start, attL_end),
                   type="misc_feature",
                   qualifiers={"note": "attL"})
    )
    features.append(
        SeqFeature(FeatureLocation(attR_start, attR_end),
                   type="misc_feature",
                   qualifiers={"note": "attR"})
    )
    features.append(
        SeqFeature(FeatureLocation(new_int_start, new_int_end, strand=int_feature_strand),
                   type="CDS",
                   qualifiers={"note": "integrase"})
    )

    att_data = {
        "attL_start": attL_start,
        "attL_end": attL_end,
        "attR_start": attR_start,
        "attR_end": attR_end,
        "att_length": attL_end - attL_start
    }

    return record, features, att_data

def annotate_mge(fasta_path: str, trna_fa_path: str, integrase_path: str, blast_path: str,
                 output_gbk: str, output_att: str) -> None:
    """Annotate MGE regions with att sites and integrase CDS.

    Args:
        fasta_path (str): Path to FASTA file of MGE regions.
        trna_fa_path (str): Path to FASTA file of extracted tRNA regions (mge_query).
        integrase_path (str): Path to integrase_hits_summary.tsv file.
        blast_path (str): Path to BLAST hit table.
        output_gbk (str): Path to output GenBank file.
        output_att (str): Path to output att site table (TSV).
    """
    logger.info("Loading input files")
    trna_records, integrase_df, blast_df = load_input_files(trna_fa_path, integrase_path, blast_path)

    if not trna_records or integrase_df.empty or blast_df.empty:
        logger.warning("Missing input data. Writing empty outputs.")
        open(output_gbk, "w").close()
        open(output_att, "w").close()
        return

    annotated = []
    att_table = []

    for record in SeqIO.parse(fasta_path, "fasta"):
        logger.info(f"Annotating {record.id}")
        core_id, trna_rec, int_row, blast_row = get_annotation_data(record.id, trna_records, integrase_df, blast_df)
        if not core_id:
            logger.warning(f"Missing data for {record.id.split(':')[0]}, skipping.")
            continue

        record, features, att_data = annotate_record(record, trna_rec, int_row, blast_row)
        record.features = features
        record.annotations["molecule_type"] = "DNA"
        annotated.append(record)

        att_data["integrase_id"] = core_id
        att_table.append(att_data)

    if annotated:
        logger.info(f"Writing annotated GenBank to {output_gbk}")
        SeqIO.write(annotated, output_gbk, "genbank")
        logger.info(f"Writing att table to {output_att}")
        pd.DataFrame(att_table).to_csv(output_att, sep="\t", index=False)
    else:
        logger.warning("No valid annotations found. Writing empty outputs.")
        open(output_gbk, "w").close()
        open(output_att, "w").close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Annotate MGE regions with att sites and integrase")
    parser.add_argument("--fasta", required=True, help="Path to MGE region FASTA")
    parser.add_argument("--trna_fa", required=True, help="Path to extracted tRNA regions (FASTA)")
    parser.add_argument("--integrase", required=True, help="Path to integrase_hits_summary.tsv")
    parser.add_argument("--blast", required=True, help="Path to mge_blast.tsv")
    parser.add_argument("--out_gbk", required=True, help="Output GenBank file path")
    parser.add_argument("--out_att", required=True, help="Output attachment site table")
    args = parser.parse_args()

    annotate_mge(args.fasta, args.trna_fa, args.integrase, args.blast, args.out_gbk, args.out_att)