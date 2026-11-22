
import base64
import html
import math
import os
import re
from collections import Counter
import flet as ft


# ==============================================================================
# BIOINFORMATICS REFERENCE DATA & TABLES
# ==============================================================================

IUPAC_DNA_MEANINGS = {
    "A": ("Adenine", {"A": 1.0}),
    "T": ("Thymine", {"T": 1.0}),
    "G": ("Guanine", {"G": 1.0}),
    "C": ("Cytosine", {"C": 1.0}),
    "R": ("Purine (A or G)", {"A": 0.5, "G": 0.5}),
    "Y": ("Pyrimidine (C or T)", {"C": 0.5, "T": 0.5}),
    "S": ("Strong (G or C)", {"G": 0.5, "C": 0.5}),
    "W": ("Weak (A or T)", {"A": 0.5, "T": 0.5}),
    "K": ("Keto (G or T)", {"G": 0.5, "T": 0.5}),
    "M": ("Amino (A or C)", {"A": 0.5, "C": 0.5}),
    "B": ("Not A (C, G, or T)", {"C": 1 / 3, "G": 1 / 3, "T": 1 / 3}),
    "D": ("Not C (A, G, or T)", {"A": 1 / 3, "G": 1 / 3, "T": 1 / 3}),
    "H": ("Not G (A, C, or T)", {"A": 1 / 3, "C": 1 / 3, "T": 1 / 3}),
    "V": ("Not T (A, C, or G)", {"A": 1 / 3, "C": 1 / 3, "G": 1 / 3}),
    "N": ("Any base (A, T, G, or C)", {"A": 0.25, "T": 0.25, "G": 0.25, "C": 0.25}),
}

IUPAC_COMPLEMENT = {
    "A": "T", "T": "A", "U": "A", "G": "C", "C": "G",
    "R": "Y", "Y": "R", "S": "S", "W": "W", "K": "M", "M": "K",
    "B": "V", "V": "B", "D": "H", "H": "D", "N": "N",
}

IUPAC_REGEX = {
    "A": "A", "T": "[TU]", "U": "[TU]", "G": "G", "C": "C",
    "R": "[AG]", "Y": "[CTU]", "S": "[GC]", "W": "[ATU]",
    "K": "[GTU]", "M": "[AC]", "B": "[CGTU]", "D": "[AGTU]",
    "H": "[ACTU]", "V": "[ACG]", "N": "[ATGCU]",
}

DNA_RESIDUE_MW = {"A": 313.21, "T": 304.20, "G": 329.21, "C": 289.18}

CODON_TABLE_STANDARD = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

CODON_TABLE_VERT_MITO = dict(CODON_TABLE_STANDARD)
CODON_TABLE_VERT_MITO.update({"AGA": "*", "AGG": "*", "ATA": "M", "TGA": "W"})

GENETIC_CODES = {
    "1 - Standard Code": (CODON_TABLE_STANDARD, {"ATG"}),
    "2 - Vertebrate Mitochondrial": (CODON_TABLE_VERT_MITO, {"ATG", "ATA", "ATT", "ATC", "GTG"}),
    "11 - Bacterial / Archaeal": (dict(CODON_TABLE_STANDARD), {"ATG", "GTG", "TTG"}),
}

AA_INFO = {
    "A": ("Ala", "Alanine", "Nonpolar", 89.09, 1.8),
    "R": ("Arg", "Arginine", "Basic (+)", 174.20, -4.5),
    "N": ("Asn", "Asparagine", "Polar", 132.12, -3.5),
    "D": ("Asp", "Aspartic acid", "Acidic (-)", 133.10, -3.5),
    "C": ("Cys", "Cysteine", "Polar", 121.16, 2.5),
    "E": ("Glu", "Glutamic acid", "Acidic (-)", 147.13, -3.5),
    "Q": ("Gln", "Glutamine", "Polar", 146.15, -3.5),
    "G": ("Gly", "Glycine", "Nonpolar", 75.07, -0.4),
    "H": ("His", "Histidine", "Basic (+)", 155.16, -3.2),
    "I": ("Ile", "Isoleucine", "Nonpolar", 131.17, 4.5),
    "L": ("Leu", "Leucine", "Nonpolar", 131.17, 3.8),
    "K": ("Lys", "Lysine", "Basic (+)", 146.19, -3.9),
    "M": ("Met", "Methionine", "Nonpolar", 149.21, 1.9),
    "F": ("Phe", "Phenylalanine", "Nonpolar", 165.19, 2.8),
    "P": ("Pro", "Proline", "Nonpolar", 115.13, -1.6),
    "S": ("Ser", "Serine", "Polar", 105.09, -0.8),
    "T": ("Thr", "Threonine", "Polar", 119.12, -0.7),
    "W": ("Trp", "Tryptophan", "Nonpolar", 204.23, -0.9),
    "Y": ("Tyr", "Tyrosine", "Polar", 181.19, -1.3),
    "V": ("Val", "Valine", "Nonpolar", 117.15, 4.2),
    "X": ("Xaa", "Unknown / Ambiguous", "Ambiguous", 110.00, 0.0),
}

PKA_N_TERM = 8.6
PKA_C_TERM = 3.6
PKA_POS = {"K": 10.8, "R": 12.5, "H": 6.5}
PKA_NEG = {"D": 3.9, "E": 4.1, "C": 8.5, "Y": 10.1}

COMMON_MOTIFS = {
    "EcoRI (GAATTC)": "GAATTC",
    "BamHI (GGATCC)": "GGATCC",
    "HindIII (AAGCTT)": "AAGCTT",
    "NotI (GCGGCCGC)": "GCGGCCGC",
    "XbaI (TCTAGA)": "TCTAGA",
    "TaqI (TCGA)": "TCGA",
    "TATA Box (TATAWAW)": "TATAWAW",
    "Kozak (GCCRCCATGG)": "GCCRCCATGG",
    "Poly-A (AATAAA)": "AATAAA",
    "CpG (CG)": "CG",
}

COLOR_PALETTES = {
    "Classic Bioinformatics": {"A": "#2563eb", "T": "#f97316", "G": "#16a34a", "C": "#dc2626", "primary": "#2563eb"},
    "Oceanic Teal & Navy": {"A": "#0b4f6c", "T": "#01baef", "G": "#20bf55", "C": "#757bc8", "primary": "#0b4f6c"},
    "Viridis Scientific": {"A": "#440154", "T": "#31688e", "G": "#35b779", "C": "#eab308", "primary": "#31688e"},
    "Warm Coral & Sunset": {"A": "#e63946", "T": "#f4a261", "G": "#2a9d8f", "C": "#457b9d", "primary": "#e63946"},
    "Emerald & Forest": {"A": "#1b4332", "T": "#2d6a4f", "G": "#52b788", "C": "#74c69d", "primary": "#2d6a4f"},
    "Royal Purple & Gold": {"A": "#5a189a", "T": "#9d4edd", "G": "#ff9e00", "C": "#ff6d00", "primary": "#5a189a"},
    "Monochrome Slate": {"A": "#1e293b", "T": "#475569", "G": "#64748b", "C": "#94a3b8", "primary": "#1e293b"},
    "Uniform Single Color": None,
}

SWATCH_COLORS = [
    "#2563eb", "#0ea5e9", "#0d9488", "#16a34a",
    "#eab308", "#f97316", "#dc2626", "#db2777",
    "#7c3aed", "#1e293b",
]

SAMPLE_FASTA = (
    ">Sample_Gene_V10_Synthetic_Contig | 1,750 bp\n"
    "CTTGGYCYYTGGATGGAAGGGGGATAACTACTGGAAACGGCAGCTAATACCGCATAACGTCGCAAGACCAAAGAGGGGGACCTTCGGGCCTCTTGCCATCGGATGTGCCCMSATGGGATTAGCTAGTAGGTGGGGGTAAAGGCTCACCTAGGCGACGATCCCTAGCTGGTCTGAGAGGATGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGYACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCAGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCTGATACTGGCAAGCTTGAGTCTCGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGGTGAAATGCGTAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACGAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGTCGACTTGGAGGTTGTGCCCTTGAGGCGTGGCTTCCGGAGCTAACGCGTTAAGTCGACCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAATGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCTGGTCTTGACATCCACGGAAGTTTTCAGAGATGAGAATGTGCCTTCGGGAACCGTGAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTTGTGAAATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCCTTTGTTGCCAGCGGTCCGGCCGGGAACTCAAAGGAGACTGCCAGTGATAAACTGGAGGAAGGTGGGGATGACGTCAAGTCATCATGGCCCTTACGACCAGGGCTACACACGTGCTACAATGGCGCATACAAAGAGAAGCGACCTCGCGAGAGCAAGCGGACCTCATAAAGTGCGTCGTAGTCCGGATTGGAGTCTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGTGGATCAGAATGCCACGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTGGGTTGCAAAAGAAGTAGGTAGCTTAACCTTCGGGAGGGCGCTTACCACTTTGTGATTCATGACTGGGGTGAAGTCGTAACAAGGTAACCGTAGGGGAACCTGCGGTTGGATCACCTCCTTACCTTGGYCYYTGGATGGAAGGGGGATAACTACTGGAAACGGCAGCTAATACCGCATAACGTCGCAAGACCAAAGAGGGGGACCTTCGGGCCTCTTGCCATCGGATGTGCCCMSATGGGATTAGCTAGTAGGTGGGGGTAAAGGCTCACCTAGGCGACGATCCCTAGCTGGTCTGAGAGGATGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGYACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGC"
)

SAMPLE_PHYLO_MULTI_FASTA = """\
>Homo_sapiens_COI
ATGTTCGCCGACCGTTGACTATTCTCTACAAACCACAAAGACATTGGAACACTATACCTATTATTCGGCGCATGAGCTGGAGTCCTAGGCACAGCTCTAAGCCTCCTTATTCGAGCCGAGCTGGGCCAGCCAGGCAACCTTCTAGGTAACGACCACATCTACAACGTTATCGTCACAGCCCATGCATTTGTAATAATCTTCTTCATAGTAATACCCATCATAATCGGAGGCTTTGGCAACTGACTAGTTCCCCTAATAATCGGTGCCCCCGATATGGCGTTTCCCCGCATAAACAACATAAGCTTCTGACTCTTACCTCCCTCTCTCCTACTCCTGCTCGCATCTGCTATAGTGGAGGCCGGAGCAGGAACAGGTTGAACAGTCTACCCTCCCTTAGCAGGGAACTACTCCCACCCTGGAGCCTCCGTAGACCTAACCATCTTCTCCTTACACCTAGCAGGTGTCTCCTCTATCTTAGGGGCCATCAATTTCATCACAACAATTATCAATATAAAACCCCCTGCCATAACCCAATACCAAACGCCCCTCTTCGTCTGATCCGTCCTAATCACAGCAGTCCTACTTCTCCTATCTCTCCCAGTCCTAGCTGCTGGCATCACTATACTACTAACAGACCGCAACCTCAACACCACCTTCTTCGACCCCGCCGGAGGAGGAGACCCCATTCTATACCAACACCTATTC
>Pan_troglodytes_COI
ATGTTCGCCGACCGTTGACTATTCTCTACAAACCACAAAGACATTGGAACACTATACCTATTATTCGGCGCATGAGCTGGAGTCCTAGGCACAGCTCTAAGCCTCCTTATTCGAGCCGAGCTGGGCCAGCCAGGCAACCTTCTAGGTAACGACCACATCTACAACGTTATCGTCACAGCCCATGCATTTGTAATAATCTTCTTCATAGTAATACCCATCATAATCGGAGGCTTTGGCAACTGACTAGTTCCCCTAATAATCGGTGCCCCCGATATGGCGTTTCCCCGCATAAACAACATAAGCTTCTGACTCTTACCTCCCTCTCTCCTACTCCTGCTCGCATCTGCTATAGTAGAGGCCGGAGCAGGAACAGGTTGAACAGTCTACCCTCCCTTAGCAGGGAACTACTCCCACCCTGGAGCCTCCGTAGACCTAACCATCTTCTCCTTACACCTAGCAGGTGTCTCCTCTATCTTAGGGGCCATCAATTTCATCACAACAATTATCAATATAAAACCCCCTGCCATAACCCAATACCAAACGCCCCTCTTCGTCTGATCCGTCCTAATCACAGCAGTCCTACTTCTCCTATCTCTCCCAGTCCTAGCTGCTGGCATCACTATACTACTAACAGACCGCAACCTCAACACCACCTTCTTCGACCCCGCCGGAGGAGGAGACCCCATTCTATACCAACACCTATTC
>Gorilla_gorilla_COI
ATGTTCGCCGACCGTTGACTATTCTCTACAAACCACAAAGACATTGGAACACTATACCTATTATTCGGCGCATGAGCTGGAGTCCTAGGCACAGCTCTAAGCCTCCTTATTCGAGCTGAACTGGGCCAGCCAGGCAACCTTCTAGGTAACGACCACATCTACAACGTTATCGTCACAGCCCATGCATTTGTAATAATCTTCTTCATAGTAATACCCATCATAATCGGAGGCTTTGGCAACTGACTAGTTCCCCTAATAATCGGTGCCCCCGATATGGCGTTTCCCCGCATAAACAACATAAGCTTCTGACTCTTACCTCCCTCTCTCCTACTCCTGCTCGCATCTGCTATAGTGGAGGCCGGAGCAGGAACAGGTTGAACAGTCTACCCTCCCTTAGCAGGGAACTACTCCCACCCTGGAGCCTCCGTAGACCTAACCATCTTCTCCTTACACCTAGCAGGTGTCTCCTCTATCTTAGGGGCCATCAATTTCATCACAACAATTATCAATATAAAACCCCCTGCCATAACCCAATACCAAACGCCCCTCTTCGTCTGATCCGTCCTAATCACAGCAGTCCTACTTCTCCTATCTCTCCCAGTCCTAGCTGCTGGCATCACTATACTACTAACAGACCGCAACCTCAACACCACCTTCTTCGACCCCGCCGGAGGAGGAGACCCCATTCTATACCAACACCTATTC
>Pongo_abelii_COI
ATGTTCGCCGACCGTTGACTATTCTCTACAAACCACAAAGACATTGGAACACTATACCTATTATTCGGCGCATGAGCTGGAGTCCTAGGCACAGCCCTAAGCCTCCTTATTCGAGCCGAGCTAGGCCAGCCAGGCAACCTTCTAGGTAACGACCACATCTACAACGTTATCGTCACAGCCCATGCATTTGTAATAATCTTCTTCATAGTAATACCCATCATAATCGGAGGCTTTGGCAACTGACTAGTTCCCCTAATAATCGGTGCCCCCGATATGGCGTTTCCCCGCATAAACAACATAAGCTTCTGACTCTTACCTCCCTCTCTCCTACTCCTGCTCGCATCTGCTATAGTAGAGGCCGGAGCAGGAACAGGTTGAACAGTCTACCCTCCCTTAGCAGGGAACTACTCCCACCCTGGAGCCTCCGTAGACCTAACCATCTTCTCCTTACACCTAGCAGGTGTCTCCTCTATCTTAGGGGCCATCAATTTCATCACAACAATTATCAATATAAAACCCCCTGCCATAACCCAATACCAAACGCCCCTCTTCGTCTGATCCGTCCTAATCACAGCAGTCCTACTTCTCCTATCTCTCCCAGTCCTAGCTGCTGGCATCACTATACTACTAACAGACCGCAACCTCAACACCACCTTCTTCGACCCCGCCGGAGGAGGAGACCCCATTCTATACCAACACCTATTC
>Mus_musculus_COI
ATGTTCATTAATCGTTGATTATTCTCAACCAATCACAAAGATATCGGAACCCTCTATCTACTATTCGGAGCCTGAGCGGGAATAGTGGGTACTGCACTAAGTATTTTAATTCGAGCAGAATTAGGTCAACCAGGTGCACTTTTAGGAGATGACCAAATTTACAATGTTATCGTAACTGCCCATGCTTTTGTTATAATTTTCTTCATAGTAATACCAATAATAATTGGAGGCTTTGGAAACTGACTTGTCCCACTAATAATCGGAGCCCCAGATATAGCATTCCCACGAATAAATAATATAAGTTTTTGACTCCTACCACCATCATTTCTCCTTCTCCTAGCATCATCAATAGTAGAAGCAGGAGCAGGAACAGGATGAACAGTCTACCCACCTCTAGCCGGAAATCTAGCCCATGCAGGAGCATCAGTAGACCTAACAATTTTCTCCCTTCATTTAGCTGGAGTGTCATCTATTTTAGGTGCAATTAATTTTATTACCACTATTATCAACATGAAACCCCCAGCCATAACACAGTATCAAACTCCACTATTTGTCTGATCCGTACTTATTACAGCCGTACTGCTCCTATTATCACTACCAGTGCTAGCCGCAGGCATTACTATACTACTAACAGACCGCAACCTAAACACAACTTTCTTTGATCCCGCTGGAGGAGGGGACCCAATTCTCTACCAGCATCTGTTC
>Gallus_gallus_COI
ATGACCTTCATCAACCGATGATTATTCTCAACCAACCACAAAGACATCGGCACTCTTTACCTAATTTTCGGCACATGGGCGGGCATAGCCGGCACAGCACTTAGCCTTCTAATTCGCGCAGAACTAGGACAGCCCGGAACTCTCTTAGGCGACGATCAAATTTACAATGTAATCGTCACAGCCCATGCTTTCGTCATAATCTTCTTTATAGTTATACCCATCATGATCGGTGGCTTCGGAAACTGACTAGTCCCACTTATAATCGGTGCCCCAGACATAGCATTCCCCCGCATAAATAACATAAGCTTCTGACTCCTCCCTCCCTCATTCCTTCTCCTACTAGCCTCATCTACAGTAGAAGCTGGGGCCGGCACAGGATGAACAGTTTACCCCCCTTTAGCCGGCAACCTAGCCCACGCTGGAGCATCAGTAGACCTAGCCATCTTCTCACTCCACTTAGCAGGTGTTTCCTCCATTCTAGGAGCAATCAACTTTATCACAACCATCATCAACATAAAACCCCCTGCCCTATCACAATACCAAACACCCCTATTCGTATGATCCGTCCTCATTACTGCCATCCTACTACTCCTCTCCCTACCCGTCCTAGCTGCTGGGATTACCATACTACTTACCGACCGCAACCTAAACACCACATTCTTCGACCCAGCTGGAGGAGGAGACCCAATCCTATACCAACACCTATTC
"""


# ==============================================================================
# PURE-PYTHON BIOINFORMATICS ENGINE
# ==============================================================================

def parse_fasta_records(raw_text: str, keep_gaps: bool = False):
    lines = raw_text.strip().splitlines()
    if not lines:
        return []
    records = []
    current_header = "Sequence_1"
    seq_chunks = []
    has_headers = any(line.strip().startswith(">") for line in lines)
    pattern = r"[^A-Za-z\-]" if keep_gaps else r"[^A-Za-z]"

    if not has_headers:
        cleaned = re.sub(pattern, "", raw_text).upper().replace("U", "T")
        return [("Manual_Input_Sequence", cleaned)] if cleaned else []

    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith(">"):
            if seq_chunks:
                records.append((current_header, "".join(seq_chunks)))
                seq_chunks = []
            current_header = stripped[1:].strip() or f"Sequence_{len(records) + 1}"
        else:
            cleaned_line = re.sub(pattern, "", stripped).upper().replace("U", "T")
            seq_chunks.append(cleaned_line)

    if seq_chunks:
        records.append((current_header, "".join(seq_chunks)))
    elif has_headers and not records:
        records.append((current_header, ""))
    return records


def format_fasta(header: str, seq: str, width: int = 70) -> str:
    clean_header = header.lstrip(">").strip() or "DNA_Analyzer_Sequence"
    wrapped = "\n".join(seq[i:i + width] for i in range(0, len(seq), width))
    return f">{clean_header}\n{wrapped}\n"


def validate_dna_sequence(seq: str):
    cleaned = re.sub(r"[^A-Za-z]", "", seq).upper().replace("U", "T")
    valid_chars = set(IUPAC_DNA_MEANINGS.keys())
    invalid = {ch for ch in cleaned if ch not in valid_chars}
    cleaned = "".join(ch for ch in cleaned if ch in valid_chars)
    return cleaned, invalid


def reverse_complement(seq: str) -> str:
    return "".join(IUPAC_COMPLEMENT.get(b, "N") for b in reversed(seq.upper()))


def complement_only(seq: str) -> str:
    return "".join(IUPAC_COMPLEMENT.get(b, "N") for b in seq.upper())


def transcribe_dna_to_rna(seq: str) -> str:
    return seq.upper().replace("T", "U")


def back_transcribe_rna_to_dna(seq: str) -> str:
    return seq.upper().replace("U", "T")


def compute_dna_statistics(seq: str) -> dict:
    n = len(seq)
    counts = Counter(seq)
    a_cnt, t_cnt, g_cnt, c_cnt = counts["A"], counts["T"], counts["G"], counts["C"]
    std_total = a_cnt + t_cnt + g_cnt + c_cnt
    ambig_total = n - std_total

    gc_strict = ((g_cnt + c_cnt) / std_total * 100.0) if std_total > 0 else 0.0
    at_strict = ((a_cnt + t_cnt) / std_total * 100.0) if std_total > 0 else 0.0

    weighted_gc_bases = 0.0
    weighted_mw = 0.0
    for base, cnt in counts.items():
        weights = IUPAC_DNA_MEANINGS.get(base, ("", {"A": 0.25, "T": 0.25, "G": 0.25, "C": 0.25}))[1]
        weighted_gc_bases += cnt * (weights.get("G", 0.0) + weights.get("C", 0.0))
        weighted_mw += cnt * sum(prob * DNA_RESIDUE_MW[b] for b, prob in weights.items())

    gc_weighted = (weighted_gc_bases / n * 100.0) if n > 0 else 0.0
    at_weighted = 100.0 - gc_weighted if n > 0 else 0.0

    ss_mw = weighted_mw + 79.0 if n > 0 else 0.0
    rc_seq = reverse_complement(seq)
    rc_mw = 0.0
    for base, cnt in Counter(rc_seq).items():
        weights = IUPAC_DNA_MEANINGS.get(base, ("", {"A": 0.25, "T": 0.25, "G": 0.25, "C": 0.25}))[1]
        rc_mw += cnt * sum(prob * DNA_RESIDUE_MW[b] for b, prob in weights.items())
    ds_mw = (ss_mw + rc_mw + 79.0) if n > 0 else 0.0

    gc_skew = ((g_cnt - c_cnt) / (g_cnt + c_cnt)) if (g_cnt + c_cnt) > 0 else 0.0
    at_skew = ((a_cnt - t_cnt) / (a_cnt + t_cnt)) if (a_cnt + t_cnt) > 0 else 0.0

    if n == 0:
        tm_basic, tm_method = 0.0, "N/A"
    elif n < 14:
        tm_basic = 2.0 * (a_cnt + t_cnt) + 4.0 * (g_cnt + c_cnt)
        tm_method = "Wallace (<14 bp)"
    else:
        tm_basic = 64.9 + 41.0 * ((g_cnt + c_cnt - 16.4) / n)
        tm_method = "Salt-Adj (>=14 bp)"

    tm_salt = (81.5 + 16.6 * math.log10(0.05) + 0.41 * gc_weighted - (675.0 / n)) if n > 0 else 0.0
    dinucs = Counter(seq[i:i + 2] for i in range(n - 1)) if n > 1 else Counter()
    cpg_count = dinucs.get("CG", 0)
    cpg_oe = ((cpg_count * n) / (c_cnt * g_cnt)) if (n > 1 and c_cnt > 0 and g_cnt > 0) else 0.0

    return {
        "length": n, "counts": counts, "std_total": std_total, "ambig_total": ambig_total,
        "gc_strict": gc_strict, "at_strict": at_strict, "gc_weighted": gc_weighted,
        "at_weighted": at_weighted, "ss_mw": ss_mw, "ds_mw": ds_mw,
        "gc_skew": gc_skew, "at_skew": at_skew, "tm_basic": tm_basic,
        "tm_method": tm_method, "tm_salt": tm_salt, "dinucs": dinucs,
        "cpg_count": cpg_count, "cpg_oe": cpg_oe,
    }


def iupac_to_regex(pattern: str) -> str:
    cleaned = re.sub(r"\s+", "", pattern).upper()
    return "".join(IUPAC_REGEX.get(ch, re.escape(ch)) for ch in cleaned)


def search_motif(seq: str, motif: str, search_reverse: bool = True):
    cleaned_motif = re.sub(r"\s+", "", motif).upper()
    if not cleaned_motif or not seq:
        return []
    compiled = re.compile(f"(?=({iupac_to_regex(cleaned_motif)}))")
    hits = []
    n = len(seq)
    for m in compiled.finditer(seq):
        s0 = m.start()
        ms = m.group(1)
        hits.append(("+", s0 + 1, s0 + len(ms), ms))
    if search_reverse:
        rc_seq = reverse_complement(seq)
        for m in compiled.finditer(rc_seq):
            rs0 = m.start()
            ms = m.group(1)
            re0 = rs0 + len(ms)
            hits.append(("-", n - re0 + 1, n - rs0, ms))
    hits.sort(key=lambda item: (item[1], 0 if item[0] == "+" else 1))
    return hits


def translate_degenerate_codon(codon: str, codon_table: dict, stop_symbol: str = "*") -> str:
    pools = [
        list(IUPAC_DNA_MEANINGS.get(b, ("", {"A": 0.25, "T": 0.25, "G": 0.25, "C": 0.25}))[1].keys())
        for b in codon
    ]
    possible = set()
    for b1 in pools[0]:
        for b2 in pools[1]:
            for b3 in pools[2]:
                aa = codon_table.get(b1 + b2 + b3, "X")
                possible.add(stop_symbol if aa == "*" else aa)
    return next(iter(possible)) if len(possible) == 1 else "X"


def translate_frame(seq: str, frame: int, codon_table: dict, stop_symbol: str = "*") -> str:
    if not seq:
        return ""
    working = seq[frame - 1:] if frame > 0 else reverse_complement(seq)[abs(frame) - 1:]
    protein = []
    for i in range(0, len(working) - 2, 3):
        codon = working[i:i + 3]
        if any(b not in "ATGC" for b in codon):
            aa = translate_degenerate_codon(codon, codon_table, stop_symbol)
        else:
            aa = codon_table.get(codon, "X")
            if aa == "*":
                aa = stop_symbol
        protein.append(aa)
    return "".join(protein)


def compute_protein_properties(protein_seq: str) -> dict:
    clean_aa = [aa for aa in protein_seq.upper() if aa in AA_INFO]
    length = len(clean_aa)
    counts = Counter(clean_aa)
    if length == 0:
        return {
            "length": 0, "mw": 0.0, "pi": 0.0, "gravy": 0.0,
            "aromaticity": 0.0, "ext_reduced": 0, "ext_oxidized": 0,
            "acidic": 0, "basic": 0, "counts": counts,
        }
    mw = sum(AA_INFO[aa][3] for aa in clean_aa) - 18.015 * (length - 1)
    gravy = sum(AA_INFO[aa][4] for aa in clean_aa) / length
    aromatic = sum(counts[aa] for aa in ("F", "W", "Y")) / length * 100.0
    ext_reduced = counts["Y"] * 1490 + counts["W"] * 5500
    ext_oxidized = ext_reduced + (counts["C"] // 2) * 125

    def net_charge(ph: float) -> float:
        q = 1.0 / (1.0 + 10.0 ** (ph - PKA_N_TERM)) - 1.0 / (1.0 + 10.0 ** (PKA_C_TERM - ph))
        for aa, pka in PKA_POS.items():
            if counts[aa]:
                q += counts[aa] / (1.0 + 10.0 ** (ph - pka))
        for aa, pka in PKA_NEG.items():
            if counts[aa]:
                q -= counts[aa] / (1.0 + 10.0 ** (pka - ph))
        return q

    lo, hi = 0.0, 14.0
    for _ in range(50):
        mid = (lo + hi) / 2.0
        if net_charge(mid) > 0:
            lo = mid
        else:
            hi = mid

    return {
        "length": length, "mw": mw, "pi": (lo + hi) / 2.0, "gravy": gravy,
        "aromaticity": aromatic, "ext_reduced": ext_reduced, "ext_oxidized": ext_oxidized,
        "acidic": counts["D"] + counts["E"], "basic": counts["K"] + counts["R"] + counts["H"],
        "counts": counts,
    }


def find_orfs(seq: str, min_bp: int = 75, codon_table: dict = None, start_codons: set = None, strand_mode: str = "Both"):
    if codon_table is None:
        codon_table = CODON_TABLE_STANDARD
    if start_codons is None:
        start_codons = {"ATG"}
    n = len(seq)
    if n < 3:
        return []
    frames = []
    if strand_mode in ("Both", "Forward (+)"):
        frames.extend([1, 2, 3])
    if strand_mode in ("Both", "Reverse (-)"):
        frames.extend([-1, -2, -3])

    rc_seq = reverse_complement(seq)
    orfs = []
    for frame in frames:
        is_fwd = frame > 0
        offset = abs(frame) - 1
        s_seq = seq if is_fwd else rc_seq
        s_len = len(s_seq)
        starts = []
        i = offset
        while i + 3 <= s_len:
            codon = s_seq[i:i + 3]
            aa = codon_table.get(codon, "X")
            if codon in start_codons:
                starts.append((i, codon))
            elif aa == "*" and starts:
                st_idx, st_c = starts[0]
                end_idx = i + 3
                nt_len = end_idx - st_idx
                if nt_len >= min_bp:
                    orf_nt = s_seq[st_idx:end_idx]
                    orf_aa = translate_frame(orf_nt, 1, codon_table, "*")
                    g_start = (st_idx + 1) if is_fwd else (n - end_idx + 1)
                    g_end = end_idx if is_fwd else (n - st_idx)
                    orfs.append({
                        "frame": f"+{frame}" if is_fwd else str(frame),
                        "strand": "+" if is_fwd else "-",
                        "start": g_start, "end": g_end, "nt_len": nt_len,
                        "aa_len": len(orf_aa.rstrip("*")),
                        "start_codon": st_c, "stop_codon": codon,
                        "nt_seq": orf_nt, "aa_seq": orf_aa,
                    })
                starts = []
            i += 3
    orfs.sort(key=lambda x: (-x["nt_len"], x["start"]))
    for idx, item in enumerate(orfs, 1):
        item["id"] = f"ORF_{idx:03d}"
    return orfs


def sliding_window_gc(seq: str, window: int = 50, step: int = 5):
    n = len(seq)
    if n == 0:
        return [], [], [], [], 0
    eff_win = min(max(1, window), n)
    eff_step = max(1, step)
    pos, gc_vals, gc_skews, at_skews = [], [], [], []
    for start in range(0, n - eff_win + 1, eff_step):
        sub = seq[start:start + eff_win]
        g, c, a, t = sub.count("G"), sub.count("C"), sub.count("A"), sub.count("T")
        gc_sum = sum(
            1.0 if ch in ("G", "C", "S") else
            0.5 if ch in ("R", "Y", "K", "M", "N") else
            (2.0 / 3.0) if ch in ("B", "V") else
            (1.0 / 3.0) if ch in ("D", "H") else 0.0
            for ch in sub
        )
        pos.append(start + eff_win // 2)
        gc_vals.append((gc_sum / eff_win) * 100.0)
        gc_skews.append(((g - c) / (g + c)) if (g + c) > 0 else 0.0)
        at_skews.append(((a - t) / (a + t)) if (a + t) > 0 else 0.0)
    if not pos:
        g, c, a, t = seq.count("G"), seq.count("C"), seq.count("A"), seq.count("T")
        pos = [n // 2]
        gc_vals = [((g + c) / n) * 100.0]
        gc_skews = [((g - c) / (g + c)) if (g + c) > 0 else 0.0]
        at_skews = [((a - t) / (a + t)) if (a + t) > 0 else 0.0]
    return pos, gc_vals, gc_skews, at_skews, eff_win


def pairwise_align(seq_a: str, seq_b: str, mode: str = "global", match: int = 2, mismatch: int = -1, gap: int = -2):
    m, n = len(seq_a), len(seq_b)
    if m == 0 or n == 0:
        return None
    score = [[0] * (n + 1) for _ in range(m + 1)]
    trace = [[0] * (n + 1) for _ in range(m + 1)]
    if mode == "global":
        for i in range(1, m + 1):
            score[i][0] = i * gap
            trace[i][0] = 2
        for j in range(1, n + 1):
            score[0][j] = j * gap
            trace[0][j] = 3

    best_score, best_pos = -10**9, (m, n)
    for i in range(1, m + 1):
        ca = seq_a[i - 1]
        row_s, prev_s, row_t = score[i], score[i - 1], trace[i]
        for j in range(1, n + 1):
            cb = seq_b[j - 1]
            s_d = prev_s[j - 1] + (match if ca == cb else mismatch)
            s_u = prev_s[j] + gap
            s_l = row_s[j - 1] + gap
            if mode == "local":
                bc, d = 0, 0
                if s_d > bc:
                    bc, d = s_d, 1
                if s_u > bc:
                    bc, d = s_u, 2
                if s_l > bc:
                    bc, d = s_l, 3
                row_s[j], row_t[j] = bc, d
                if bc >= best_score:
                    best_score, best_pos = bc, (i, j)
            else:
                if s_d >= s_u and s_d >= s_l:
                    row_s[j], row_t[j] = s_d, 1
                elif s_u >= s_l:
                    row_s[j], row_t[j] = s_u, 2
                else:
                    row_s[j], row_t[j] = s_l, 3

    i, j = (m, n) if mode == "global" else best_pos
    final_score = score[m][n] if mode == "global" else best_score
    aln_a, aln_b, mid = [], [], []
    while i > 0 or j > 0:
        if mode == "local" and score[i][j] == 0:
            break
        t = trace[i][j]
        if t == 1:
            ca, cb = seq_a[i - 1], seq_b[j - 1]
            aln_a.append(ca)
            aln_b.append(cb)
            mid.append("|" if ca == cb else ".")
            i -= 1
            j -= 1
        elif t == 2:
            aln_a.append(seq_a[i - 1])
            aln_b.append("-")
            mid.append(" ")
            i -= 1
        elif t == 3:
            aln_a.append("-")
            aln_b.append(seq_b[j - 1])
            mid.append(" ")
            j -= 1
        else:
            break
    aln_a.reverse()
    aln_b.reverse()
    mid.reverse()
    aln_len = len(aln_a)
    matches = sum(1 for c in mid if c == "|")
    mismatches = sum(1 for c in mid if c == ".")
    gaps = aln_len - matches - mismatches
    return {
        "score": final_score, "aligned_len": aln_len, "matches": matches,
        "mismatches": mismatches, "gaps": gaps,
        "identity": (matches / aln_len * 100.0) if aln_len > 0 else 0.0,
        "aln_a": "".join(aln_a), "midline": "".join(mid), "aln_b": "".join(aln_b),
        "start_a": i + 1, "start_b": j + 1,
    }


# ==============================================================================
# PHYLOGENETIC TREE ENGINE
# ==============================================================================

class PhyloNode:
    def __init__(self, name: str = "", branch_length: float = 0.0, children=None):
        self.name = name
        self.branch_length = max(0.0, float(branch_length))
        self.children = list(children) if children else []

    @property
    def is_leaf(self) -> bool:
        return len(self.children) == 0

    def leaves(self):
        if self.is_leaf:
            return [self]
        out = []
        for ch in self.children:
            out.extend(ch.leaves())
        return out

    def leaf_count(self) -> int:
        return len(self.leaves())

    def to_newick(self, is_root: bool = True) -> str:
        safe = re.sub(r"[^A-Za-z0-9_.\-]", "_", self.name) if self.name else ""
        if self.is_leaf:
            return f"{safe}:{self.branch_length:.5f}"
        inner = ",".join(ch.to_newick(is_root=False) for ch in self.children)
        return f"({inner}){safe};" if is_root else f"({inner}){safe}:{self.branch_length:.5f}"


def progressive_msa(records):
    if not records:
        return []
    if len({len(seq) for _, seq in records}) == 1:
        return [(hdr, seq) for hdr, seq in records]
    sorted_idx = sorted(range(len(records)), key=lambda i: len(records[i][1]))
    c_idx = sorted_idx[len(sorted_idx) // 2]
    c_seq = records[c_idx][1][:900]
    pw = {}
    for idx, (_, seq) in enumerate(records):
        if idx != c_idx:
            res = pairwise_align(c_seq, seq[:900], mode="global", match=2, mismatch=-1, gap=-2)
            pw[idx] = (res["aln_a"], res["aln_b"])
    c_len = len(c_seq)
    max_ins = [0] * (c_len + 1)
    ops = {}
    for idx, (ac, aq) in pw.items():
        ins = [[] for _ in range(c_len + 1)]
        ach = ["-"] * c_len
        cp = 0
        for cc, qc in zip(ac, aq):
            if cc == "-":
                ins[cp].append(qc)
            else:
                ach[cp] = qc
                cp += 1
        for p in range(c_len + 1):
            max_ins[p] = max(max_ins[p], len(ins[p]))
        ops[idx] = (ins, ach)

    out_rec = []
    for idx, (hdr, _) in enumerate(records):
        if idx == c_idx:
            s = "".join("-" * max_ins[p] + c_seq[p] for p in range(c_len)) + "-" * max_ins[c_len]
            out_rec.append((hdr, s))
        else:
            ins, ach = ops[idx]
            s = "".join("".join(ins[p]) + "-" * (max_ins[p] - len(ins[p])) + ach[p] for p in range(c_len))
            s += "".join(ins[c_len]) + "-" * (max_ins[c_len] - len(ins[c_len]))
            out_rec.append((hdr, s))
    return out_rec


def compute_evolutionary_distance(seq1: str, seq2: str, model: str) -> float:
    ti_pairs = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
    valid = {"A", "T", "G", "C"}
    comp = diffs = ti = tv = 0
    for b1, b2 in zip(seq1, seq2):
        if b1 in valid and b2 in valid:
            comp += 1
            if b1 != b2:
                diffs += 1
                if (b1, b2) in ti_pairs:
                    ti += 1
                else:
                    tv += 1
    if comp == 0:
        return 1.0
    p = diffs / comp
    if p == 0.0:
        return 0.0
    if "p-distance" in model:
        return p
    if "Jukes-Cantor" in model:
        return max(0.0, -0.75 * math.log(1.0 - (4.0 / 3.0) * min(p, 0.74)))
    P, Q = ti / comp, tv / comp
    t1, t2 = 1.0 - 2.0 * P - Q, 1.0 - 2.0 * Q
    if t1 <= 0.01 or t2 <= 0.01:
        return max(0.0, -0.75 * math.log(1.0 - (4.0 / 3.0) * min(p, 0.74)))
    return max(0.0, -0.5 * math.log(t1 * math.sqrt(t2)))


def build_distance_matrix(aligned_records, model: str):
    n = len(aligned_records)
    seen, names = {}, []
    for hdr, _ in aligned_records:
        base = hdr.split()[0][:24]
        if base in seen:
            seen[base] += 1
            names.append(f"{base}_{seen[base]}")
        else:
            seen[base] = 1
            names.append(base)
    mat = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            d = compute_evolutionary_distance(aligned_records[i][1], aligned_records[j][1], model)
            mat[i][j] = mat[j][i] = d
    return names, mat


def build_nj_tree(names, dist_matrix) -> PhyloNode:
    n = len(names)
    if n == 1:
        return PhyloNode(names[0], 0.0)
    if n == 2:
        d = dist_matrix[0][1] / 2.0
        return PhyloNode("Root", 0.0, [PhyloNode(names[0], d), PhyloNode(names[1], d)])
    clusters = [PhyloNode(nm, 0.0) for nm in names]
    D = [r[:] for r in dist_matrix]
    while len(clusters) > 2:
        m = len(clusters)
        r_sum = [sum(D[i][j] for j in range(m)) for i in range(m)]
        min_q, bi, bj = float("inf"), 0, 1
        for i in range(m):
            for j in range(i + 1, m):
                q = (m - 2) * D[i][j] - r_sum[i] - r_sum[j]
                if q < min_q:
                    min_q, bi, bj = q, i, j
        i, j = bi, bj
        dij = D[i][j]
        li = max(0.0, 0.5 * dij + (r_sum[i] - r_sum[j]) / (2.0 * (m - 2)))
        lj = max(0.0, dij - li)
        clusters[i].branch_length, clusters[j].branch_length = li, lj
        node = PhyloNode("", 0.0, [clusters[i], clusters[j]])
        rem = [k for k in range(m) if k not in (i, j)]
        nd = [max(0.0, 0.5 * (D[i][k] + D[j][k] - dij)) for k in rem]
        clusters = [clusters[k] for k in rem] + [node]
        nm = len(clusters)
        nD = [[0.0] * nm for _ in range(nm)]
        for a, oa in enumerate(rem):
            for b, ob in enumerate(rem):
                nD[a][b] = D[oa][ob]
        for idx, dv in enumerate(nd):
            nD[idx][nm - 1] = nD[nm - 1][idx] = dv
        D = nD
    dl = max(0.0, D[0][1]) / 2.0
    clusters[0].branch_length = clusters[1].branch_length = dl
    return PhyloNode("Root", 0.0, [clusters[0], clusters[1]])


def build_agglomerative_tree(names, dist_matrix, linkage: str) -> PhyloNode:
    n = len(names)
    if n == 1:
        return PhyloNode(names[0], 0.0)
    clusters = [PhyloNode(nm, 0.0) for nm in names]
    heights, sizes = [0.0] * n, [1] * n
    D = [r[:] for r in dist_matrix]
    while len(clusters) > 1:
        m = len(clusters)
        md, bi, bj = float("inf"), 0, 1
        for i in range(m):
            for j in range(i + 1, m):
                if D[i][j] < md:
                    md, bi, bj = D[i][j], i, j
        i, j = bi, bj
        nh = md / 2.0
        clusters[i].branch_length = max(0.0, nh - heights[i])
        clusters[j].branch_length = max(0.0, nh - heights[j])
        node = PhyloNode("", 0.0, [clusters[i], clusters[j]])
        ns = sizes[i] + sizes[j]
        rem = [k for k in range(m) if k not in (i, j)]
        nd = []
        for k in rem:
            dik, djk = D[i][k], D[j][k]
            if linkage == "UPGMA":
                nd.append((sizes[i] * dik + sizes[j] * djk) / ns)
            elif linkage == "WPGMA":
                nd.append(0.5 * (dik + djk))
            elif linkage == "Single-Linkage":
                nd.append(min(dik, djk))
            else:
                nd.append(max(dik, djk))
        clusters = [clusters[k] for k in rem] + [node]
        heights = [heights[k] for k in rem] + [nh]
        sizes = [sizes[k] for k in rem] + [ns]
        nm = len(clusters)
        nD = [[0.0] * nm for _ in range(nm)]
        for a, oa in enumerate(rem):
            for b, ob in enumerate(rem):
                nD[a][b] = D[oa][ob]
        for idx, dv in enumerate(nd):
            nD[idx][nm - 1] = nD[nm - 1][idx] = dv
        D = nD
    return clusters[0]


def fitch_parsimony_score(root: PhyloNode, seq_map: dict) -> int:
    leaves = root.leaves()
    if not leaves:
        return 0
    slen = len(seq_map[leaves[0].name])
    steps = 0

    def post(node: PhyloNode, site: int):
        nonlocal steps
        if node.is_leaf:
            ch = seq_map[node.name][site]
            return {"A", "T", "G", "C"} if ch not in "ATGC" else {ch}
        csets = [post(c, site) for c in node.children]
        inter = set.intersection(*csets)
        if inter:
            return inter
        steps += 1
        return set.union(*csets)

    for s in range(slen):
        post(root, s)
    return steps


def build_maximum_parsimony_tree(names, aligned_records, dist_matrix):
    seq_map = {nm: seq for nm, (_, seq) in zip(names, aligned_records)}
    best_tree = build_nj_tree(names, dist_matrix)
    best_sc = fitch_parsimony_score(best_tree, seq_map)

    def clone(node: PhyloNode) -> PhyloNode:
        return PhyloNode(node.name, node.branch_length, [clone(c) for c in node.children])

    improved, passes = True, 0
    while improved and passes < 6:
        improved = False
        passes += 1
        cand = clone(best_tree)
        edges = []

        def collect(node: PhyloNode):
            for idx, ch in enumerate(node.children):
                if not ch.is_leaf and len(node.children) == 2 and len(ch.children) == 2:
                    edges.append((node, idx, ch))
                collect(ch)

        collect(cand)
        for parent, c_idx, ch in edges:
            o_idx = 1 - c_idx
            for s_idx in (0, 1):
                parent.children[o_idx], ch.children[s_idx] = ch.children[s_idx], parent.children[o_idx]
                sc = fitch_parsimony_score(cand, seq_map)
                if sc < best_sc:
                    best_sc = sc
                    best_tree = clone(cand)
                    improved = True
                else:
                    parent.children[o_idx], ch.children[s_idx] = ch.children[s_idx], parent.children[o_idx]
    return best_tree, best_sc


# ==============================================================================
# CRISP VECTOR SVG CHART & PHYLOGENETIC TREE RENDERERS (100% PURE PYTHON)
# ==============================================================================

def svg_to_data_uri(svg_str: str) -> str:
    b64 = base64.b64encode(svg_str.encode("utf-8")).decode("ascii")
    return f"data:image/svg+xml;base64,{b64}"


def _svg_bar_elements(x0, y0, w, h, labels, values, colors, ylabel, title, show_labels=True, show_grid=True):
    els = [
        f'<text x="{x0 + w/2:.1f}" y="{y0 - 14:.1f}" text-anchor="middle" font-family="sans-serif" font-size="15" font-weight="bold" fill="#0f172a">{html.escape(title)}</text>',
        f'<text x="{x0 - 44:.1f}" y="{y0 + h/2:.1f}" text-anchor="middle" transform="rotate(-90 {x0 - 44:.1f},{y0 + h/2:.1f})" font-family="sans-serif" font-size="12" font-weight="bold" fill="#334155">{html.escape(ylabel)}</text>',
    ]
    max_v = max(values) if values and max(values) > 0 else 10
    top_v = max_v * 1.18
    if show_grid:
        for step in range(5):
            gy = y0 + h - (h * step / 4.0)
            gv = top_v * step / 4.0
            els.append(f'<line x1="{x0}" y1="{gy:.1f}" x2="{x0+w}" y2="{gy:.1f}" stroke="#cbd5e1" stroke-dasharray="4,4" stroke-width="1"/>')
            els.append(f'<text x="{x0 - 8}" y="{gy + 4:.1f}" text-anchor="end" font-family="sans-serif" font-size="11" fill="#475569">{int(round(gv))}</text>')
    els.append(f'<line x1="{x0}" y1="{y0}" x2="{x0}" y2="{y0+h}" stroke="#1e293b" stroke-width="1.5"/>')
    els.append(f'<line x1="{x0}" y1="{y0+h}" x2="{x0+w}" y2="{y0+h}" stroke="#1e293b" stroke-width="1.5"/>')

    n = max(1, len(labels))
    slot_w = w / n
    bar_w = min(slot_w * 0.62, 95)
    total_v = sum(values)
    for idx, (lbl, val, col) in enumerate(zip(labels, values, colors)):
        bh = (val / top_v) * h if top_v > 0 else 0
        bx = x0 + idx * slot_w + (slot_w - bar_w) / 2.0
        by = y0 + h - bh
        els.append(f'<rect x="{bx:.1f}" y="{by:.1f}" width="{bar_w:.1f}" height="{bh:.1f}" fill="{col}" stroke="#0f172a" stroke-width="1" rx="3"/>')
        els.append(f'<text x="{bx + bar_w/2:.1f}" y="{y0 + h + 18:.1f}" text-anchor="middle" font-family="sans-serif" font-size="{11 if n <= 6 else 9}" font-weight="bold" fill="#0f172a">{html.escape(lbl)}</text>')
        if show_labels:
            pct = (val / total_v * 100.0) if total_v > 0 else 0.0
            lbl_txt = f"{val:,} ({pct:.1f}%)" if n <= 5 else f"{val}"
            els.append(f'<text x="{bx + bar_w/2:.1f}" y="{by - 6:.1f}" text-anchor="middle" font-family="sans-serif" font-size="{11 if n <= 6 else 9}" font-weight="bold" fill="#0f172a">{lbl_txt}</text>')
    return "\n".join(els)


def _svg_pie_elements(cx, cy, r, labels, values, colors, title, show_labels=True):
    els = [
        f'<text x="{cx:.1f}" y="{cy - r - 24:.1f}" text-anchor="middle" font-family="sans-serif" font-size="15" font-weight="bold" fill="#0f172a">{html.escape(title)}</text>'
    ]
    total = sum(values)
    if total <= 0:
        els.append(f'<text x="{cx}" y="{cy}" text-anchor="middle" font-family="sans-serif" font-size="13" fill="#475569">No data to display</text>')
        return "\n".join(els)

    angle = -math.pi / 2.0
    non_zero = [(l, v, c) for l, v, c in zip(labels, values, colors) if v > 0]
    for idx, (lbl, val, col) in enumerate(non_zero):
        frac = val / total
        sweep = frac * 2.0 * math.pi
        mid = angle + sweep / 2.0
        ox = math.cos(mid) * 5.0
        oy = math.sin(mid) * 5.0
        if len(non_zero) == 1:
            els.append(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r:.1f}" fill="{col}" stroke="#ffffff" stroke-width="2"/>')
        else:
            x1 = cx + ox + r * math.cos(angle)
            y1 = cy + oy + r * math.sin(angle)
            x2 = cx + ox + r * math.cos(angle + sweep)
            y2 = cy + oy + r * math.sin(angle + sweep)
            large = 1 if sweep > math.pi else 0
            d = f"M {cx+ox:.1f},{cy+oy:.1f} L {x1:.1f},{y1:.1f} A {r:.1f},{r:.1f} 0 {large},1 {x2:.1f},{y2:.1f} Z"
            els.append(f'<path d="{d}" fill="{col}" stroke="#ffffff" stroke-width="2"/>')
        if show_labels and frac >= 0.04:
            tx = cx + ox + (r * 0.62) * math.cos(mid)
            ty = cy + oy + (r * 0.62) * math.sin(mid)
            els.append(f'<text x="{tx:.1f}" y="{ty + 4:.1f}" text-anchor="middle" font-family="sans-serif" font-size="12" font-weight="bold" fill="#ffffff">{frac*100:.1f}%</text>')
        angle += sweep

    # Legend below pie
    ly = cy + r + 24
    col_w = 175
    start_x = cx - min(len(non_zero), 2) * col_w / 2.0 + 10
    for idx, (lbl, val, col) in enumerate(non_zero):
        lx = start_x + (idx % 2) * col_w
        row_y = ly + (idx // 2) * 22
        els.append(f'<rect x="{lx:.1f}" y="{row_y - 11:.1f}" width="13" height="13" rx="2" fill="{col}"/>')
        els.append(f'<text x="{lx + 18:.1f}" y="{row_y:.1f}" font-family="sans-serif" font-size="11.5" font-weight="bold" fill="#1e293b">{html.escape(lbl)}: {val:,} ({val/total*100:.1f}%)</text>')
    return "\n".join(els)


def _svg_line_elements(x0, y0, w, h, xs, ys, color, ylabel, xlabel, title, y_min=None, y_max=None, ref_y=None, ref_label="", show_grid=True):
    els = [
        f'<text x="{x0 + w/2:.1f}" y="{y0 - 14:.1f}" text-anchor="middle" font-family="sans-serif" font-size="14" font-weight="bold" fill="#0f172a">{html.escape(title)}</text>',
        f'<text x="{x0 - 44:.1f}" y="{y0 + h/2:.1f}" text-anchor="middle" transform="rotate(-90 {x0 - 44:.1f},{y0 + h/2:.1f})" font-family="sans-serif" font-size="12" font-weight="bold" fill="#334155">{html.escape(ylabel)}</text>',
        f'<text x="{x0 + w/2:.1f}" y="{y0 + h + 34:.1f}" text-anchor="middle" font-family="sans-serif" font-size="12" font-weight="bold" fill="#334155">{html.escape(xlabel)}</text>',
    ]
    if not xs or not ys:
        return "\n".join(els)
    lo_y = min(ys) if y_min is None else y_min
    hi_y = max(ys) if y_max is None else y_max
    if hi_y <= lo_y:
        hi_y = lo_y + 1.0
    lo_x, hi_x = min(xs), max(xs)
    if hi_x <= lo_x:
        hi_x = lo_x + 1.0

    if show_grid:
        for step in range(5):
            gy = y0 + h - (h * step / 4.0)
            gv = lo_y + (hi_y - lo_y) * step / 4.0
            els.append(f'<line x1="{x0}" y1="{gy:.1f}" x2="{x0+w}" y2="{gy:.1f}" stroke="#cbd5e1" stroke-dasharray="4,4" stroke-width="1"/>')
            els.append(f'<text x="{x0 - 8}" y="{gy + 4:.1f}" text-anchor="end" font-family="sans-serif" font-size="10.5" fill="#475569">{gv:.1f}</text>')
        for step in range(5):
            gx = x0 + (w * step / 4.0)
            xv = lo_x + (hi_x - lo_x) * step / 4.0
            els.append(f'<text x="{gx:.1f}" y="{y0 + h + 16:.1f}" text-anchor="middle" font-family="sans-serif" font-size="10.5" fill="#475569">{int(round(xv))}</text>')

    els.append(f'<line x1="{x0}" y1="{y0}" x2="{x0}" y2="{y0+h}" stroke="#1e293b" stroke-width="1.5"/>')
    els.append(f'<line x1="{x0}" y1="{y0+h}" x2="{x0+w}" y2="{y0+h}" stroke="#1e293b" stroke-width="1.5"/>')

    pts = []
    for x, y in zip(xs, ys):
        px = x0 + ((x - lo_x) / (hi_x - lo_x)) * w
        py = y0 + h - ((y - lo_y) / (hi_y - lo_y)) * h
        pts.append(f"{px:.1f},{py:.1f}")
    els.append(f'<polyline fill="none" stroke="{color}" stroke-width="2.2" points="{" ".join(pts)}"/>')

    if ref_y is not None and lo_y <= ref_y <= hi_y:
        ry = y0 + h - ((ref_y - lo_y) / (hi_y - lo_y)) * h
        els.append(f'<line x1="{x0}" y1="{ry:.1f}" x2="{x0+w}" y2="{ry:.1f}" stroke="#dc2626" stroke-dasharray="6,4" stroke-width="1.5"/>')
        if ref_label:
            els.append(f'<text x="{x0 + w - 6}" y="{ry - 6:.1f}" text-anchor="end" font-family="sans-serif" font-size="11" font-weight="bold" fill="#dc2626">{html.escape(ref_label)}</text>')
    return "\n".join(els)


def render_viz_svg(
    seq: str,
    header: str,
    mode: str,
    base_colors: dict,
    primary_color: str,
    custom_title: str = "",
    win_size: int = 50,
    step_size: int = 5,
    show_labels: bool = True,
    show_grid: bool = True,
    gcode_key: str = "1 - Standard Code",
) -> str:
    bases = ["A", "T", "G", "C"]
    counts = [seq.count(b) for b in bases]
    cols4 = [base_colors[b] for b in bases]

    if mode == "Nucleotide Bar Chart (A, T, G, C)":
        W, H = 760, 480
        body = _svg_bar_elements(
            80, 55, 640, 365, bases, counts, cols4, "Base Count (bp)",
            custom_title or f"Nucleotide Base Composition — {header[:34]}",
            show_labels, show_grid,
        )
    elif mode == "Nucleotide Pie Chart (Base %)":
        W, H = 760, 480
        body = _svg_pie_elements(
            380, 220, 145, bases, counts, cols4,
            custom_title or f"Nucleotide Base Distribution — {header[:34]}",
            show_labels,
        )
    elif mode == "Bar + Pie Combined Dashboard":
        W, H = 760, 780
        b_svg = _svg_bar_elements(
            80, 55, 640, 280, bases, counts, cols4, "Base Count (bp)",
            custom_title or f"Nucleotide Counts (Bar Chart) — {header[:30]}",
            show_labels, show_grid,
        )
        p_svg = _svg_pie_elements(
            380, 540, 130, bases, counts, cols4,
            "Nucleotide Share (Pie Chart)", show_labels,
        )
        body = b_svg + "\n" + p_svg
    elif mode == "Purine/Pyrimidine & GC/AT (Bar + Pie)":
        W, H = 760, 780
        pur = seq.count("A") + seq.count("G")
        pyr = seq.count("C") + seq.count("T")
        gc_c = seq.count("G") + seq.count("C")
        at_c = seq.count("A") + seq.count("T")
        b_svg = _svg_bar_elements(
            80, 55, 640, 280,
            ["Purines (A+G)", "Pyrimidines (C+T)", "GC Pairs", "AT Pairs"],
            [pur, pyr, gc_c, at_c], cols4, "Count (bp)",
            custom_title or "Structural Base Classes (Bar Chart)",
            show_labels, show_grid,
        )
        p_svg = _svg_pie_elements(
            380, 540, 130, ["GC Pairs", "AT Pairs"], [gc_c, at_c],
            [base_colors["G"], base_colors["A"]], "GC vs. AT Ratio (Pie Chart)", show_labels,
        )
        body = b_svg + "\n" + p_svg
    elif mode == "Sliding Window GC Content (%)":
        W, H = 760, 460
        pos, gc_vals, _, _, eff_w = sliding_window_gc(seq, win_size, step_size)
        mean_gc = sum(gc_vals) / len(gc_vals) if gc_vals else 0.0
        body = _svg_line_elements(
            80, 55, 640, 340, pos, gc_vals, primary_color, "GC Content (%)",
            "Sequence Position (bp)",
            custom_title or f"Sliding Window GC Content ({eff_w} bp window)",
            0.0, 100.0, mean_gc, f"Mean GC = {mean_gc:.2f}%", show_grid,
        )
    elif mode == "GC & AT Skew Profile":
        W, H = 760, 680
        pos, _, gc_sk, at_sk, eff_w = sliding_window_gc(seq, win_size, step_size)
        s1 = _svg_line_elements(
            80, 50, 640, 240, pos, gc_sk, base_colors["G"], "GC Skew",
            "Position (bp)", custom_title or f"GC Skew (G-C)/(G+C) [{eff_w} bp]",
            -1.0, 1.0, 0.0, "", show_grid,
        )
        s2 = _svg_line_elements(
            80, 380, 640, 240, pos, at_sk, base_colors["T"], "AT Skew",
            "Position (bp)", f"AT Skew (A-T)/(A+T) [{eff_w} bp]",
            -1.0, 1.0, 0.0, "", show_grid,
        )
        body = s1 + "\n" + s2
    elif mode == "Cumulative GC Skew (OriC Finder)":
        W, H = 760, 460
        cum, run = [], 0.0
        for ch in seq:
            if ch == "G":
                run += 1.0
            elif ch == "C":
                run -= 1.0
            cum.append(run)
        xs = list(range(1, len(seq) + 1))
        body = _svg_line_elements(
            80, 55, 640, 340, xs, cum, primary_color, "Cumulative (G - C)",
            "Position (bp)", custom_title or "Cumulative GC Skew Profile (OriC / Ter)",
            None, None, None, "", show_grid,
        )
    elif mode == "Dinucleotide Frequencies (Bar Chart)":
        W, H = 760, 480
        dinucs = Counter(seq[i:i + 2] for i in range(len(seq) - 1) if seq[i] in "ATGC" and seq[i + 1] in "ATGC")
        all16 = [b1 + b2 for b1 in "ATGC" for b2 in "ATGC"]
        vals = [dinucs.get(d, 0) for d in all16]
        dcols = [base_colors[d[0]] for d in all16]
        body = _svg_bar_elements(
            75, 55, 650, 365, all16, vals, dcols, "Frequency",
            custom_title or "All 16 Canonical Dinucleotides (Bar Chart)",
            show_labels, show_grid,
        )
    else:  # Amino Acid Class (Bar + Pie)
        W, H = 760, 780
        ctable, _ = GENETIC_CODES[gcode_key]
        pep = translate_frame(seq, 1, ctable, "").replace("X", "")
        cc = Counter(AA_INFO[aa][2] for aa in pep if aa in AA_INFO)
        clabels = ["Nonpolar", "Polar", "Basic (+)", "Acidic (-)"]
        cvals = [cc.get(c, 0) for c in clabels]
        b_svg = _svg_bar_elements(
            80, 55, 640, 280, clabels, cvals, cols4, "Residue Count",
            custom_title or "Frame +1 Amino Acid Classes (Bar)",
            show_labels, show_grid,
        )
        p_svg = _svg_pie_elements(
            380, 540, 130, clabels, cvals, cols4,
            "Amino Acid Class Share (Pie)", show_labels,
        )
        body = b_svg + "\n" + p_svg

    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">'
        f'<rect width="100%" height="100%" fill="#ffffff" rx="8"/>'
        f'{body}</svg>'
    )


def render_dotplot_svg(seq_a: str, seq_b: str) -> str:
    a_sub, b_sub = seq_a[:500], seq_b[:500]
    k = 6 if min(len(a_sub), len(b_sub)) >= 180 else 4
    idx_b = {}
    for j in range(len(b_sub) - k + 1):
        idx_b.setdefault(b_sub[j:j + k], []).append(j)
    W, H = 560, 560
    x0, y0, size = 65, 50, 450
    la, lb = max(1, len(a_sub)), max(1, len(b_sub))
    dots = []
    for i in range(len(a_sub) - k + 1):
        km = a_sub[i:i + k]
        if km in idx_b:
            px = x0 + (i / la) * size
            for j in idx_b[km]:
                py = y0 + (j / lb) * size
                dots.append(f'<circle cx="{px:.1f}" cy="{py:.1f}" r="1.8" fill="#2563eb"/>')
                if len(dots) >= 4500:
                    break
        if len(dots) >= 4500:
            break
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">'
        f'<rect width="100%" height="100%" fill="#ffffff"/>'
        f'<text x="{W/2}" y="28" text-anchor="middle" font-family="sans-serif" font-size="14" font-weight="bold" fill="#0f172a">Sequence Dot Plot (k={k}, {len(dots)} hits)</text>'
        f'<rect x="{x0}" y="{y0}" width="{size}" height="{size}" fill="#f8fafc" stroke="#1e293b" stroke-width="1.5"/>'
        f'<text x="{x0 + size/2}" y="{y0 + size + 28}" text-anchor="middle" font-family="sans-serif" font-size="12" font-weight="bold" fill="#334155">Main Sequence A (1..{len(a_sub)} bp)</text>'
        f'<text x="20" y="{y0 + size/2}" text-anchor="middle" transform="rotate(-90 20,{y0 + size/2})" font-family="sans-serif" font-size="12" font-weight="bold" fill="#334155">Reference Sequence B (1..{len(b_sub)} bp)</text>'
        f'{"".join(dots)}</svg>'
    )


def render_phylo_svg(tree: PhyloNode, method: str, style: str, show_len: bool = True) -> str:
    leaves = tree.leaves()
    n_leaves = max(1, len(leaves))
    y_coords = {lf: float(i) for i, lf in enumerate(leaves)}

    def assign_y(nd: PhyloNode) -> float:
        if nd.is_leaf:
            return y_coords[nd]
        yv = sum(assign_y(c) for c in nd.children) / len(nd.children)
        y_coords[nd] = yv
        return yv

    assign_y(tree)
    is_clado = "Cladogram" in style

    def max_depth(nd: PhyloNode) -> int:
        return 0 if nd.is_leaf else 1 + max(max_depth(c) for c in nd.children)

    t_levels = max(1, max_depth(tree))
    x_coords = {}

    def assign_x_phylo(nd: PhyloNode, cx: float):
        x_coords[nd] = cx
        for c in nd.children:
            assign_x_phylo(c, cx + (c.branch_length if c.branch_length > 0 else 0.002))

    def assign_x_clado(nd: PhyloNode, lvl: int):
        if nd.is_leaf:
            x_coords[nd] = float(t_levels)
        else:
            x_coords[nd] = float(lvl)
            for c in nd.children:
                assign_x_clado(c, lvl + 1)

    assign_x_clado(tree, 0) if is_clado else assign_x_phylo(tree, 0.0)
    max_x = max(x_coords.values()) if x_coords else 1.0
    if max_x <= 0:
        max_x = 1.0

    els = []
    if "Rectangular" in style:
        W = 760
        H = max(420, 65 * n_leaves + 90)
        x0, y0, pw, ph = 40, 50, 480, H - 95

        def sx(val):
            return x0 + (val / max_x) * pw

        def sy(val):
            return y0 + (val / max(1, n_leaves - 1)) * ph if n_leaves > 1 else y0 + ph / 2.0

        def draw_rect(nd: PhyloNode):
            nx0, ny0 = sx(x_coords[nd]), sy(y_coords[nd])
            if not nd.is_leaf:
                cys = [sy(y_coords[c]) for c in nd.children]
                els.append(f'<line x1="{nx0:.1f}" y1="{min(cys):.1f}" x2="{nx0:.1f}" y2="{max(cys):.1f}" stroke="#1e40af" stroke-width="2.2"/>')
                for c in nd.children:
                    nx1, ny1 = sx(x_coords[c]), sy(y_coords[c])
                    els.append(f'<line x1="{nx0:.1f}" y1="{ny1:.1f}" x2="{nx1:.1f}" y2="{ny1:.1f}" stroke="#1e40af" stroke-width="2.2"/>')
                    if show_len and not is_clado and c.branch_length > 0:
                        els.append(f'<text x="{(nx0+nx1)/2:.1f}" y="{ny1 - 5:.1f}" text-anchor="middle" font-family="sans-serif" font-size="10" fill="#475569">{c.branch_length:.4f}</text>')
                    draw_rect(c)
            else:
                els.append(f'<circle cx="{nx0:.1f}" cy="{ny0:.1f}" r="4.5" fill="#ea580c"/>')
                els.append(f'<text x="{nx0 + 10:.1f}" y="{ny0 + 4:.1f}" font-family="sans-serif" font-size="12.5" font-weight="bold" fill="#0f172a">{html.escape(nd.name)}</text>')

        draw_rect(tree)
        els.append(f'<text x="{W/2}" y="28" text-anchor="middle" font-family="sans-serif" font-size="14" font-weight="bold" fill="#0f172a">{html.escape(method)} — {html.escape(style)}</text>')
    elif "Circular" in style:
        W, H = 720, 720
        cx, cy, max_r = 360, 370, 210
        angles = {nd: (2.0 * math.pi * y_coords[nd]) / n_leaves for nd in x_coords}

        def draw_circ(nd: PhyloNode):
            r0 = (x_coords[nd] / max_x) * max_r
            t0 = angles[nd]
            if not nd.is_leaf:
                cts = [angles[c] for c in nd.children]
                tmin, tmax = min(cts), max(cts)
                pts = [f"{cx + r0*math.cos(tmin + (tmax-tmin)*k/24):.1f},{cy + r0*math.sin(tmin + (tmax-tmin)*k/24):.1f}" for k in range(25)]
                els.append(f'<polyline fill="none" stroke="#1e40af" stroke-width="2" points="{" ".join(pts)}"/>')
                for c in nd.children:
                    r1 = (x_coords[c] / max_x) * max_r
                    t1 = angles[c]
                    els.append(f'<line x1="{cx + r0*math.cos(t1):.1f}" y1="{cy + r0*math.sin(t1):.1f}" x2="{cx + r1*math.cos(t1):.1f}" y2="{cy + r1*math.sin(t1):.1f}" stroke="#1e40af" stroke-width="2"/>')
                    draw_circ(c)
            else:
                lx, ly = cx + r0 * math.cos(t0), cy + r0 * math.sin(t0)
                tx, ty = cx + (r0 + 14) * math.cos(t0), cy + (r0 + 14) * math.sin(t0)
                anc = "start" if math.cos(t0) >= 0 else "end"
                els.append(f'<circle cx="{lx:.1f}" cy="{ly:.1f}" r="4.5" fill="#ea580c"/>')
                els.append(f'<text x="{tx:.1f}" y="{ty + 4:.1f}" text-anchor="{anc}" font-family="sans-serif" font-size="11.5" font-weight="bold" fill="#0f172a">{html.escape(nd.name)}</text>')

        draw_circ(tree)
        els.append(f'<text x="{W/2}" y="28" text-anchor="middle" font-family="sans-serif" font-size="14" font-weight="bold" fill="#0f172a">{html.escape(method)} — Circular Layout</text>')
    else:  # Unrooted Radial
        W, H = 720, 680
        pos2d = {}

        def layout_unrooted(nd: PhyloNode, x: float, y: float, a0: float, a1: float):
            pos2d[nd] = (x, y)
            if nd.is_leaf:
                return
            tl = nd.leaf_count()
            ca = a0
            for c in nd.children:
                span = (a1 - a0) * (c.leaf_count() / max(1, tl))
                ma = ca + span / 2.0
                ln = c.branch_length if c.branch_length > 0 else 0.005
                layout_unrooted(c, x + ln * math.cos(ma), y + ln * math.sin(ma), ca, ca + span)
                ca += span

        layout_unrooted(tree, 0.0, 0.0, 0.0, 2.0 * math.pi)
        xs = [p[0] for p in pos2d.values()]
        ys = [p[1] for p in pos2d.values()]
        min_x, max_ux = min(xs), max(xs)
        min_y, max_uy = min(ys), max(ys)
        sp = max(max_ux - min_x, max_uy - min_y, 0.01)

        def map_pt(pt):
            return 130 + ((pt[0] - min_x) / sp) * 440, 110 + ((pt[1] - min_y) / sp) * 440

        def draw_unrooted(nd: PhyloNode):
            x0, y0 = map_pt(pos2d[nd])
            for c in nd.children:
                x1, y1 = map_pt(pos2d[c])
                els.append(f'<line x1="{x0:.1f}" y1="{y0:.1f}" x2="{x1:.1f}" y2="{y1:.1f}" stroke="#1e40af" stroke-width="2.2"/>')
                draw_unrooted(c)
            if nd.is_leaf:
                els.append(f'<circle cx="{x0:.1f}" cy="{y0:.1f}" r="4.5" fill="#ea580c"/>')
                els.append(f'<text x="{x0 + 8:.1f}" y="{y0 + 4:.1f}" font-family="sans-serif" font-size="12" font-weight="bold" fill="#0f172a">{html.escape(nd.name)}</text>')

        draw_unrooted(tree)
        els.append(f'<text x="{W/2}" y="28" text-anchor="middle" font-family="sans-serif" font-size="14" font-weight="bold" fill="#0f172a">{html.escape(method)} — Unrooted Radial</text>')

    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">'
        f'<rect width="100%" height="100%" fill="#ffffff" rx="8"/>'
        f'{"".join(els)}</svg>'
    )


# ==============================================================================
# FLET ANDROID MOBILE APPLICATION (ALL 7 MODULES + FULL-SCREEN GRAPH STUDIO)
# ==============================================================================

def main(page: ft.Page):
    page.title = "DNA Analyzer V10"
    page.theme_mode = ft.ThemeMode.LIGHT
    page.bgcolor = "#f1f5f9"
    page.padding = 0

    state = {
        "header": "Sequence_1",
        "records": [],
        "active_tab": 0,
        "in_viz_studio": False,
        "base_colors": {"A": "#2563eb", "T": "#f97316", "G": "#16a34a", "C": "#dc2626"},
        "primary_color": "#2563eb",
        "last_viz_svg": "",
        "last_phylo_tree": None,
        "last_phylo_svg": "",
        "last_newick": "",
        "detected_orfs": [],
        "pending_export_bytes": b"",
    }

    status_txt = ft.Text("Ready — Enter or upload a FASTA sequence.", size=12, color="#334155")
    len_badge = ft.Container(
        content=ft.Text("0 bp", size=12, weight=ft.FontWeight.BOLD, color="#ffffff"),
        bgcolor="#1e40af",
        padding=ft.Padding(10, 4, 10, 4),
        border_radius=12,
    )

    def notify(msg: str):
        status_txt.value = msg
        try:
            sb = ft.SnackBar(content=ft.Text(msg), duration=2200)
            if hasattr(page, "show_dialog"):
                page.show_dialog(sb)
            elif hasattr(page, "open"):
                page.open(sb)
            else:
                page.overlay.append(sb)
                sb.open = True
        except Exception:
            pass
        page.update()

    def copy_text(text: str, label: str = "Text"):
        if not text:
            notify("Nothing to copy.")
            return
        try:
            if hasattr(page, "set_clipboard"):
                page.set_clipboard(text)
            else:
                page.clipboard = text
            notify(f"Copied {label} to clipboard.")
        except Exception:
            notify(f"{label} ready in box.")

    # Native Android File Pickers (Flet 1.x Async Service API)
    main_picker = ft.FilePicker()
    phylo_picker = ft.FilePicker()
    save_picker = ft.FilePicker()
    page.services.extend([main_picker, phylo_picker, save_picker])

    async def pick_main_fasta(_e=None):
        try:
            files = await main_picker.pick_files(
                allow_multiple=False,
                with_data=True,
                file_type=ft.FilePickerFileType.ANY,
            )
        except Exception as ex:
            notify(f"File picker error: {ex}")
            return
        if not files:
            return
        f = files[0]
        raw_str = ""
        if getattr(f, "bytes", None):
            raw_str = f.bytes.decode("utf-8", errors="replace")
        elif getattr(f, "path", None) and os.path.exists(f.path):
            with open(f.path, "r", encoding="utf-8", errors="replace") as fp:
                raw_str = fp.read()
        recs = parse_fasta_records(raw_str)
        if not recs or not recs[0][1]:
            notify("No valid DNA bases found in file.")
            return
        state["records"] = recs
        state["header"] = recs[0][0]
        seq_input.value = recs[0][1]
        update_seq_len()
        if len(recs) >= 3:
            phylo_input.value = raw_str.strip()
        refresh_active_module()
        notify(f"Loaded '{f.name}' ({len(recs)} record(s), {len(recs[0][1]):,} bp).")

    async def pick_phylo_fasta(_e=None):
        try:
            files = await phylo_picker.pick_files(
                allow_multiple=False,
                with_data=True,
                file_type=ft.FilePickerFileType.ANY,
            )
        except Exception as ex:
            notify(f"File picker error: {ex}")
            return
        if not files:
            return
        f = files[0]
        raw_str = ""
        if getattr(f, "bytes", None):
            raw_str = f.bytes.decode("utf-8", errors="replace")
        elif getattr(f, "path", None) and os.path.exists(f.path):
            with open(f.path, "r", encoding="utf-8", errors="replace") as fp:
                raw_str = fp.read()
        recs = parse_fasta_records(raw_str, keep_gaps=True)
        if len(recs) < 3:
            notify("Need >= 3 FASTA sequences for a phylogenetic tree.")
            return
        phylo_input.value = raw_str.strip()
        run_build_tree(None)

    async def export_text_or_bytes(filename: str, data_bytes: bytes):
        state["pending_export_bytes"] = data_bytes
        try:
            res = await save_picker.save_file(file_name=filename, src_bytes=data_bytes, file_type=ft.FilePickerFileType.ANY)
            if res:
                notify(f"Exported {filename}")
        except Exception:
            notify(f"Use Copy button to share {filename}.")

    # --------------------------------------------------------------------------
    # SEQUENCE INPUT & VALIDATION HELPERS
    # --------------------------------------------------------------------------
    def update_seq_len():
        recs = parse_fasta_records(seq_input.value or "")
        s = recs[0][1] if recs else ""
        len_badge.content.value = f"{len(s):,} bp"

    def get_clean_seq(warn: bool = True) -> str:
        raw = (seq_input.value or "").strip()
        if not raw:
            if warn:
                notify("Enter or upload a DNA sequence first.")
            return ""
        recs = parse_fasta_records(raw)
        if recs:
            hdr, seq = recs[0]
            if raw.lstrip().startswith(">"):
                state["header"] = hdr
        else:
            seq = raw
        cleaned, invalid = validate_dna_sequence(seq)
        if invalid:
            notify(f"Ignored non-IUPAC chars: {', '.join(sorted(invalid))}")
        if not cleaned and warn:
            notify("No valid DNA/IUPAC bases found.")
        return cleaned

    # Main Sequence Input Box (Starts completely empty!)
    seq_input = ft.TextField(
        label="DNA / FASTA Sequence Input (Starts Empty)",
        hint_text="Paste DNA or FASTA sequence here, or tap Upload FASTA / Sample...",
        multiline=True,
        min_lines=4,
        max_lines=6,
        text_style=ft.TextStyle(font_family="monospace", size=13),
        bgcolor="#ffffff",
        border_radius=8,
        on_change=lambda _e: (update_seq_len(), page.update()),
    )

    def load_sample_dna(_e):
        recs = parse_fasta_records(SAMPLE_FASTA)
        state["records"] = recs
        state["header"] = recs[0][0]
        seq_input.value = recs[0][1]
        update_seq_len()
        refresh_active_module()
        notify(f"Loaded sample sequence ({len(recs[0][1]):,} bp).")

    def clear_all(_e):
        seq_input.value = ""
        state["records"] = []
        state["header"] = "Sequence_1"
        update_seq_len()
        stats_summary_col.controls.clear()
        stats_table_col.controls.clear()
        tool_out_field.value = ""
        motif_hits_col.controls.clear()
        protein_out_field.value = ""
        protein_props_col.controls.clear()
        orf_list_col.controls.clear()
        orf_detail_field.value = ""
        aln_out_field.value = ""
        notify("Cleared input and results.")

    # --------------------------------------------------------------------------
    # MODULE 1: DNA STATISTICS
    # --------------------------------------------------------------------------
    stats_summary_col = ft.Column(spacing=6)
    stats_table_col = ft.Column(spacing=4)

    def run_analyze_dna(_e=None):
        seq = get_clean_seq(warn=True)
        if not seq:
            return
        st = compute_dna_statistics(seq)
        metrics = [
            ("Record Header", state["header"]),
            ("Total Length", f"{st['length']:,} bp (Standard: {st['std_total']:,}, IUPAC: {st['ambig_total']:,})"),
            ("GC / AT Content (Strict)", f"GC: {st['gc_strict']:.2f}%   |   AT: {st['at_strict']:.2f}%"),
            ("GC / AT Content (IUPAC)", f"GC: {st['gc_weighted']:.2f}%   |   AT: {st['at_weighted']:.2f}%"),
            ("GC Skew / AT Skew", f"GC Skew: {st['gc_skew']:+.4f}   |   AT Skew: {st['at_skew']:+.4f}"),
            ("CpG Islands (Obs/Exp)", f"Count: {st['cpg_count']:,}   |   O/E Ratio: {st['cpg_oe']:.4f}"),
            ("ssDNA Molecular Weight", f"{st['ss_mw']:,.2f} Da ({st['ss_mw']/1000:,.2f} kDa)"),
            ("dsDNA Molecular Weight", f"{st['ds_mw']:,.2f} Da ({st['ds_mw']/1000:,.2f} kDa)"),
            ("Melting Temp (Tm)", f"{st['tm_basic']:.2f} °C [{st['tm_method']}] | Salt: {st['tm_salt']:.2f} °C"),
        ]
        stats_summary_col.controls = [
            ft.Container(
                content=ft.Column([
                    ft.Text(k, size=11, color="#475569", weight=ft.FontWeight.BOLD),
                    ft.Text(v, size=13, color="#0f172a", weight=ft.FontWeight.W_600),
                ], spacing=1),
                bgcolor="#f8fafc",
                padding=8,
                border_radius=6,
                border=ft.Border.all(1, "#e2e8f0"),
            )
            for k, v in metrics
        ]

        rows = []
        ordered = ["A", "T", "G", "C"] + [b for b in IUPAC_DNA_MEANINGS if b not in ("A", "T", "G", "C")]
        for b in ordered:
            cnt = st["counts"].get(b, 0)
            if b in ("A", "T", "G", "C") or cnt > 0:
                pct = (cnt / st["length"] * 100.0) if st["length"] > 0 else 0.0
                rows.append(
                    ft.Container(
                        content=ft.Row([
                            ft.Text(b, width=32, weight=ft.FontWeight.BOLD, color="#1e40af"),
                            ft.Text(IUPAC_DNA_MEANINGS[b][0], expand=True, size=12),
                            ft.Text(f"{cnt:,}", width=70, text_align=ft.TextAlign.RIGHT, size=12),
                            ft.Text(f"{pct:.2f}%", width=65, text_align=ft.TextAlign.RIGHT, weight=ft.FontWeight.BOLD, size=12),
                        ]),
                        padding=ft.Padding(8, 6, 8, 6),
                        bgcolor="#ffffff",
                        border=ft.Border.all(1, "#e2e8f0"),
                        border_radius=6,
                    )
                )
        stats_table_col.controls = rows
        notify(f"DNA statistics calculated ({st['length']:,} bp).")

    async def export_stats_csv(_e):
        seq = get_clean_seq(warn=True)
        if not seq:
            return
        st = compute_dna_statistics(seq)
        lines = [
            "Metric,Value",
            f"Header,{state['header']}",
            f"Length_bp,{st['length']}",
            f"GC_Strict_Pct,{st['gc_strict']:.4f}",
            f"GC_IUPAC_Pct,{st['gc_weighted']:.4f}",
            f"GC_Skew,{st['gc_skew']:.4f}",
            f"AT_Skew,{st['at_skew']:.4f}",
            f"ssDNA_MW_Da,{st['ss_mw']:.2f}",
            f"dsDNA_MW_Da,{st['ds_mw']:.2f}",
            f"Tm_C,{st['tm_basic']:.2f}",
        ]
        await export_text_or_bytes("dna_statistics.csv", "\n".join(lines).encode("utf-8"))

    view_stats = ft.Column([
        ft.Row([
            ft.Button("Calculate Stats", icon=ft.Icons.ANALYTICS, on_click=run_analyze_dna),
            ft.OutlinedButton("Export CSV", icon=ft.Icons.DOWNLOAD, on_click=export_stats_csv),
        ], wrap=True),
        ft.Text("Sequence Properties & Thermodynamics", weight=ft.FontWeight.BOLD, size=14),
        stats_summary_col,
        ft.Text("Nucleotide & IUPAC Base Composition", weight=ft.FontWeight.BOLD, size=14),
        stats_table_col,
    ], spacing=10)

    # --------------------------------------------------------------------------
    # MODULE 2: SEQUENCE TOOLS & IUPAC MOTIF FINDER
    # --------------------------------------------------------------------------
    tool_out_field = ft.TextField(
        label="Transformed Output",
        multiline=True,
        min_lines=4,
        max_lines=6,
        text_style=ft.TextStyle(font_family="monospace", size=12),
        bgcolor="#ffffff",
    )
    motif_field = ft.TextField(label="IUPAC Motif Pattern", value="GAATTC", width=180, bgcolor="#ffffff")
    motif_rev_chk = ft.Checkbox(label="Reverse Strand (-)", value=True)
    motif_hits_col = ft.Column(spacing=4)

    def run_transform(op: str):
        seq = get_clean_seq(warn=True)
        if not seq:
            return
        if op == "rev":
            tool_out_field.value = seq[::-1]
            lbl = "Reversed (3' -> 5')"
        elif op == "comp":
            tool_out_field.value = complement_only(seq)
            lbl = "Complement (IUPAC)"
        elif op == "revcomp":
            tool_out_field.value = reverse_complement(seq)
            lbl = "Reverse Complement"
        elif op == "rna":
            tool_out_field.value = transcribe_dna_to_rna(seq)
            lbl = "Transcribed mRNA (T -> U)"
        elif op == "cdna":
            tool_out_field.value = back_transcribe_rna_to_dna(seq)
            lbl = "Back-Transcribed cDNA"
        elif op == "fasta":
            tool_out_field.value = format_fasta(state["header"], seq, 70)
            lbl = "Formatted FASTA"
        else:
            tool_out_field.value = "".join(b for b in seq if b in "ATGC")
            lbl = "Stripped Ambiguous Bases"
        notify(f"Applied: {lbl}")

    def run_motif_search(_e=None):
        seq = get_clean_seq(warn=True)
        if not seq:
            return
        pat = (motif_field.value or "").strip()
        if not pat:
            notify("Enter an IUPAC motif first.")
            return
        hits = search_motif(seq, pat, bool(motif_rev_chk.value))
        motif_hits_col.controls = [
            ft.Container(
                content=ft.Text(f"#{i}  [{st}]  {s:,}..{e:,} bp   →   {ms}", font_family="monospace", size=12),
                bgcolor="#ffffff",
                padding=6,
                border_radius=6,
                border=ft.Border.all(1, "#e2e8f0"),
            )
            for i, (st, s, e, ms) in enumerate(hits[:120], 1)
        ]
        notify(f"Found {len(hits)} match(es) for '{pat.upper()}'.")

    def send_tool_to_main(_e):
        txt = (tool_out_field.value or "").strip()
        if not txt:
            return
        recs = parse_fasta_records(txt)
        seq_input.value = recs[0][1] if recs else txt
        update_seq_len()
        notify("Sent transformed sequence to Main Input.")

    view_tools = ft.Column([
        ft.Text("Sequence Transformations", weight=ft.FontWeight.BOLD, size=14),
        ft.Row([
            ft.Button("Reverse", on_click=lambda _e: run_transform("rev")),
            ft.Button("Complement", on_click=lambda _e: run_transform("comp")),
            ft.Button("Rev-Complement", on_click=lambda _e: run_transform("revcomp")),
            ft.Button("DNA → RNA", on_click=lambda _e: run_transform("rna")),
            ft.Button("RNA → DNA", on_click=lambda _e: run_transform("cdna")),
            ft.Button("Format FASTA", on_click=lambda _e: run_transform("fasta")),
            ft.Button("ATGC Only", on_click=lambda _e: run_transform("strip")),
        ], wrap=True),
        tool_out_field,
        ft.Row([
            ft.OutlinedButton("Copy Output", icon=ft.Icons.COPY, on_click=lambda _e: copy_text(tool_out_field.value, "Output")),
            ft.OutlinedButton("Send to Main Input", icon=ft.Icons.UPLOAD, on_click=send_tool_to_main),
        ], wrap=True),
        ft.Divider(),
        ft.Text("IUPAC Motif & Restriction Site Finder", weight=ft.FontWeight.BOLD, size=14),
        ft.Row([
            motif_field,
            ft.Dropdown(
                label="Preset",
                width=185,
                options=[ft.dropdown.Option(k) for k in COMMON_MOTIFS.keys()],
                on_select=lambda e: (setattr(motif_field, "value", COMMON_MOTIFS[e.control.value]), run_motif_search()),
            ),
        ], wrap=True),
        ft.Row([
            motif_rev_chk,
            ft.Button("Search Motif", icon=ft.Icons.SEARCH, on_click=run_motif_search),
        ]),
        motif_hits_col,
    ], spacing=10)

    # --------------------------------------------------------------------------
    # MODULE 3: PROTEIN TRANSLATION & PHYSICOCHEMICAL PROPERTIES
    # --------------------------------------------------------------------------
    gcode_dd = ft.Dropdown(
        label="Genetic Code",
        value="1 - Standard Code",
        width=240,
        options=[ft.dropdown.Option(k) for k in GENETIC_CODES.keys()],
    )
    frame_dd = ft.Dropdown(
        label="Reading Frame",
        value="Frame +1",
        width=165,
        options=[ft.dropdown.Option(x) for x in ["Frame +1", "Frame +2", "Frame +3", "Frame -1", "Frame -2", "Frame -3", "All 6 Frames"]],
    )
    protein_out_field = ft.TextField(
        label="Translated Protein Sequence",
        multiline=True,
        min_lines=4,
        max_lines=6,
        text_style=ft.TextStyle(font_family="monospace", size=12),
        bgcolor="#ffffff",
    )
    protein_props_col = ft.Column(spacing=5)

    def run_translate(_e=None):
        seq = get_clean_seq(warn=True)
        if not seq:
            return
        ctable, _ = GENETIC_CODES[gcode_dd.value]
        fchoice = frame_dd.value
        if fchoice == "All 6 Frames":
            blocks = []
            for f in (1, 2, 3, -1, -2, -3):
                lbl = f"+{f}" if f > 0 else str(f)
                pep = translate_frame(seq, f, ctable, "*")
                blocks.append(f">Frame {lbl} ({len(pep)} aa)\n{pep}")
            protein_out_field.value = "\n\n".join(blocks)
            target_pep = translate_frame(seq, 1, ctable, "*")
            analyzed_lbl = "Frame +1"
        else:
            fnum = int(fchoice.replace("Frame ", ""))
            target_pep = translate_frame(seq, fnum, ctable, "*")
            protein_out_field.value = target_pep
            analyzed_lbl = fchoice

        props = compute_protein_properties(target_pep.replace("*", ""))
        stops = target_pep.count("*")
        items = [
            ("Analyzed Frame", f"{analyzed_lbl} ({gcode_dd.value})"),
            ("Residue Length", f"{props['length']:,} aa (Stop codons: {stops})"),
            ("Molecular Weight", f"{props['mw']:,.2f} Da ({props['mw']/1000:,.2f} kDa)"),
            ("Theoretical pI", f"{props['pi']:.2f}"),
            ("GRAVY Hydropathy", f"{props['gravy']:+.3f} ({'Hydrophobic' if props['gravy'] > 0 else 'Hydrophilic'})"),
            ("Aromaticity (F+W+Y)", f"{props['aromaticity']:.2f}%"),
            ("Ext. Coeff (280 nm)", f"{props['ext_reduced']:,} M⁻¹cm⁻¹ (red) | {props['ext_oxidized']:,} (ox)"),
            ("Charged Residues", f"Acidic (D+E): {props['acidic']} | Basic (K+R+H): {props['basic']}"),
        ]
        protein_props_col.controls = [
            ft.Container(
                content=ft.Text(f"{k}: {v}", size=12, weight=ft.FontWeight.W_600),
                bgcolor="#ffffff",
                padding=7,
                border_radius=6,
                border=ft.Border.all(1, "#e2e8f0"),
            )
            for k, v in items
        ]
        notify(f"Protein translation complete ({analyzed_lbl}).")

    async def export_protein_faa(_e):
        txt = (protein_out_field.value or "").strip()
        if not txt:
            notify("Translate a sequence first.")
            return
        out = txt if txt.startswith(">") else format_fasta(f"{state['header']}_{frame_dd.value}", txt, 70)
        await export_text_or_bytes("translated_protein.faa", out.encode("utf-8"))

    view_protein = ft.Column([
        ft.Row([gcode_dd, frame_dd], wrap=True),
        ft.Row([
            ft.Button("Translate", icon=ft.Icons.TRANSLATE, on_click=run_translate),
            ft.OutlinedButton("Copy Protein", icon=ft.Icons.COPY, on_click=lambda _e: copy_text(protein_out_field.value, "Protein")),
            ft.OutlinedButton("Export .FAA", icon=ft.Icons.DOWNLOAD, on_click=export_protein_faa),
        ], wrap=True),
        protein_out_field,
        ft.Text("Physicochemical Protein Properties", weight=ft.FontWeight.BOLD, size=14),
        protein_props_col,
    ], spacing=10)

    # --------------------------------------------------------------------------
    # MODULE 4: ORF FINDER
    # --------------------------------------------------------------------------
    orf_min_field = ft.TextField(label="Min Length (bp)", value="75", width=130, bgcolor="#ffffff")
    orf_start_dd = ft.Dropdown(
        label="Start Codons",
        value="ATG only",
        width=185,
        options=[ft.dropdown.Option("ATG only"), ft.dropdown.Option("ATG + GTG/TTG")],
    )
    orf_strand_dd = ft.Dropdown(
        label="Strands",
        value="Both",
        width=145,
        options=[ft.dropdown.Option("Both"), ft.dropdown.Option("Forward (+)"), ft.dropdown.Option("Reverse (-)")],
    )
    orf_list_col = ft.Column(spacing=5)
    orf_detail_field = ft.TextField(
        label="Selected ORF Nucleotide & Protein Sequence",
        multiline=True,
        min_lines=4,
        max_lines=6,
        text_style=ft.TextStyle(font_family="monospace", size=11),
        bgcolor="#ffffff",
    )

    def show_orf_detail(orf_obj: dict):
        orf_detail_field.value = (
            f">{orf_obj['id']} | Frame {orf_obj['frame']} | {orf_obj['start']}..{orf_obj['end']} | "
            f"{orf_obj['nt_len']} nt / {orf_obj['aa_len']} aa\n"
            f"[Nucleotide]\n{orf_obj['nt_seq']}\n\n[Protein]\n{orf_obj['aa_seq']}"
        )
        page.update()

    def run_find_orfs(_e=None):
        seq = get_clean_seq(warn=True)
        if not seq:
            return
        ctable, def_starts = GENETIC_CODES[gcode_dd.value]
        starts = {"ATG"} if orf_start_dd.value == "ATG only" else (def_starts | {"ATG", "GTG", "TTG"})
        try:
            min_bp = max(3, int(orf_min_field.value or 75))
        except ValueError:
            min_bp = 75
        orfs = find_orfs(seq, min_bp, ctable, starts, orf_strand_dd.value)
        state["detected_orfs"] = orfs
        orf_list_col.controls = [
            ft.Container(
                content=ft.Row([
                    ft.Text(f"{o['id']} (Frame {o['frame']})", weight=ft.FontWeight.BOLD, size=12, expand=True),
                    ft.Text(f"{o['start']}..{o['end']} | {o['nt_len']} nt ({o['aa_len']} aa)", size=11, color="#334155"),
                ]),
                bgcolor="#ffffff",
                padding=9,
                border_radius=6,
                border=ft.Border.all(1, "#cbd5e1"),
                on_click=lambda _e, item=o: show_orf_detail(item),
            )
            for o in orfs[:80]
        ]
        if orfs:
            show_orf_detail(orfs[0])
        else:
            orf_detail_field.value = f"No ORFs >= {min_bp} bp found."
        notify(f"ORF scan complete: {len(orfs)} ORF(s) found.")

    async def export_orfs_fasta(_e):
        if not state["detected_orfs"]:
            notify("Scan ORFs first.")
            return
        out = "".join(
            format_fasta(f"{o['id']} Frame={o['frame']} {o['start']}..{o['end']}", o["nt_seq"], 70)
            for o in state["detected_orfs"]
        )
        await export_text_or_bytes("detected_orfs.fasta", out.encode("utf-8"))

    view_orf = ft.Column([
        ft.Row([orf_min_field, orf_start_dd, orf_strand_dd], wrap=True),
        ft.Row([
            ft.Button("Scan ORFs", icon=ft.Icons.FIND_IN_PAGE, on_click=run_find_orfs),
            ft.OutlinedButton("Export ORFs (.FASTA)", icon=ft.Icons.DOWNLOAD, on_click=export_orfs_fasta),
        ], wrap=True),
        ft.Text("Detected ORFs (Tap any card to inspect)", weight=ft.FontWeight.BOLD, size=13),
        orf_list_col,
        orf_detail_field,
        ft.OutlinedButton("Copy Selected ORF", icon=ft.Icons.COPY, on_click=lambda _e: copy_text(orf_detail_field.value, "ORF")),
    ], spacing=10)

    # --------------------------------------------------------------------------
    # MODULE 5: FULL-SCREEN VISUALIZATION & GRAPH STUDIO (BAR, PIE, COLORS)
    # --------------------------------------------------------------------------
    viz_mode_dd = ft.Dropdown(
        label="Chart / Graph Type",
        value="Nucleotide Bar Chart (A, T, G, C)",
        width=290,
        options=[
            ft.dropdown.Option("Nucleotide Bar Chart (A, T, G, C)"),
            ft.dropdown.Option("Nucleotide Pie Chart (Base %)"),
            ft.dropdown.Option("Bar + Pie Combined Dashboard"),
            ft.dropdown.Option("Purine/Pyrimidine & GC/AT (Bar + Pie)"),
            ft.dropdown.Option("Sliding Window GC Content (%)"),
            ft.dropdown.Option("GC & AT Skew Profile"),
            ft.dropdown.Option("Cumulative GC Skew (OriC Finder)"),
            ft.dropdown.Option("Dinucleotide Frequencies (Bar Chart)"),
            ft.dropdown.Option("Amino Acid Class (Bar + Pie)"),
        ],
    )
    viz_palette_dd = ft.Dropdown(
        label="Graph Color Theme",
        value="Classic Bioinformatics",
        width=240,
        options=[ft.dropdown.Option(k) for k in COLOR_PALETTES.keys()],
    )
    viz_color_target_dd = ft.Dropdown(
        label="Target to Recolor",
        value="Primary / All Bars",
        width=175,
        options=[
            ft.dropdown.Option("Primary / All Bars"),
            ft.dropdown.Option("Base A"),
            ft.dropdown.Option("Base T"),
            ft.dropdown.Option("Base G"),
            ft.dropdown.Option("Base C"),
        ],
    )
    viz_title_field = ft.TextField(label="Custom Chart Title (Optional)", width=270, bgcolor="#ffffff")
    viz_win_field = ft.TextField(label="Window (bp)", value="50", width=105, bgcolor="#ffffff")
    viz_step_field = ft.TextField(label="Step (bp)", value="5", width=95, bgcolor="#ffffff")
    viz_labels_chk = ft.Checkbox(label="Show Values & %", value=True)
    viz_grid_chk = ft.Checkbox(label="Show Gridlines", value=True)
    viz_badge_txt = ft.Text("", size=12, weight=ft.FontWeight.BOLD, color="#1e40af")
    viz_image = ft.Image(src="", fit=ft.BoxFit.CONTAIN, width=720)

    def redraw_studio_chart(_e=None):
        seq = get_clean_seq(warn=False)
        if not seq:
            return
        try:
            w_val = max(5, int(viz_win_field.value or 50))
        except ValueError:
            w_val = 50
        try:
            s_val = max(1, int(viz_step_field.value or 5))
        except ValueError:
            s_val = 5
        svg = render_viz_svg(
            seq=seq,
            header=state["header"],
            mode=viz_mode_dd.value,
            base_colors=state["base_colors"],
            primary_color=state["primary_color"],
            custom_title=viz_title_field.value or "",
            win_size=w_val,
            step_size=s_val,
            show_labels=bool(viz_labels_chk.value),
            show_grid=bool(viz_grid_chk.value),
            gcode_key=gcode_dd.value,
        )
        state["last_viz_svg"] = svg
        viz_image.src = svg_to_data_uri(svg)
        page.update()

    def on_palette_select(e):
        pal = COLOR_PALETTES.get(e.control.value)
        if pal:
            state["base_colors"] = {"A": pal["A"], "T": pal["T"], "G": pal["G"], "C": pal["C"]}
            state["primary_color"] = pal["primary"]
        else:
            for b in ("A", "T", "G", "C"):
                state["base_colors"][b] = state["primary_color"]
        redraw_studio_chart()

    def apply_swatch_color(hex_col: str):
        target = viz_color_target_dd.value
        if target == "Primary / All Bars":
            state["primary_color"] = hex_col
            viz_palette_dd.value = "Uniform Single Color"
            for b in ("A", "T", "G", "C"):
                state["base_colors"][b] = hex_col
        elif target.startswith("Base "):
            base = target.replace("Base ", "")
            state["base_colors"][base] = hex_col
        redraw_studio_chart()

    viz_mode_dd.on_select = redraw_studio_chart
    viz_palette_dd.on_select = on_palette_select
    viz_labels_chk.on_change = redraw_studio_chart
    viz_grid_chk.on_change = redraw_studio_chart

    async def export_chart_svg(_e):
        if not state["last_viz_svg"]:
            notify("Generate a chart first.")
            return
        await export_text_or_bytes("dna_analyzer_chart.svg", state["last_viz_svg"].encode("utf-8"))

    def open_viz_studio(_e=None):
        seq = get_clean_seq(warn=True)
        if not seq:
            return
        st = compute_dna_statistics(seq)
        viz_badge_txt.value = f"{state['header'][:28]}  |  {st['length']:,} bp  |  GC: {st['gc_strict']:.2f}%"
        state["in_viz_studio"] = True
        main_shell.visible = False
        viz_studio_shell.visible = True
        redraw_studio_chart()
        notify(f"Graph Studio opened — {viz_mode_dd.value}")

    def close_viz_studio(_e=None):
        state["in_viz_studio"] = False
        viz_studio_shell.visible = False
        main_shell.visible = True
        page.update()

    swatch_row = ft.Row(
        [
            ft.Container(
                width=28,
                height=28,
                bgcolor=c,
                border_radius=14,
                border=ft.Border.all(1.5, "#0f172a"),
                on_click=lambda _e, col=c: apply_swatch_color(col),
            )
            for c in SWATCH_COLORS
        ],
        wrap=True,
        spacing=6,
    )

    viz_studio_shell = ft.Container(
        visible=False,
        expand=True,
        bgcolor="#f1f5f9",
        padding=10,
        content=ft.Column(
            scroll=ft.ScrollMode.AUTO,
            spacing=10,
            controls=[
                ft.Container(
                    bgcolor="#0f172a",
                    padding=10,
                    border_radius=10,
                    content=ft.Row([
                        ft.Button(
                            "← Back",
                            bgcolor="#2563eb",
                            color="#ffffff",
                            on_click=close_viz_studio,
                        ),
                        ft.Text("GRAPH & VISUALIZATION STUDIO", color="#ffffff", weight=ft.FontWeight.BOLD, size=14, expand=True),
                        ft.IconButton(icon=ft.Icons.DOWNLOAD, icon_color="#ffffff", tooltip="Save Chart (.SVG)", on_click=export_chart_svg),
                    ]),
                ),
                viz_badge_txt,
                # Full-height spacious chart container (sequence input & other panels hidden!)
                ft.Container(
                    content=viz_image,
                    bgcolor="#ffffff",
                    padding=6,
                    border_radius=10,
                    border=ft.Border.all(1, "#cbd5e1"),
                    alignment=ft.Alignment(0, 0),
                ),
                # Graph & Color Editing Panel below chart
                ft.Container(
                    bgcolor="#ffffff",
                    padding=12,
                    border_radius=10,
                    border=ft.Border.all(1, "#cbd5e1"),
                    content=ft.Column([
                        ft.Text("Graph Type & Color Customization", weight=ft.FontWeight.BOLD, size=14),
                        ft.Row([viz_mode_dd, viz_palette_dd], wrap=True),
                        ft.Text("Tap a Color Swatch for Selected Target:", size=12, weight=ft.FontWeight.W_600),
                        ft.Row([viz_color_target_dd, swatch_row], wrap=True),
                        ft.Row([viz_title_field, viz_win_field, viz_step_field], wrap=True),
                        ft.Row([viz_labels_chk, viz_grid_chk], wrap=True),
                        ft.Row([
                            ft.Button("Update Graph", icon=ft.Icons.REFRESH, on_click=redraw_studio_chart),
                            ft.OutlinedButton("Save Chart (.SVG)", icon=ft.Icons.DOWNLOAD, on_click=export_chart_svg),
                        ], wrap=True),
                    ], spacing=8),
                ),
            ],
        ),
    )

    # --------------------------------------------------------------------------
    # MODULE 6: PAIRWISE ALIGNMENT & DOT PLOT
    # --------------------------------------------------------------------------
    seq_b_field = ft.TextField(
        label="Reference Sequence B (Compared against Main Input)",
        multiline=True,
        min_lines=3,
        max_lines=5,
        text_style=ft.TextStyle(font_family="monospace", size=12),
        bgcolor="#ffffff",
    )
    aln_mode_dd = ft.Dropdown(
        label="Alignment Mode",
        value="Global (Needleman-Wunsch)",
        width=230,
        options=[ft.dropdown.Option("Global (Needleman-Wunsch)"), ft.dropdown.Option("Local (Smith-Waterman)")],
    )
    aln_summary_txt = ft.Text("Enter Sequence B and tap Run Alignment.", size=12, weight=ft.FontWeight.BOLD, color="#1e40af")
    aln_out_field = ft.TextField(
        label="Pairwise Alignment Output",
        multiline=True,
        min_lines=5,
        max_lines=8,
        text_style=ft.TextStyle(font_family="monospace", size=11),
        bgcolor="#ffffff",
    )
    dotplot_img = ft.Image(src="", fit=ft.BoxFit.CONTAIN, width=420, visible=False)

    def run_alignment(_e=None):
        seq_a = get_clean_seq(warn=True)
        if not seq_a:
            return
        recs_b = parse_fasta_records(seq_b_field.value or "")
        seq_b = recs_b[0][1] if recs_b else ""
        if not seq_b:
            notify("Enter Reference Sequence B first.")
            return
        sa, sb = seq_a[:650], seq_b[:650]
        mode = "global" if "Global" in aln_mode_dd.value else "local"
        res = pairwise_align(sa, sb, mode=mode)
        if not res:
            return
        aln_summary_txt.value = (
            f"Score: {res['score']} | Identity: {res['identity']:.2f}% ({res['matches']}/{res['aligned_len']}) | "
            f"Mismatches: {res['mismatches']} | Gaps: {res['gaps']}"
        )
        blocks = []
        for i in range(0, res["aligned_len"], 45):
            blocks.append(
                f"A: {res['aln_a'][i:i+45]}\n"
                f"   {res['midline'][i:i+45]}\n"
                f"B: {res['aln_b'][i:i+45]}"
            )
        aln_out_field.value = "\n\n".join(blocks)
        notify(f"Alignment complete ({res['identity']:.1f}% identity).")

    def run_dotplot(_e=None):
        seq_a = get_clean_seq(warn=True)
        if not seq_a:
            return
        recs_b = parse_fasta_records(seq_b_field.value or "")
        seq_b = recs_b[0][1] if recs_b else seq_a
        dotplot_img.src = svg_to_data_uri(render_dotplot_svg(seq_a, seq_b))
        dotplot_img.visible = True
        notify("Generated sequence dot plot.")

    view_compare = ft.Column([
        ft.Row([
            ft.OutlinedButton("Use RevComp of Main", on_click=lambda _e: (setattr(seq_b_field, "value", reverse_complement(get_clean_seq(True))), page.update())),
            aln_mode_dd,
        ], wrap=True),
        seq_b_field,
        ft.Row([
            ft.Button("Run Alignment", icon=ft.Icons.COMPARE_ARROWS, on_click=run_alignment),
            ft.OutlinedButton("Generate Dot Plot", icon=ft.Icons.SCATTER_PLOT, on_click=run_dotplot),
        ], wrap=True),
        aln_summary_txt,
        aln_out_field,
        dotplot_img,
    ], spacing=10)

    # --------------------------------------------------------------------------
    # MODULE 7: PHYLOGENETIC TREE CONSTRUCTION
    # --------------------------------------------------------------------------
    tree_method_dd = ft.Dropdown(
        label="Tree Method",
        value="Neighbor-Joining (NJ)",
        width=225,
        options=[
            ft.dropdown.Option("Neighbor-Joining (NJ)"),
            ft.dropdown.Option("UPGMA (Average Linkage)"),
            ft.dropdown.Option("WPGMA (Weighted Linkage)"),
            ft.dropdown.Option("Single-Linkage (Minimum)"),
            ft.dropdown.Option("Complete-Linkage (Maximum)"),
            ft.dropdown.Option("Maximum Parsimony (Fitch + NNI)"),
        ],
    )
    tree_model_dd = ft.Dropdown(
        label="Distance Model",
        value="Kimura 2-Parameter (K2P)",
        width=215,
        options=[
            ft.dropdown.Option("Kimura 2-Parameter (K2P)"),
            ft.dropdown.Option("Jukes-Cantor (JC69)"),
            ft.dropdown.Option("p-distance (Uncorrected)"),
        ],
    )
    tree_style_dd = ft.Dropdown(
        label="Tree Style",
        value="Rectangular Phylogram",
        width=215,
        options=[
            ft.dropdown.Option("Rectangular Phylogram"),
            ft.dropdown.Option("Rectangular Cladogram"),
            ft.dropdown.Option("Circular / Radial Phylogram"),
            ft.dropdown.Option("Unrooted Radial Layout"),
        ],
    )
    tree_len_chk = ft.Checkbox(label="Show Branch Lengths", value=True)
    phylo_input = ft.TextField(
        label="Multi-FASTA Taxa Input (>= 3 Sequences)",
        multiline=True,
        min_lines=5,
        max_lines=7,
        text_style=ft.TextStyle(font_family="monospace", size=11),
        bgcolor="#ffffff",
    )
    phylo_img = ft.Image(src="", fit=ft.BoxFit.CONTAIN, width=720, visible=False)
    phylo_meta_field = ft.TextField(
        label="Newick String & Evolutionary Distance Matrix",
        multiline=True,
        min_lines=4,
        max_lines=6,
        text_style=ft.TextStyle(font_family="monospace", size=11),
        bgcolor="#ffffff",
    )

    def redraw_phylo(_e=None):
        if state["last_phylo_tree"] is None:
            return
        svg = render_phylo_svg(state["last_phylo_tree"], tree_method_dd.value, tree_style_dd.value, bool(tree_len_chk.value))
        state["last_phylo_svg"] = svg
        phylo_img.src = svg_to_data_uri(svg)
        phylo_img.visible = True
        page.update()

    tree_style_dd.on_select = redraw_phylo
    tree_len_chk.on_change = redraw_phylo

    def run_build_tree(_e=None):
        recs = parse_fasta_records(phylo_input.value or "", keep_gaps=True)
        if len(recs) < 3:
            notify("Provide >= 3 FASTA sequences or tap 'Sample Taxa'.")
            return
        recs = recs[:40]
        aligned = progressive_msa(recs)
        names, dmat = build_distance_matrix(aligned, tree_model_dd.value)
        method = tree_method_dd.value
        extra = ""
        if "Neighbor-Joining" in method:
            tree = build_nj_tree(names, dmat)
        elif "UPGMA" in method:
            tree = build_agglomerative_tree(names, dmat, "UPGMA")
        elif "WPGMA" in method:
            tree = build_agglomerative_tree(names, dmat, "WPGMA")
        elif "Single-Linkage" in method:
            tree = build_agglomerative_tree(names, dmat, "Single-Linkage")
        elif "Complete-Linkage" in method:
            tree = build_agglomerative_tree(names, dmat, "Complete-Linkage")
        else:
            tree, steps = build_maximum_parsimony_tree(names, aligned, dmat)
            extra = f" | Parsimony: {steps} steps"

        state["last_phylo_tree"] = tree
        state["last_newick"] = tree.to_newick(True)
        lines = [
            f"Method: {method}{extra} | Model: {tree_model_dd.value}",
            f"Newick: {state['last_newick']}",
            "-" * 48,
            "Pairwise Distance Matrix:",
        ]
        for i, nm in enumerate(names):
            lines.append(f"{nm[:14].ljust(15)} " + " ".join(f"{dmat[i][j]:.4f}" for j in range(len(names))))
        phylo_meta_field.value = "\n".join(lines)
        redraw_phylo()
        notify(f"Phylogenetic tree built ({len(names)} taxa){extra}.")

    async def export_newick(_e):
        if not state["last_newick"]:
            notify("Construct a tree first.")
            return
        await export_text_or_bytes("phylogenetic_tree.nwk", (state["last_newick"] + "\n").encode("utf-8"))

    async def export_tree_svg(_e):
        if not state["last_phylo_svg"]:
            notify("Construct a tree first.")
            return
        await export_text_or_bytes("phylogenetic_tree.svg", state["last_phylo_svg"].encode("utf-8"))

    view_phylo = ft.Column([
        ft.Row([
            ft.OutlinedButton("Upload Multi-FASTA", icon=ft.Icons.FOLDER_OPEN, on_click=pick_phylo_fasta),
            ft.OutlinedButton("Sample Taxa", on_click=lambda _e: (setattr(phylo_input, "value", SAMPLE_PHYLO_MULTI_FASTA.strip()), run_build_tree())),
        ], wrap=True),
        phylo_input,
        ft.Row([tree_method_dd, tree_model_dd, tree_style_dd, tree_len_chk], wrap=True),
        ft.Row([
            ft.Button("Construct Tree", icon=ft.Icons.ACCOUNT_TREE, on_click=run_build_tree),
            ft.OutlinedButton("Copy Newick", icon=ft.Icons.COPY, on_click=lambda _e: copy_text(state["last_newick"], "Newick")),
            ft.OutlinedButton("Export .NWK", on_click=export_newick),
            ft.OutlinedButton("Save Tree (.SVG)", on_click=export_tree_svg),
        ], wrap=True),
        ft.Container(content=phylo_img, bgcolor="#ffffff", padding=6, border_radius=8, border=ft.Border.all(1, "#cbd5e1")),
        phylo_meta_field,
    ], spacing=10)

    # --------------------------------------------------------------------------
    # MODULE SWITCHER & MAIN MOBILE SHELL
    # --------------------------------------------------------------------------
    module_views = [view_stats, view_tools, view_protein, view_orf, view_compare, view_phylo]
    module_names = ["1. DNA Stats", "2. Tools", "3. Protein", "4. ORFs", "5. Graphs", "6. Align", "7. Phylo Tree"]
    module_host = ft.Container(content=view_stats, padding=ft.Padding(0, 6, 0, 20))
    tab_buttons = []

    def refresh_active_module():
        idx = state["active_tab"]
        if idx == 0 and get_clean_seq(False):
            run_analyze_dna()
        elif idx == 2 and get_clean_seq(False):
            run_translate()
        elif idx == 3 and get_clean_seq(False):
            run_find_orfs()

    def select_tab(idx: int):
        if idx == 4:  # "5. Graphs" -> Opens dedicated Full-Screen Graph Studio
            open_viz_studio()
            return
        state["active_tab"] = idx
        view_map = {0: view_stats, 1: view_tools, 2: view_protein, 3: view_orf, 5: view_compare, 6: view_phylo}
        module_host.content = view_map[idx]
        for b_idx, btn in enumerate(tab_buttons):
            active = b_idx == idx
            btn.bgcolor = "#1e40af" if active else "#e2e8f0"
            btn.content.color = "#ffffff" if active else "#0f172a"
        refresh_active_module()
        page.update()

    for i, label in enumerate(module_names):
        active = i == 0
        btn = ft.Container(
            content=ft.Text(label, size=12, weight=ft.FontWeight.BOLD, color="#ffffff" if active else "#0f172a"),
            bgcolor="#1e40af" if active else "#e2e8f0",
            padding=ft.Padding(12, 8, 12, 8),
            border_radius=18,
            on_click=lambda _e, tidx=i: select_tab(tidx),
        )
        tab_buttons.append(btn)

    async def export_full_txt_report(_e):
        seq = get_clean_seq(warn=True)
        if not seq:
            return
        st = compute_dna_statistics(seq)
        report = (
            f"DNA ANALYZER V10 - REPORT\n"
            f"Header: {state['header']}\n"
            f"Length: {st['length']:,} bp | GC: {st['gc_strict']:.2f}% | Tm: {st['tm_basic']:.2f} C\n\n"
            f"{format_fasta(state['header'], seq, 70)}"
        )
        await export_text_or_bytes("DNA_Analyzer_Report.txt", report.encode("utf-8"))

    main_shell = ft.Container(
        visible=True,
        expand=True,
        padding=12,
        content=ft.Column(
            scroll=ft.ScrollMode.AUTO,
            spacing=10,
            controls=[
                # App Header Banner
                ft.Container(
                    bgcolor="#0f172a",
                    padding=12,
                    border_radius=10,
                    content=ft.Row([
                        ft.Column([
                            ft.Text("DNA ANALYZER V10", color="#ffffff", size=17, weight=ft.FontWeight.BOLD),
                            ft.Text("Bioinformatics & Phylogenetic Toolkit", color="#94a3b8", size=11),
                        ], spacing=1, expand=True),
                        len_badge,
                    ]),
                ),
                # Sequence Input Card (Starts Empty)
                ft.Container(
                    bgcolor="#ffffff",
                    padding=10,
                    border_radius=10,
                    border=ft.Border.all(1, "#cbd5e1"),
                    content=ft.Column([
                        seq_input,
                        ft.Row([
                            ft.Button("Analyze", icon=ft.Icons.PLAY_ARROW, on_click=lambda _e: select_tab(0)),
                            ft.Button("Graphs", icon=ft.Icons.BAR_CHART, bgcolor="#0d9488", color="#ffffff", on_click=open_viz_studio),
                            ft.OutlinedButton("Upload FASTA", icon=ft.Icons.UPLOAD_FILE, on_click=pick_main_fasta),
                            ft.OutlinedButton("Sample", on_click=load_sample_dna),
                            ft.OutlinedButton("Export Report", on_click=export_full_txt_report),
                            ft.OutlinedButton("Clear", on_click=clear_all),
                        ], wrap=True, spacing=6),
                    ], spacing=8),
                ),
                # Scrollable Module Selector Pills (All 7 Modules)
                ft.Row(controls=tab_buttons, scroll=ft.ScrollMode.AUTO, spacing=6),
                # Active Module Content
                module_host,
                # Status Footer
                ft.Container(
                    content=status_txt,
                    bgcolor="#e2e8f0",
                    padding=8,
                    border_radius=6,
                ),
            ],
        ),
    )

    page.add(ft.Stack([main_shell, viz_studio_shell], expand=True))


if __name__ == "__main__":
    ft.run(main)