import argparse
import math
from pathlib import Path

import pandas as pd
from Bio import SeqIO, pairwise2
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.Polypeptide import is_aa


AA3_TO_1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Compute AlphaFold-derived per-position structural features and merge them into "
            "the PRNP mutation feature table."
        )
    )
    parser.add_argument(
        "--af_structure",
        required=True,
        help="Path to AlphaFold structure file (.pdb or .cif/.mmcif).",
    )
    parser.add_argument(
        "--reference_fasta",
        default="prnp_reference_protein.fasta",
        help="Reference protein FASTA used by mutation generation.",
    )
    parser.add_argument(
        "--tier2_csv",
        default="prnp_mutation_features_tier2.csv",
        help="Tier-2 mutation feature CSV to enrich.",
    )
    parser.add_argument(
        "--output_csv",
        default="prnp_mutation_features_tier2_structural.csv",
        help="Output CSV path with structural features merged by position.",
    )
    parser.add_argument(
        "--position_csv",
        default="prnp_alphafold_position_features.csv",
        help="Optional output CSV for per-position structural features.",
    )
    parser.add_argument(
        "--chain",
        default=None,
        help="Optional chain ID. If omitted, first amino-acid chain is used.",
    )
    return parser.parse_args()


def get_parser(path):
    suffix = path.suffix.lower()
    if suffix in {".cif", ".mmcif"}:
        return MMCIFParser(QUIET=True)
    return PDBParser(QUIET=True)


def select_chain(structure, requested_chain=None):
    model = next(structure.get_models())

    if requested_chain is not None:
        if requested_chain not in model:
            raise ValueError(f"Requested chain not present in structure: {requested_chain}")
        return model[requested_chain]

    for chain in model:
        for residue in chain.get_residues():
            if is_aa(residue, standard=True):
                return chain

    raise ValueError("No amino-acid chain found in first model.")


def extract_chain_residues(chain):
    residues = []

    for residue in chain.get_residues():
        if not is_aa(residue, standard=True):
            continue

        aa3 = residue.get_resname().upper()
        aa1 = AA3_TO_1.get(aa3)
        if aa1 is None:
            continue

        seq_id = residue.id[1]
        ins_code = residue.id[2].strip() if isinstance(residue.id[2], str) else ""

        ca = residue["CA"] if "CA" in residue else None
        if ca is None:
            coord = (math.nan, math.nan, math.nan)
            plddt = math.nan
        else:
            xyz = ca.get_coord()
            coord = (float(xyz[0]), float(xyz[1]), float(xyz[2]))
            # AlphaFold PDB stores pLDDT in B-factor.
            plddt = float(ca.get_bfactor())

        residues.append(
            {
                "aa": aa1,
                "seq_id": int(seq_id),
                "ins_code": ins_code,
                "coord": coord,
                "plddt": plddt,
            }
        )

    if not residues:
        raise ValueError("No standard amino-acid residues found in selected chain.")

    return residues


def align_reference_to_structure(reference_seq, structure_seq):
    alignment = pairwise2.align.globalms(
        reference_seq,
        structure_seq,
        2,
        -1,
        -10,
        -0.5,
        one_alignment_only=True,
    )
    if not alignment:
        raise ValueError("Could not align reference sequence to AlphaFold structure sequence.")

    aligned_ref, aligned_struct = alignment[0].seqA, alignment[0].seqB
    ref_to_struct = {}

    ref_pos = 0
    struct_pos = 0
    matched_aa = 0

    for ref_char, struct_char in zip(aligned_ref, aligned_struct):
        if ref_char != "-":
            ref_pos += 1
        if struct_char != "-":
            struct_pos += 1

        if ref_char != "-" and struct_char != "-":
            ref_to_struct[ref_pos] = struct_pos
            if ref_char == struct_char:
                matched_aa += 1

    identity = matched_aa / max(1, len(reference_seq))
    return ref_to_struct, identity


def finite_coord(coord):
    return all(not math.isnan(v) for v in coord)


def distance(c1, c2):
    return math.sqrt((c1[0] - c2[0]) ** 2 + (c1[1] - c2[1]) ** 2 + (c1[2] - c2[2]) ** 2)


def nanmean(values):
    finite = [v for v in values if not math.isnan(v)]
    if not finite:
        return math.nan
    return sum(finite) / len(finite)


def nanmin(values):
    finite = [v for v in values if not math.isnan(v)]
    if not finite:
        return math.nan
    return min(finite)


def compute_position_features(reference_seq, residues, ref_to_struct):
    features = []

    for ref_pos in range(1, len(reference_seq) + 1):
        row = {
            "position": ref_pos,
            "af2_mapped": 0,
            "af2_ref_aa": reference_seq[ref_pos - 1],
            "af2_struct_aa": None,
            "af2_struct_seq_id": math.nan,
            "af2_plddt": math.nan,
            "af2_plddt_window_mean_5": math.nan,
            "af2_plddt_window_min_5": math.nan,
            "af2_neighbor_count_8A": math.nan,
            "af2_neighbor_count_12A": math.nan,
            "af2_mean_neighbor_distance_8A": math.nan,
            "af2_min_neighbor_distance": math.nan,
            "af2_is_buried_proxy": math.nan,
            "af2_is_exposed_proxy": math.nan,
            "af2_is_low_confidence": math.nan,
            "af2_is_high_confidence": math.nan,
        }

        struct_pos = ref_to_struct.get(ref_pos)
        if struct_pos is None:
            features.append(row)
            continue

        residue = residues[struct_pos - 1]
        row["af2_mapped"] = 1
        row["af2_struct_aa"] = residue["aa"]
        row["af2_struct_seq_id"] = residue["seq_id"]
        row["af2_plddt"] = residue["plddt"]

        # Local confidence context around the mapped residue in structure space.
        left = max(1, struct_pos - 2)
        right = min(len(residues), struct_pos + 2)
        local_plddt = [residues[i - 1]["plddt"] for i in range(left, right + 1)]
        row["af2_plddt_window_mean_5"] = nanmean(local_plddt)
        row["af2_plddt_window_min_5"] = nanmin(local_plddt)

        base_coord = residue["coord"]
        if not finite_coord(base_coord):
            row["af2_is_low_confidence"] = 1 if (not math.isnan(row["af2_plddt"]) and row["af2_plddt"] < 70) else 0
            row["af2_is_high_confidence"] = 1 if (not math.isnan(row["af2_plddt"]) and row["af2_plddt"] >= 90) else 0
            features.append(row)
            continue

        dists = []
        for idx, other in enumerate(residues, start=1):
            if idx == struct_pos:
                continue
            if not finite_coord(other["coord"]):
                continue
            dists.append(distance(base_coord, other["coord"]))

        within_8 = [d for d in dists if d <= 8.0]
        within_12 = [d for d in dists if d <= 12.0]

        row["af2_neighbor_count_8A"] = len(within_8)
        row["af2_neighbor_count_12A"] = len(within_12)
        row["af2_mean_neighbor_distance_8A"] = nanmean(within_8)
        row["af2_min_neighbor_distance"] = min(dists) if dists else math.nan

        # Contact-density proxies for buried/exposed status.
        row["af2_is_buried_proxy"] = 1 if len(within_8) >= 14 else 0
        row["af2_is_exposed_proxy"] = 1 if len(within_8) <= 6 else 0

        plddt = row["af2_plddt"]
        row["af2_is_low_confidence"] = 1 if (not math.isnan(plddt) and plddt < 70) else 0
        row["af2_is_high_confidence"] = 1 if (not math.isnan(plddt) and plddt >= 90) else 0

        features.append(row)

    return pd.DataFrame(features)


def main():
    args = parse_args()

    af_path = Path(args.af_structure)
    ref_path = Path(args.reference_fasta)
    tier2_path = Path(args.tier2_csv)

    if not af_path.exists():
        raise FileNotFoundError(f"AlphaFold structure file not found: {af_path}")
    if not ref_path.exists():
        raise FileNotFoundError(f"Reference FASTA not found: {ref_path}")
    if not tier2_path.exists():
        raise FileNotFoundError(f"Tier-2 CSV not found: {tier2_path}")

    reference_seq = str(next(SeqIO.parse(str(ref_path), "fasta")).seq)

    parser = get_parser(af_path)
    structure = parser.get_structure("af", str(af_path))
    chain = select_chain(structure, requested_chain=args.chain)
    residues = extract_chain_residues(chain)
    structure_seq = "".join(r["aa"] for r in residues)

    ref_to_struct, seq_identity = align_reference_to_structure(reference_seq, structure_seq)
    position_df = compute_position_features(reference_seq, residues, ref_to_struct)

    tier2_df = pd.read_csv(tier2_path)
    if "position" not in tier2_df.columns:
        raise ValueError("Tier-2 CSV does not contain required column: position")

    merged = tier2_df.merge(position_df, on="position", how="left")
    merged.to_csv(args.output_csv, index=False)
    position_df.to_csv(args.position_csv, index=False)

    mapped = int(position_df["af2_mapped"].sum())
    low_conf = int(position_df["af2_is_low_confidence"].fillna(0).sum())

    print("Reference length:", len(reference_seq))
    print("Structure residues used:", len(residues))
    print("Reference-to-structure AA identity:", round(seq_identity, 4))
    print("Mapped positions:", mapped)
    print("Low-confidence positions (pLDDT < 70):", low_conf)
    print("Tier-2 rows:", len(tier2_df))
    print("Merged rows:", len(merged))
    print("Merged columns:", merged.shape[1])
    print("Saved merged CSV:", args.output_csv)
    print("Saved per-position CSV:", args.position_csv)


if __name__ == "__main__":
    main()
