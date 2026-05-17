import argparse
import csv
import json
import shutil
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd


DEFAULT_FOLDX = r"D:\FoldX\foldx_20261231.exe"
DEFAULT_ROTABASE = r"D:\FoldX\rotabase.txt"
DEFAULT_PDB = "AF-P04156-F1-model_v4.pdb"
DEFAULT_MUTATIONS = "ddg_inputs/ddg_label_template.csv"
DEFAULT_RUN_DIR = "foldx_runs/prnp_ddg"
DEFAULT_RESULTS = "foldx_ddg_results.csv"
DEFAULT_MERGED = "prnp_mutation_features_ddg.csv"


ENERGY_COLUMNS = [
    "total energy",
    "Backbone Hbond",
    "Sidechain Hbond",
    "Van der Waals",
    "Electrostatics",
    "Solvation Polar",
    "Solvation Hydrophobic",
    "Van der Waals clashes",
    "entropy sidechain",
    "entropy mainchain",
    "sloop_entropy",
    "mloop_entropy",
    "cis_bond",
    "torsional clash",
    "backbone clash",
    "helix dipole",
    "water bridge",
    "disulfide",
    "electrostatic kon",
    "partial covalent bonds",
    "energy Ionisation",
    "Entropy Complex",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Run FoldX BuildModel in resumable batches for all PRNP single-point mutations, "
            "then merge DDG labels into the mutation feature table."
        )
    )
    parser.add_argument("--foldx", default=DEFAULT_FOLDX)
    parser.add_argument("--rotabase", default=DEFAULT_ROTABASE)
    parser.add_argument("--pdb", default=DEFAULT_PDB)
    parser.add_argument("--mutation_csv", default=DEFAULT_MUTATIONS)
    parser.add_argument("--run_dir", default=DEFAULT_RUN_DIR)
    parser.add_argument("--results_csv", default=DEFAULT_RESULTS)
    parser.add_argument("--merged_output_csv", default=DEFAULT_MERGED)
    parser.add_argument("--source_feature_csv", default="prnp_mutation_features_tier2_structural.csv")
    parser.add_argument("--batch_size", type=int, default=50)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--repair_timeout", type=int, default=900)
    parser.add_argument("--batch_timeout", type=int, default=900)
    parser.add_argument("--limit", type=int, default=None, help="Optional mutation limit for testing.")
    parser.add_argument("--force_repair", action="store_true")
    return parser.parse_args()


def ensure_exists(path, label):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"{label} not found: {path}")
    return path


def foldx_mutation_from_row(row):
    if "foldx_mutation" in row and pd.notna(row["foldx_mutation"]):
        return str(row["foldx_mutation"]).strip()
    foldx_position = row.get("foldx_position", row.get("af2_struct_seq_id", row["position"]))
    return f"{row['wt_aa']}A{int(float(foldx_position))}{row['mut_aa']};"


def load_mutations(path, limit=None):
    df = pd.read_csv(path)
    required = {"mutation", "position", "wt_aa", "mut_aa"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"Mutation CSV missing required columns: {missing}")
    df = df.copy()
    df["foldx_mutation"] = df.apply(foldx_mutation_from_row, axis=1)
    if limit is not None:
        df = df.head(limit)
    if "foldx_position" not in df.columns:
        if "af2_struct_seq_id" in df.columns:
            df["foldx_position"] = pd.to_numeric(df["af2_struct_seq_id"], errors="coerce").fillna(df["position"]).astype(int)
        else:
            df["foldx_position"] = df["position"].astype(int)
    return df[["mutation", "position", "foldx_position", "wt_aa", "mut_aa", "foldx_mutation"]]


def load_completed(results_csv):
    path = Path(results_csv)
    if not path.exists():
        return set()
    try:
        df = pd.read_csv(path)
    except pd.errors.EmptyDataError:
        return set()
    if "mutation" not in df.columns or "status" not in df.columns:
        return set()
    return set(df.loc[df["status"] == "ok", "mutation"].astype(str))


def append_rows(path, rows):
    if not rows:
        return
    path = Path(path)
    fieldnames = [
        "mutation",
        "foldx_mutation",
        "position",
        "foldx_position",
        "wt_aa",
        "mut_aa",
        "ddg",
        "status",
        "batch_id",
        "foldx_pdb_row",
        "elapsed_seconds",
        "error",
    ] + [f"foldx_{col.replace(' ', '_')}" for col in ENERGY_COLUMNS if col != "total energy"]

    exists = path.exists() and path.stat().st_size > 0
    with path.open("a", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        if not exists:
            writer.writeheader()
        for row in rows:
            writer.writerow(row)


def run_command(command, cwd, timeout):
    started = time.time()
    result = subprocess.run(
        command,
        cwd=str(cwd),
        text=True,
        capture_output=True,
        timeout=timeout,
    )
    return result, time.time() - started


def repair_structure(args, run_dir):
    pdb_path = ensure_exists(args.pdb, "Input PDB").resolve()
    repaired_dir = (run_dir / "repaired").resolve()
    repaired_dir.mkdir(parents=True, exist_ok=True)
    local_pdb = repaired_dir / pdb_path.name
    repaired_pdb = repaired_dir / f"{pdb_path.stem}_Repair.pdb"

    if repaired_pdb.exists() and not args.force_repair:
        return repaired_pdb

    shutil.copy2(pdb_path, local_pdb)
    command = [
        str(Path(args.foldx).resolve()),
        "--command=RepairPDB",
        f"--pdb={local_pdb.name}",
        f"--pdb-dir={repaired_dir.resolve()}",
        f"--output-dir={repaired_dir.resolve()}",
        f"--rotabaseLocation={Path(args.rotabase).resolve()}",
        "--screen=0",
    ]
    result, elapsed = run_command(command, Path(args.foldx).resolve().parent, args.repair_timeout)
    log = repaired_dir / "repair_stdout_stderr.log"
    log.write_text(result.stdout + "\n\nSTDERR\n" + result.stderr, encoding="utf-8")

    if not repaired_pdb.exists():
        raise RuntimeError(
            f"FoldX RepairPDB did not produce {repaired_pdb}. See {log}. Elapsed: {elapsed:.1f}s"
        )
    return repaired_pdb


def parse_dif_fxout(path, batch_df, batch_id, elapsed_seconds):
    lines = Path(path).read_text(encoding="utf-8", errors="replace").splitlines()
    header_idx = None
    for idx, line in enumerate(lines):
        if line.startswith("Pdb\t"):
            header_idx = idx
            break
    if header_idx is None:
        raise ValueError(f"No energy table header found in {path}")

    table_lines = [line for line in lines[header_idx:] if line.strip()]
    header = table_lines[0].split("\t")
    data_rows = [line.split("\t") for line in table_lines[1:]]
    if len(data_rows) < len(batch_df):
        raise ValueError(f"Expected {len(batch_df)} DDG rows, found {len(data_rows)} in {path}")

    rows = []
    for mutation_row, values in zip(batch_df.to_dict("records"), data_rows):
        record = dict(zip(header, values))
        out = {
            "mutation": mutation_row["mutation"],
            "foldx_mutation": mutation_row["foldx_mutation"],
            "position": int(mutation_row["position"]),
            "foldx_position": int(mutation_row["foldx_position"]),
            "wt_aa": mutation_row["wt_aa"],
            "mut_aa": mutation_row["mut_aa"],
            "ddg": float(record["total energy"]),
            "status": "ok",
            "batch_id": batch_id,
            "foldx_pdb_row": record.get("Pdb", ""),
            "elapsed_seconds": round(elapsed_seconds, 3),
            "error": "",
        }
        for col in ENERGY_COLUMNS:
            if col == "total energy":
                continue
            out[f"foldx_{col.replace(' ', '_')}"] = record.get(col, "")
        rows.append(out)
    return rows


def run_batch(args, repaired_pdb, batch_id, batch_df):
    batch_dir = (Path(args.run_dir) / "batches" / f"batch_{batch_id:05d}").resolve()
    batch_dir.mkdir(parents=True, exist_ok=True)
    mutation_file = batch_dir / f"individual_list_batch_{batch_id:05d}.txt"
    mutation_file.write_text("\n".join(batch_df["foldx_mutation"].tolist()) + "\n", encoding="utf-8")

    command = [
        str(Path(args.foldx).resolve()),
        "--command=BuildModel",
        f"--pdb={repaired_pdb.name}",
        f"--pdb-dir={repaired_pdb.parent}",
        f"--mutant-file={mutation_file.resolve()}",
        f"--output-dir={batch_dir.resolve()}",
        f"--rotabaseLocation={Path(args.rotabase).resolve()}",
        "--numberOfRuns=1",
        "--out-pdb=0",
        "--screen=0",
    ]
    try:
        result, elapsed = run_command(command, Path(args.foldx).resolve().parent, args.batch_timeout)
        (batch_dir / "foldx_stdout_stderr.log").write_text(
            result.stdout + "\n\nSTDERR\n" + result.stderr,
            encoding="utf-8",
        )
        dif_path = batch_dir / f"Dif_{repaired_pdb.stem}.fxout"
        if not dif_path.exists():
            raise RuntimeError("FoldX did not produce Dif fxout")
        if "PROBLEM" in result.stdout or "There was a problem" in result.stdout:
            raise RuntimeError("FoldX reported a BuildModel problem")
        rows = parse_dif_fxout(dif_path, batch_df, batch_id, elapsed)
        return batch_id, rows
    except Exception as exc:
        rows = []
        for mutation_row in batch_df.to_dict("records"):
            rows.append(
                {
                    "mutation": mutation_row["mutation"],
                    "foldx_mutation": mutation_row["foldx_mutation"],
                    "position": int(mutation_row["position"]),
                    "foldx_position": int(mutation_row["foldx_position"]),
                    "wt_aa": mutation_row["wt_aa"],
                    "mut_aa": mutation_row["mut_aa"],
                    "ddg": "",
                    "status": "failed",
                    "batch_id": batch_id,
                    "foldx_pdb_row": "",
                    "elapsed_seconds": "",
                    "error": str(exc),
                }
            )
        return batch_id, rows


def make_batches(df, batch_size):
    for idx, start in enumerate(range(0, len(df), batch_size), start=1):
        yield idx, df.iloc[start : start + batch_size].copy()


def merge_results(source_feature_csv, results_csv, output_csv):
    base = pd.read_csv(source_feature_csv)
    ddg = pd.read_csv(results_csv)
    ok = ddg.loc[ddg["status"] == "ok"].copy()
    ok = ok.sort_values(["mutation", "batch_id"]).drop_duplicates("mutation", keep="last")
    merge_cols = ["mutation", "ddg", "foldx_mutation"] + [
        col for col in ok.columns if col.startswith("foldx_") and col != "foldx_mutation"
    ]
    merged = base.merge(ok[merge_cols], on="mutation", how="left")
    merged["ddg_label_available"] = merged["ddg"].notna().astype(int)
    merged["destabilizing_ddg_label"] = (pd.to_numeric(merged["ddg"], errors="coerce") >= 1.0).astype("Int64")
    merged.to_csv(output_csv, index=False)
    return merged


def write_manifest(args, run_dir, repaired_pdb, total, pending):
    manifest = {
        "foldx": str(Path(args.foldx).resolve()),
        "rotabase": str(Path(args.rotabase).resolve()),
        "input_pdb": str(Path(args.pdb).resolve()),
        "repaired_pdb": str(repaired_pdb.resolve()),
        "mutation_csv": str(Path(args.mutation_csv).resolve()),
        "batch_size": args.batch_size,
        "workers": args.workers,
        "total_mutations_requested": int(total),
        "pending_at_start": int(pending),
        "ddg_sign_convention": "FoldX BuildModel Dif total energy; positive values usually indicate destabilization.",
        "alpha_fold_caveat": "DDG values are computed on an AlphaFold model and should be treated as approximate computational labels.",
    }
    (run_dir / "foldx_run_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")


def main():
    args = parse_args()
    ensure_exists(args.foldx, "FoldX executable")
    ensure_exists(args.rotabase, "FoldX rotabase")
    run_dir = Path(args.run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)

    mutation_df = load_mutations(args.mutation_csv, args.limit)
    repaired_pdb = repair_structure(args, run_dir)
    completed = load_completed(args.results_csv)
    pending_df = mutation_df.loc[~mutation_df["mutation"].astype(str).isin(completed)].copy()
    write_manifest(args, run_dir, repaired_pdb, len(mutation_df), len(pending_df))

    if pending_df.empty:
        print("No pending mutations. Merging existing DDG results.")
    else:
        print("Total mutations:", len(mutation_df))
        print("Completed mutations:", len(completed))
        print("Pending mutations:", len(pending_df))
        print("Batch size:", args.batch_size)
        print("Workers:", args.workers)

        batches = list(make_batches(pending_df, args.batch_size))
        if args.workers <= 1:
            for batch_id, batch_df in batches:
                _, rows = run_batch(args, repaired_pdb, batch_id, batch_df)
                append_rows(args.results_csv, rows)
                ok_count = sum(row["status"] == "ok" for row in rows)
                print(f"Batch {batch_id:05d}: {ok_count}/{len(rows)} ok")
        else:
            with ThreadPoolExecutor(max_workers=args.workers) as executor:
                futures = {
                    executor.submit(run_batch, args, repaired_pdb, batch_id, batch_df): batch_id
                    for batch_id, batch_df in batches
                }
                for future in as_completed(futures):
                    batch_id, rows = future.result()
                    append_rows(args.results_csv, rows)
                    ok_count = sum(row["status"] == "ok" for row in rows)
                    print(f"Batch {batch_id:05d}: {ok_count}/{len(rows)} ok")

    merged = merge_results(args.source_feature_csv, args.results_csv, args.merged_output_csv)
    available = int(merged["ddg_label_available"].sum())
    print("DDG labels available:", available)
    print("Merged dataset:", args.merged_output_csv)


if __name__ == "__main__":
    main()
