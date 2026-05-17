import argparse
from pathlib import Path

import pandas as pd


DEFAULT_INPUT = "prnp_mutation_features_tier2_structural.csv"
DEFAULT_OUTPUT_DIR = "ddg_inputs"
DEFAULT_MERGED_OUTPUT = "prnp_mutation_features_ddg.csv"


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Prepare PRNP mutation lists for external DDG tools and optionally merge "
            "computed DDG labels back into the mutation feature table."
        )
    )
    parser.add_argument("--input_csv", default=DEFAULT_INPUT)
    parser.add_argument("--output_dir", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--chain", default="A")
    parser.add_argument(
        "--ddg_results_csv",
        default=None,
        help=(
            "Optional CSV containing mutation-level DDG results. It must contain a "
            "`mutation` column and one DDG value column."
        ),
    )
    parser.add_argument(
        "--ddg_column",
        default="ddg",
        help="Column name in --ddg_results_csv containing DDG values.",
    )
    parser.add_argument("--merged_output_csv", default=DEFAULT_MERGED_OUTPUT)
    return parser.parse_args()


def require_columns(df, columns):
    missing = [col for col in columns if col not in df.columns]
    if missing:
        raise ValueError(f"Input is missing required columns: {missing}")


def foldx_position(row):
    if "af2_struct_seq_id" in row and pd.notna(row["af2_struct_seq_id"]):
        return int(float(row["af2_struct_seq_id"]))
    return int(row["position"])


def foldx_mutation(row, chain):
    return f"{row.wt_aa}{chain}{foldx_position(row)}{row.mut_aa};"


def rosetta_mutation_line(row):
    return f"{row.wt_aa} {foldx_position(row)} {row.mut_aa}"


def write_rosetta_mutfiles(df, output_dir):
    mutfile_dir = output_dir / "rosetta_mutfiles"
    mutfile_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for row in df.itertuples(index=False):
        path = mutfile_dir / f"{row.mutation}.mutfile"
        position = int(getattr(row, "foldx_position", getattr(row, "position")))
        text = f"total 1\n1\n{row.wt_aa} {position} {row.mut_aa}\n"
        path.write_text(text, encoding="utf-8")
        rows.append({"mutation": row.mutation, "rosetta_mutfile": str(path)})
    return pd.DataFrame(rows)


def write_ddg_template(df, output_dir, chain):
    base_columns = ["mutation", "position", "wt_aa", "mut_aa"]
    if "af2_struct_seq_id" in df.columns:
        base_columns.append("af2_struct_seq_id")
    template = df[base_columns].copy()
    template["foldx_position"] = template.apply(foldx_position, axis=1)
    template["chain"] = chain
    template["foldx_mutation"] = template.apply(lambda row: foldx_mutation(row, chain), axis=1)
    template["rosetta_mutation"] = template.apply(rosetta_mutation_line, axis=1)
    template["ddg"] = pd.NA
    template["ddg_tool"] = pd.NA
    template["ddg_units"] = "kcal/mol"
    template["ddg_sign_convention"] = "positive values usually indicate destabilization"
    template["ddg_notes"] = pd.NA
    return template


def write_foldx_files(template, output_dir):
    foldx_dir = output_dir / "foldx"
    foldx_dir.mkdir(parents=True, exist_ok=True)
    mutation_lines = template["foldx_mutation"].tolist()
    (foldx_dir / "individual_list.txt").write_text(
        "\n".join(mutation_lines) + "\n",
        encoding="utf-8",
    )
    template[["mutation", "foldx_mutation"]].to_csv(
        foldx_dir / "foldx_mutation_mapping.csv",
        index=False,
    )


def merge_ddg(base_df, ddg_results_csv, ddg_column):
    ddg_df = pd.read_csv(ddg_results_csv)
    require_columns(ddg_df, ["mutation", ddg_column])
    ddg_df = ddg_df[["mutation", ddg_column]].rename(columns={ddg_column: "ddg"})
    merged = base_df.merge(ddg_df, on="mutation", how="left", validate="one_to_one")
    merged["ddg_label_available"] = merged["ddg"].notna().astype(int)
    merged["destabilizing_ddg_label"] = (pd.to_numeric(merged["ddg"], errors="coerce") >= 1.0).astype("Int64")
    return merged


def write_protocol(output_dir, chain):
    protocol = f"""# DDG Generation Protocol

This folder contains mutation-list files for external protein stability engines.

## Inputs

- Structure: `AF-P04156-F1-model_v4.pdb`
- Chain assumed for mutation commands: `{chain}`
- Mutation count: 4,655 single-point substitutions

## FoldX

Generated file:

`foldx/individual_list.txt`

Typical FoldX workflow:

1. Repair the AlphaFold PDB with FoldX `RepairPDB`.
2. Run `BuildModel` with `--mutant-file foldx/individual_list.txt`.
3. Parse FoldX output energies into a CSV with at least:
   - `mutation`
   - `ddg`

Important: AlphaFold structures are predicted models. DDG values from FoldX on AlphaFold structures should be treated as approximate computational labels, not experimental measurements.

## Rosetta

Generated folder:

`rosetta_mutfiles/`

Each mutation has a one-mutation Rosetta mutfile. Rosetta protocols vary by installation and scoring function; keep protocol metadata with any generated DDG values.

## Merge Back

After external DDG calculation, create a CSV with:

`mutation,ddg`

Then run:

`python prepare_ddg_inputs.py --ddg_results_csv your_ddg_results.csv --ddg_column ddg`

This will create `prnp_mutation_features_ddg.csv`.
"""
    (output_dir / "DDG_PROTOCOL.md").write_text(protocol, encoding="utf-8")


def main():
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    base_df = pd.read_csv(args.input_csv)
    require_columns(base_df, ["mutation", "position", "wt_aa", "mut_aa"])

    template = write_ddg_template(base_df, output_dir, args.chain)
    write_foldx_files(template, output_dir)
    rosetta_map = write_rosetta_mutfiles(
        template[["mutation", "position", "foldx_position", "wt_aa", "mut_aa"]],
        output_dir,
    )
    template = template.merge(rosetta_map, on="mutation", how="left")
    template.to_csv(output_dir / "ddg_label_template.csv", index=False)
    write_protocol(output_dir, args.chain)

    if args.ddg_results_csv:
        merged = merge_ddg(base_df, args.ddg_results_csv, args.ddg_column)
        merged.to_csv(args.merged_output_csv, index=False)
        print("Merged DDG labels:", args.merged_output_csv)
        print("DDG labels available:", int(merged["ddg_label_available"].sum()))
    else:
        empty = base_df.copy()
        empty["ddg"] = pd.NA
        empty["ddg_label_available"] = 0
        empty["destabilizing_ddg_label"] = pd.NA
        empty.to_csv(args.merged_output_csv, index=False)
        print("Created DDG-ready dataset with empty labels:", args.merged_output_csv)

    print("Wrote DDG inputs to:", output_dir)
    print("Mutation count:", len(base_df))


if __name__ == "__main__":
    main()
