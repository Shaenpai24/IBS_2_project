# DDG Generation Protocol

This folder contains mutation-list files for external protein stability engines.

## Inputs

- Structure: `AF-P04156-F1-model_v4.pdb`
- Chain assumed for mutation commands: `A`
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
