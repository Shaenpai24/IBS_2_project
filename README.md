# PRNP Mutational Features & AlphaFold Integration

This repository contains a reproducible pipeline for exhaustive in-silico single-residue mutagenesis of the human prion protein (PRNP), engineering of sequence-based (Tier-2) features, and integration of AlphaFold-derived structural proxies. It also includes an interactive HTML viewer for exploring mutation-level features alongside the AlphaFold structure.

Purpose
-------
- Generate all single amino-acid substitutions for a canonical PRNP reference sequence and compute engineered features useful for downstream modeling (variant effect prediction, pathogenicity prioritization).
- Augment sequence-level features with per-residue structural proxies derived from an AlphaFold prediction to capture local confidence, packing, and exposure.
- Provide interactive tools for rapid inspection and demonstration (local static viewer and GitHub Pages deployment artifacts).

Repository layout
-----------------
- `prnp_single_point_mutations.csv` — Frozen table of all generated single-point mutants (245 positions × 19 substitutions = 4,655 rows). Columns: `position`,`original_aa`,`mutated_aa`,`mutation`,`original_sequence`,`mutated_sequence`.
- `prnp_mutation_features_tier2.csv` — Tier-2 mutation-level engineered features derived from sequence, conservation, substitution matrices, and context windows (32 columns).
- `prnp_alphafold_position_features.csv` — Per-position structural proxy features computed from an AlphaFold PDB (245 rows, ~16 columns; columns prefixed with `af2_`).
- `prnp_mutation_features_tier2_structural.csv` — Merge of Tier-2 mutation features with AlphaFold-derived per-position proxies (4,655 rows, ~47 columns).
- `AF-P04156-F1-model_v4.pdb` — AlphaFold prediction used to derive structural proxies.
- `generate_prnp_alphafold_features.py` — Script to parse the AlphaFold structure, compute per-residue metrics (pLDDT windows, neighbor counts, exposure proxies), and merge into mutation-level features.
- `prnp_viewer.html` — Single-file interactive dashboard (CSV viewer + embedded 3Dmol.js visualiser) for exploring variants and highlighting residues on the AlphaFold model.
- `.github/workflows/deploy-pages.yml`, `index.html`, `.nojekyll` — Files to enable GitHub Pages hosting for the viewer.

Key numbers
-----------
- Reference length: 245 amino acids (human PRNP canonical sequence).
- Variants: 245 positions × 19 possible substitutions = 4,655 single-residue mutants (all rows present in `prnp_single_point_mutations.csv`).

How features were generated
--------------------------

1) In-silico mutagenesis
- A canonical reference protein sequence (FASTA) was selected as the authoritative PRNP reference.
- For each sequence position, all possible single amino-acid substitutions except the wild-type residue were generated programmatically and saved to `prnp_single_point_mutations.csv`.

2) Tier-2 sequence-derived features
- The `prnp_mutation_features_tier2.csv` table contains engineered features computed for each mutant. These include (non-exhaustive):
  - Conservation / alignment-derived scores (position-wise frequency, entropy).
  - Physicochemical change metrics between wild-type and mutant residues (hydrophobicity delta, charge change).
  - Substitution matrix scores (BLOSUM/other) for the given replacement.
  - Local sequence context windows mean/variance (sliding window statistics of properties like pI, hydropathy).
  - Global sequence properties that can affect folding (isoelectric point, molecular weight differences) when mutated.

Rationale: These features capture evolutionary signal, biochemical perturbation magnitude, and local context — all of which are informative for variant effect models.

3) AlphaFold-derived structural proxies
- Source: an AlphaFold predicted structure file `AF-P04156-F1-model_v4.pdb` (a local AlphaFold/ColabFold output saved in the repo).
- Per-residue metrics computed by `generate_prnp_alphafold_features.py` include:
  - `af2_plddt`: raw pLDDT score from AlphaFold per residue.
  - `af2_plddt_window_mean_5` / `af2_plddt_window_min_5`: rolling-window summary statistics to capture local confidence neighborhoods (window size 5 residues).
  - `af2_neighbor_count_8A`, `af2_neighbor_count_12A`: number of non-hydrogen atoms/residues within 8Å or 12Å of the residue Cα — a proxy for packing density.
  - `af2_mean_neighbor_distance_8A`, `af2_min_neighbor_distance`: short-range distance statistics capturing how close the residue is to nearby atoms.
  - `af2_is_buried_proxy`, `af2_is_exposed_proxy`: coarse exposure/burial indicators derived from neighbor counts and mean distances.
  - `af2_is_low_confidence`, `af2_is_high_confidence`: boolean flags for very low or very high pLDDT thresholds.

Rationale: AlphaFold outputs a per-residue confidence metric (pLDDT) and an atomic model which can be used to approximate packing and exposure. These structural proxies help distinguish residues that are buried and structurally constrained from flexible/exposed sites where mutations may be more tolerated.

4) Merging to mutation level
- Position-level structural proxies were joined with `prnp_mutation_features_tier2.csv` on the reference position to produce `prnp_mutation_features_tier2_structural.csv`. Each mutation row therefore contains both sequence-derived mutation features and structural context features for the mutated position.

Feature list and interpretation
-------------------------------
Below are the main feature groups and why they matter.

- Sequence / conservation features
  - `conservation_freq`, `entropy`: High conservation often implies functional/structural importance; mutations at conserved sites are more likely deleterious.

- Physicochemical delta features
  - `hydropathy_delta`, `charge_delta`, `polarity_delta`: These quantify biochemical property shifts that can disrupt local packing or interactions.

- Substitution matrix scores
  - `blosum62_score` (or equivalent): Higher scores indicate biochemically conservative substitutions; low/negative scores often correspond to disruptive changes.

- Context window statistics
  - `window_mean_hydropathy_5`, `window_var_hydropathy_11`: Capture whether the mutation sits in a hydrophobic core or polar loop; context matters more than single-residue changes in many cases.

- AlphaFold structural proxies (prefixed `af2_`)
  - `af2_plddt`: Low values suggest low confidence (intrinsic disorder or prediction uncertainty); mutations in low-pLDDT regions may be harder to interpret structurally.
  - `af2_plddt_window_mean_5`: Smooths pLDDT across neighbors — identifies locally confident regions.
  - `af2_neighbor_count_8A` / `af2_neighbor_count_12A`: High neighbour counts suggest buried/core residues; mutations here often destabilize structure.
  - `af2_mean_neighbor_distance_8A`: Smaller mean distances indicate tight packing; large changes in side-chain volume at such sites are likely disruptive.
  - `af2_is_buried_proxy` / `af2_is_exposed_proxy`: Useful categorical features for models and human interpretation.

Practical notes & limitations
---------------------------
- AlphaFold is an excellent predictor for likely folded regions but not an absolute ground truth: pLDDT signals confidence, not experimental stability. Use proxies as features, not definitive labels.
- Structural features are derived from a *single-model* prediction; large conformational changes or oligomerization contexts may not be represented.
- The neighbor-count exposure proxy is geometry-based and does not replace rigorous solvent-accessible surface area (SASA) measured from experimental structures.
- When using these features for downstream modeling (e.g., pathogenicity classifiers), consider calibrating thresholds and validating on known pathogenic/benign variants.

Reproducing the pipeline
------------------------
Minimal environment (suggested): Python 3.10+ with a virtual environment.

Install packages
```
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install biopython pandas numpy 3Dmoljs
```

Generate AlphaFold structural proxies (if you have a PDB/mmCIF)
```
python generate_prnp_alphafold_features.py \
  --af_structure AF-P04156-F1-model_v4.pdb \
  --reference_fasta prnp_reference_protein.fasta \
  --tier2_csv prnp_mutation_features_tier2.csv \
  --output_csv prnp_mutation_features_tier2_structural.csv \
  --position_csv prnp_alphafold_position_features.csv \
  --chain A
```

Quick demo (local viewer)
```
python -m http.server 8000
# then open http://127.0.0.1:8000/prnp_viewer.html
```

Suggested next steps
--------------------
- Add `requirements.txt` or `pyproject.toml` to lock dependencies and ease reproducibility.
- Optionally run selective mutant structural modeling (ColabFold or Rosetta ddG) for a hand-picked shortlist of variants where structural detail matters.
- Build a feature registry that documents exact column names, types, ranges, and expected interpretation for each feature for downstream users.

Contact / Attribution
---------------------
If you use this work in a publication, please cite the repository and describe which AlphaFold model/version and inputs were used to compute structural proxies.

---
Generated by the project assistant — update as desired for publishing or internal reporting.
# BIO_PROJECT
