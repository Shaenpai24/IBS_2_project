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
- `train_mutation_impact_models.py` — Baseline mutation-impact modeling script. It builds transparent proxy targets from biochemical, conservation, and AlphaFold structural features, trains Random Forest classifier/regressor models, and exports ranked mutation-impact predictions.
- `validate_proxy_leakage.py` — Audits circularity in the heuristic proxy targets and runs grouped ablation experiments for biochemical-only, conservation-only, structural-only, combined, and formula-feature-removed models.
- `prepare_ddg_inputs.py` — Exports mutation lists/templates for external ΔΔG engines such as FoldX or Rosetta and can merge completed ΔΔG results into `prnp_mutation_features_ddg.csv`.
- `run_foldx_ddg_batches.py` — Resumable batched FoldX runner. It repairs the AlphaFold structure once, runs BuildModel in checkpointed mutation batches, skips completed mutations on resume, and merges parsed ΔΔG values into `prnp_mutation_features_ddg.csv`.
- `train_ddg_models_and_hotspots.py` — Trains biologically grounded ΔΔG regression and destabilization classification models, then generates mutation explanations, residue hotspot rankings, heatmaps, and structural sensitivity maps.
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

5) Baseline mutation-impact modeling
- `train_mutation_impact_models.py` creates interpretable proxy scores:
  - `chemical_impact_score`
  - `evolutionary_importance_score`
  - `structural_sensitivity_score`
  - `mutation_severity_score`
  - `destabilization_proxy_score`
- The script then trains:
  - a Random Forest classifier for high-severity mutation probability
  - a Random Forest regressor for mutation severity score
- The train/test split is grouped by residue position, so all mutations from a held-out position remain together in the test set.
- Outputs are saved under `results_mutation_impact/`, including ranked predictions, model metrics, feature importance tables, and serialized model files.

Important: these are heuristic proxy targets derived from existing features, not experimentally validated pathogenicity or ΔΔG labels. They are useful for ranking, exploration, and model plumbing, but should be calibrated against curated variant labels or stability measurements before biological claims are made.

6) Proxy leakage validation
- `validate_proxy_leakage.py` identifies the features directly used in heuristic target construction and runs ablation models with grouped cross-validation by residue position.
- The output folder `results_proxy_leakage/` contains:
  - `proxy_leakage_report.md`
  - `ablation_summary.csv`
  - `ablation_fold_metrics.csv`
  - `proxy_feature_target_correlations.csv`
  - `proxy_audit_manifest.json`
- This stage is intended to quantify circularity and proxy bias before extending the platform with biological labels.

7) ΔΔG label preparation
- `prepare_ddg_inputs.py` creates external-tool input files under `ddg_inputs/`:
  - `foldx/individual_list.txt`
  - `foldx/foldx_mutation_mapping.csv`
  - `rosetta_mutfiles/*.mutfile`
  - `ddg_label_template.csv`
  - `DDG_PROTOCOL.md`
- If a completed external ΔΔG result CSV is supplied, it merges values into `prnp_mutation_features_ddg.csv`.
- Without completed external ΔΔG values, the script creates a DDG-ready dataset with empty DDG label columns so downstream code has a stable schema.

8) FoldX ΔΔG generation and biological modeling
- FoldX executable used locally: `D:\FoldX\foldx_20261231.exe`.
- The AlphaFold structure is repaired once with FoldX `RepairPDB`; BuildModel is then run in resumable batches.
- Mutation commands use the AlphaFold/FoldX residue ID (`af2_struct_seq_id` / `foldx_position`), not only the biological reference position. This matters because the reference sequence and PDB numbering diverge after the signal peptide region.
- Final merged biological-label dataset: `prnp_mutation_features_ddg.csv`.
- Valid FoldX ΔΔG labels currently available: 4,636 / 4,655 mutations.
- Missing labels: 19 `V121*` mutations. The mutation table has wild-type `V` at reference position 121, while the AlphaFold structure has `M` at mapped structure residue 129. FoldX correctly rejects `VA129*` because residue 129 is methionine in this structure.
- Main downstream outputs are in `results_ddg_models/`:
  - `mutation_ddg_predictions_with_explanations.csv`
  - `residue_hotspot_rankings.csv`
  - `structural_sensitivity_map.csv`
  - `heatmap_observed_ddg.csv`
  - `heatmap_predicted_ddg.csv`
  - `heatmap_destabilizing_probability.csv`
  - `heatmap_predicted_ddg.png`
  - `heatmap_destabilizing_probability.png`
  - `ddg_cross_validation_summary.csv`
  - `ddg_holdout_metrics.csv`

Important: FoldX ΔΔG labels are computational labels generated from an AlphaFold model. They are more biologically grounded than engineered heuristic proxy targets, but they are still not experimental stability measurements.

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

Train baseline mutation-impact models
```
python train_mutation_impact_models.py
```

Audit proxy leakage and feature-family bias
```
python validate_proxy_leakage.py
```

Prepare external ΔΔG mutation inputs
```
python prepare_ddg_inputs.py
```

Run FoldX ΔΔG generation in resumable batches
```
python run_foldx_ddg_batches.py --batch_size 50 --workers 1
```

Merge completed ΔΔG labels
```
python prepare_ddg_inputs.py --ddg_results_csv your_ddg_results.csv --ddg_column ddg
```

Train ΔΔG models and generate hotspot/heatmap outputs
```
python train_ddg_models_and_hotspots.py
```

Main outputs
```
results_mutation_impact/mutation_impact_ranked.csv
results_mutation_impact/mutation_impact_features_with_proxy_targets.csv
results_mutation_impact/model_metrics.json
results_mutation_impact/feature_importance_high_severity.csv
results_mutation_impact/feature_importance_severity_score.csv
results_mutation_impact/random_forest_high_severity.joblib
results_mutation_impact/random_forest_severity_score.joblib
results_proxy_leakage/proxy_leakage_report.md
ddg_inputs/ddg_label_template.csv
prnp_mutation_features_ddg.csv
foldx_ddg_results.csv
results_ddg_models/mutation_ddg_predictions_with_explanations.csv
results_ddg_models/residue_hotspot_rankings.csv
results_ddg_models/heatmap_predicted_ddg.png
```

Quick demo (local viewer)
```
python -m http.server 8000
# then open http://127.0.0.1:8000/prnp_viewer.html
```

Suggested next steps
--------------------
- Add `requirements.txt` or `pyproject.toml` to lock dependencies and ease reproducibility.
- Replace or calibrate proxy targets with curated PRNP pathogenic/benign labels and/or experimentally measured or physics-derived ΔΔG values.
- Add mutation-aware ESM2 embeddings by embedding wild-type and mutant sequences and using embedding deltas as model features.
- Optionally run selective mutant structural modeling (ColabFold or Rosetta ddG) for a hand-picked shortlist of variants where structural detail matters.
- Build a feature registry that documents exact column names, types, ranges, and expected interpretation for each feature for downstream users.

Contact / Attribution
---------------------
If you use this work in a publication, please cite the repository and describe which AlphaFold model/version and inputs were used to compute structural proxies.

---
Generated by the project assistant — update as desired for publishing or internal reporting.
# BIO_PROJECT
