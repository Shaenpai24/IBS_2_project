import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    f1_score,
    mean_absolute_error,
    mean_squared_error,
    r2_score,
    roc_auc_score,
)
from sklearn.model_selection import GroupKFold

from train_mutation_impact_models import add_proxy_targets, feature_matrix


DEFAULT_INPUT = "prnp_mutation_features_tier2_structural.csv"
DEFAULT_OUTPUT_DIR = "results_proxy_leakage"
RANDOM_STATE = 42


PROXY_FORMULA_FEATURES = {
    "chemical": [
        "grantham",
        "blosum62",
        "delta_hydropathy",
        "delta_volume",
        "delta_charge",
        "delta_polarity_flag",
    ],
    "conservation": [
        "conservation_wt_freq",
        "position_entropy",
    ],
    "structural": [
        "af2_neighbor_count_8A",
        "af2_is_buried_proxy",
        "is_globular_domain_region",
        "af2_plddt",
        "delta_volume",
    ],
}


FEATURE_GROUPS = {
    "biochemical_only": [
        "grantham",
        "blosum62",
        "delta_hydropathy",
        "delta_volume",
        "delta_charge",
        "wt_is_polar",
        "mut_is_polar",
        "delta_polarity_flag",
        "delta_instability_index",
        "delta_gravy",
        "delta_pI",
        "delta_net_charge",
        "delta_aromaticity",
        "delta_aliphatic_index",
        "delta_hydrophobic_stretch_count",
        "delta_octapeptide_repeat_count",
        "delta_n_glyco_motif_count",
    ],
    "conservation_only": [
        "conservation_wt_freq",
        "conservation_mut_freq",
        "conservation_delta_freq",
        "position_entropy",
        "position_observation_count",
    ],
    "structural_only": [
        "af2_mapped",
        "af2_struct_seq_id",
        "af2_plddt",
        "af2_plddt_window_mean_5",
        "af2_plddt_window_min_5",
        "af2_neighbor_count_8A",
        "af2_neighbor_count_12A",
        "af2_mean_neighbor_distance_8A",
        "af2_min_neighbor_distance",
        "af2_is_buried_proxy",
        "af2_is_exposed_proxy",
        "af2_is_low_confidence",
        "af2_is_high_confidence",
        "is_signal_peptide_region",
        "is_octapeptide_repeat_region",
        "is_globular_domain_region",
        "is_gpi_anchor_region",
    ],
}


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Audit circularity and proxy bias in heuristic PRNP mutation severity targets. "
            "Runs grouped cross-validated ablations by biochemical, conservation, and "
            "structural feature families."
        )
    )
    parser.add_argument("--input_csv", default=DEFAULT_INPUT)
    parser.add_argument("--output_dir", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--n_estimators", type=int, default=300)
    return parser.parse_args()


def available_columns(df, columns):
    return [col for col in columns if col in df.columns]


def build_group_matrix(df, group_name):
    if group_name == "combined_all_numeric":
        return feature_matrix(df)
    if group_name == "combined_no_formula_features":
        direct = set(sum(PROXY_FORMULA_FEATURES.values(), []))
        return feature_matrix(df).drop(columns=[col for col in direct if col in df.columns])

    columns = available_columns(df, FEATURE_GROUPS[group_name])
    if not columns:
        raise ValueError(f"No available columns for group: {group_name}")

    X = df[columns].apply(pd.to_numeric, errors="coerce")
    X = X.replace([np.inf, -np.inf], np.nan)
    X = X.fillna(X.median(numeric_only=True)).fillna(0)
    return X


def evaluate_group(X, y_class, y_reg, groups, folds, n_estimators):
    n_splits = min(folds, int(pd.Series(groups).nunique()))
    splitter = GroupKFold(n_splits=n_splits)

    class_rows = []
    reg_rows = []

    for fold, (train_idx, test_idx) in enumerate(splitter.split(X, y_class, groups=groups), start=1):
        X_train = X.iloc[train_idx]
        X_test = X.iloc[test_idx]
        y_class_train = y_class.iloc[train_idx]
        y_class_test = y_class.iloc[test_idx]
        y_reg_train = y_reg.iloc[train_idx]
        y_reg_test = y_reg.iloc[test_idx]

        clf = RandomForestClassifier(
            n_estimators=n_estimators,
            random_state=RANDOM_STATE + fold,
            class_weight="balanced",
            n_jobs=-1,
            min_samples_leaf=2,
        )
        reg = RandomForestRegressor(
            n_estimators=n_estimators,
            random_state=RANDOM_STATE + fold,
            n_jobs=-1,
            min_samples_leaf=2,
        )

        clf.fit(X_train, y_class_train)
        reg.fit(X_train, y_reg_train)

        pred = clf.predict(X_test)
        prob = clf.predict_proba(X_test)[:, 1]
        reg_pred = reg.predict(X_test)

        class_rows.append(
            {
                "fold": fold,
                "accuracy": accuracy_score(y_class_test, pred),
                "f1": f1_score(y_class_test, pred, zero_division=0),
                "roc_auc": roc_auc_score(y_class_test, prob),
                "average_precision": average_precision_score(y_class_test, prob),
            }
        )
        reg_rows.append(
            {
                "fold": fold,
                "mae": mean_absolute_error(y_reg_test, reg_pred),
                "rmse": float(np.sqrt(mean_squared_error(y_reg_test, reg_pred))),
                "r2": r2_score(y_reg_test, reg_pred),
            }
        )

    return pd.DataFrame(class_rows), pd.DataFrame(reg_rows)


def summarize(metric_df, group_name, target_type, feature_count):
    metric_cols = [col for col in metric_df.columns if col != "fold"]
    row = {
        "feature_group": group_name,
        "target_type": target_type,
        "feature_count": feature_count,
    }
    for col in metric_cols:
        row[f"{col}_mean"] = metric_df[col].mean()
        row[f"{col}_std"] = metric_df[col].std(ddof=0)
    return row


def correlation_report(df):
    target_cols = [
        "chemical_impact_score",
        "evolutionary_importance_score",
        "structural_sensitivity_score",
        "mutation_severity_score",
        "destabilization_proxy_score",
        "high_severity_label",
        "destabilizing_label",
    ]
    direct_features = sorted(set(sum(PROXY_FORMULA_FEATURES.values(), [])))
    rows = []
    for feature in direct_features:
        if feature not in df.columns:
            continue
        for target in target_cols:
            corr = pd.to_numeric(df[feature], errors="coerce").corr(
                pd.to_numeric(df[target], errors="coerce"), method="spearman"
            )
            rows.append(
                {
                    "feature": feature,
                    "target": target,
                    "spearman_correlation": corr,
                    "directly_used_in_proxy_formula": True,
                }
            )
    return pd.DataFrame(rows).sort_values(
        "spearman_correlation", key=lambda s: s.abs(), ascending=False
    )


def write_leakage_report(df, ablation_summary, output_path):
    direct = sorted(set(sum(PROXY_FORMULA_FEATURES.values(), [])))
    table = ablation_summary.to_csv(index=False)
    with open(output_path, "w", encoding="utf-8") as handle:
        handle.write("# Proxy Leakage and Bias Audit\n\n")
        handle.write("## Direct Formula Features\n\n")
        handle.write(
            "The current heuristic targets are explicitly constructed from these input features. "
            "High model performance against proxy targets is therefore expected and should not be "
            "interpreted as biological validation.\n\n"
        )
        for group, features in PROXY_FORMULA_FEATURES.items():
            handle.write(f"- {group}: {', '.join(features)}\n")
        handle.write("\n## Circularity Assessment\n\n")
        handle.write(
            "The `mutation_severity_score`, `destabilization_proxy_score`, and derived binary "
            "labels are circular with the features listed above. Models trained to predict these "
            "targets measure how well an algorithm can reconstruct the engineered formula, not "
            "whether it predicts experimental pathogenicity or protein stability.\n\n"
        )
        handle.write("## Ablation Summary\n\n")
        handle.write("```csv\n")
        handle.write(table)
        handle.write("```\n")
        handle.write("\n\n## Interpretation\n\n")
        handle.write(
            "- If a single feature group performs very well, the proxy target is dominated by that group.\n"
            "- If `combined_no_formula_features` drops sharply, the original performance was mostly formula leakage.\n"
            "- The next scientifically grounded step is to replace these targets with external ΔΔG or curated variant labels.\n"
        )


def main():
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.input_csv)
    df = add_proxy_targets(df)
    groups = df["position"]

    group_names = [
        "biochemical_only",
        "conservation_only",
        "structural_only",
        "combined_all_numeric",
        "combined_no_formula_features",
    ]

    fold_rows = []
    summary_rows = []
    for group_name in group_names:
        X = build_group_matrix(df, group_name)
        class_df, reg_df = evaluate_group(
            X,
            df["high_severity_label"],
            df["mutation_severity_score"],
            groups,
            args.folds,
            args.n_estimators,
        )
        class_df.insert(0, "feature_group", group_name)
        reg_df.insert(0, "feature_group", group_name)
        class_df.insert(1, "target_type", "high_severity_classification")
        reg_df.insert(1, "target_type", "severity_score_regression")
        fold_rows.extend([class_df, reg_df])
        summary_rows.append(
            summarize(class_df.drop(columns=["feature_group", "target_type"]), group_name, "high_severity_classification", X.shape[1])
        )
        summary_rows.append(
            summarize(reg_df.drop(columns=["feature_group", "target_type"]), group_name, "severity_score_regression", X.shape[1])
        )

    fold_metrics = pd.concat(fold_rows, ignore_index=True)
    ablation_summary = pd.DataFrame(summary_rows)
    correlations = correlation_report(df)

    fold_metrics.to_csv(output_dir / "ablation_fold_metrics.csv", index=False)
    ablation_summary.to_csv(output_dir / "ablation_summary.csv", index=False)
    correlations.to_csv(output_dir / "proxy_feature_target_correlations.csv", index=False)

    manifest = {
        "input_csv": args.input_csv,
        "rows": int(len(df)),
        "folds": int(args.folds),
        "grouped_by": "position",
        "proxy_formula_features": PROXY_FORMULA_FEATURES,
        "feature_groups": FEATURE_GROUPS,
        "warning": (
            "This audit evaluates proxy target circularity. It does not validate biological "
            "pathogenicity or experimental stability."
        ),
    }
    with open(output_dir / "proxy_audit_manifest.json", "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)
    write_leakage_report(df, ablation_summary, output_dir / "proxy_leakage_report.md")

    print("Rows:", len(df))
    print("Wrote:", output_dir)
    print("Ablation groups:", ", ".join(group_names))


if __name__ == "__main__":
    main()
