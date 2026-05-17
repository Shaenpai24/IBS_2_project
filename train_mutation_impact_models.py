import argparse
import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.inspection import permutation_importance
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    classification_report,
    confusion_matrix,
    mean_absolute_error,
    mean_squared_error,
    precision_recall_fscore_support,
    r2_score,
    roc_auc_score,
)
from sklearn.model_selection import GroupShuffleSplit
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


DEFAULT_INPUT = "prnp_mutation_features_tier2_structural.csv"
DEFAULT_OUTPUT_DIR = "results_mutation_impact"
RANDOM_STATE = 42


METADATA_COLUMNS = {
    "mutation",
    "wt_aa",
    "mut_aa",
    "local_window_wt",
    "af2_ref_aa",
    "af2_struct_aa",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Train baseline mutation-impact models for the PRNP mutation atlas. "
            "Targets are transparent proxy labels derived from biochemical, "
            "evolutionary, and AlphaFold structural features."
        )
    )
    parser.add_argument("--input_csv", default=DEFAULT_INPUT)
    parser.add_argument("--output_dir", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--test_size", type=float, default=0.2)
    parser.add_argument("--n_estimators", type=int, default=600)
    return parser.parse_args()


def minmax(series):
    values = pd.to_numeric(series, errors="coerce").astype(float)
    lo = values.min()
    hi = values.max()
    if not np.isfinite(lo) or not np.isfinite(hi) or hi == lo:
        return pd.Series(np.zeros(len(values)), index=values.index)
    return (values - lo) / (hi - lo)


def inverse_minmax(series):
    return 1.0 - minmax(series)


def require_columns(df, columns):
    missing = [col for col in columns if col not in df.columns]
    if missing:
        raise ValueError(f"Input dataset is missing required columns: {missing}")


def add_proxy_targets(df):
    required = [
        "grantham",
        "blosum62",
        "delta_hydropathy",
        "delta_volume",
        "delta_charge",
        "delta_polarity_flag",
        "conservation_wt_freq",
        "position_entropy",
        "af2_neighbor_count_8A",
        "af2_plddt",
        "af2_is_buried_proxy",
        "is_globular_domain_region",
    ]
    require_columns(df, required)

    work = df.copy()
    abs_charge = pd.to_numeric(work["delta_charge"], errors="coerce").abs()

    chemical = (
        0.35 * minmax(work["grantham"])
        + 0.20 * inverse_minmax(work["blosum62"])
        + 0.18 * minmax(pd.to_numeric(work["delta_hydropathy"], errors="coerce").abs())
        + 0.18 * minmax(pd.to_numeric(work["delta_volume"], errors="coerce").abs())
        + 0.06 * minmax(abs_charge)
        + 0.03 * minmax(work["delta_polarity_flag"])
    )

    evolutionary = (
        0.70 * minmax(work["conservation_wt_freq"])
        + 0.30 * inverse_minmax(work["position_entropy"])
    )

    structural = (
        0.35 * minmax(work["af2_neighbor_count_8A"])
        + 0.20 * minmax(work["af2_is_buried_proxy"])
        + 0.20 * minmax(work["is_globular_domain_region"])
        + 0.15 * minmax(work["af2_plddt"])
        + 0.10 * minmax(pd.to_numeric(work["delta_volume"], errors="coerce").abs())
    )

    severity = 100.0 * (
        0.42 * chemical
        + 0.30 * evolutionary
        + 0.28 * structural
    )

    destabilization = 100.0 * (
        0.35 * structural
        + 0.25 * minmax(pd.to_numeric(work["delta_volume"], errors="coerce").abs())
        + 0.20 * minmax(pd.to_numeric(work["delta_hydropathy"], errors="coerce").abs())
        + 0.12 * inverse_minmax(work["blosum62"])
        + 0.08 * minmax(work["grantham"])
    )

    work["chemical_impact_score"] = chemical.clip(0, 1) * 100
    work["evolutionary_importance_score"] = evolutionary.clip(0, 1) * 100
    work["structural_sensitivity_score"] = structural.clip(0, 1) * 100
    work["mutation_severity_score"] = severity.clip(0, 100)
    work["destabilization_proxy_score"] = destabilization.clip(0, 100)

    high_cutoff = work["mutation_severity_score"].quantile(0.75)
    moderate_cutoff = work["mutation_severity_score"].quantile(0.50)
    work["severity_class"] = np.select(
        [
            work["mutation_severity_score"] >= high_cutoff,
            work["mutation_severity_score"] >= moderate_cutoff,
        ],
        ["high", "moderate"],
        default="low",
    )
    work["high_severity_label"] = (work["mutation_severity_score"] >= high_cutoff).astype(int)
    work["destabilizing_label"] = (
        work["destabilization_proxy_score"]
        >= work["destabilization_proxy_score"].quantile(0.75)
    ).astype(int)

    return work


def feature_matrix(df):
    numeric = df.select_dtypes(include=[np.number]).copy()
    target_cols = [
        "chemical_impact_score",
        "evolutionary_importance_score",
        "structural_sensitivity_score",
        "mutation_severity_score",
        "destabilization_proxy_score",
        "high_severity_label",
        "destabilizing_label",
    ]
    numeric = numeric.drop(columns=[col for col in target_cols if col in numeric.columns])
    numeric = numeric.replace([np.inf, -np.inf], np.nan)
    numeric = numeric.fillna(numeric.median(numeric_only=True))
    numeric = numeric.fillna(0)
    return numeric


def grouped_split(X, y, groups, test_size):
    splitter = GroupShuffleSplit(n_splits=1, test_size=test_size, random_state=RANDOM_STATE)
    train_idx, test_idx = next(splitter.split(X, y, groups=groups))
    return train_idx, test_idx


def classifier_metrics(model, X_test, y_test):
    pred = model.predict(X_test)
    prob = model.predict_proba(X_test)[:, 1]
    precision, recall, f1, _ = precision_recall_fscore_support(
        y_test, pred, average="binary", zero_division=0
    )
    return {
        "accuracy": float(accuracy_score(y_test, pred)),
        "precision": float(precision),
        "recall": float(recall),
        "f1": float(f1),
        "roc_auc": float(roc_auc_score(y_test, prob)),
        "average_precision": float(average_precision_score(y_test, prob)),
        "confusion_matrix": confusion_matrix(y_test, pred).tolist(),
        "classification_report": classification_report(y_test, pred, zero_division=0),
    }


def regressor_metrics(model, X_test, y_test):
    pred = model.predict(X_test)
    rmse = float(np.sqrt(mean_squared_error(y_test, pred)))
    return {
        "mae": float(mean_absolute_error(y_test, pred)),
        "rmse": rmse,
        "r2": float(r2_score(y_test, pred)),
    }


def save_feature_importance(model, X_test, y_test, output_csv):
    result = permutation_importance(
        model,
        X_test,
        y_test,
        n_repeats=12,
        random_state=RANDOM_STATE,
        n_jobs=-1,
    )
    importance = pd.DataFrame(
        {
            "feature": X_test.columns,
            "importance_mean": result.importances_mean,
            "importance_std": result.importances_std,
        }
    ).sort_values("importance_mean", ascending=False)
    importance.to_csv(output_csv, index=False)
    return importance


def main():
    args = parse_args()
    input_path = Path(args.input_csv)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(input_path)
    df = add_proxy_targets(df)
    X = feature_matrix(df)
    groups = df["position"]

    train_idx, test_idx = grouped_split(
        X,
        df["high_severity_label"],
        groups=groups,
        test_size=args.test_size,
    )

    X_train = X.iloc[train_idx]
    X_test = X.iloc[test_idx]

    clf = Pipeline(
        [
            ("scaler", StandardScaler()),
            (
                "model",
                RandomForestClassifier(
                    n_estimators=args.n_estimators,
                    random_state=RANDOM_STATE,
                    class_weight="balanced",
                    n_jobs=-1,
                    min_samples_leaf=2,
                ),
            ),
        ]
    )
    reg = Pipeline(
        [
            ("scaler", StandardScaler()),
            (
                "model",
                RandomForestRegressor(
                    n_estimators=args.n_estimators,
                    random_state=RANDOM_STATE,
                    n_jobs=-1,
                    min_samples_leaf=2,
                ),
            ),
        ]
    )

    y_clf = df["high_severity_label"]
    y_reg = df["mutation_severity_score"]

    clf.fit(X_train, y_clf.iloc[train_idx])
    reg.fit(X_train, y_reg.iloc[train_idx])

    metrics = {
        "input_csv": str(input_path),
        "rows": int(len(df)),
        "features": int(X.shape[1]),
        "train_rows": int(len(train_idx)),
        "test_rows": int(len(test_idx)),
        "grouped_by": "position",
        "target_note": (
            "Targets are heuristic proxy labels/scores derived from existing mutation features. "
            "Replace or calibrate them with experimental pathogenicity or ddG labels when available."
        ),
        "high_severity_classifier": classifier_metrics(clf, X_test, y_clf.iloc[test_idx]),
        "severity_regressor": regressor_metrics(reg, X_test, y_reg.iloc[test_idx]),
    }

    clf_prob = clf.predict_proba(X)[:, 1]
    reg_score = reg.predict(X)
    ranked = df.copy()
    ranked["predicted_high_severity_probability"] = clf_prob
    ranked["predicted_mutation_severity_score"] = np.clip(reg_score, 0, 100)
    ranked["predicted_destabilization_probability"] = (
        ranked["predicted_high_severity_probability"] * 0.45
        + minmax(ranked["destabilization_proxy_score"]) * 0.55
    ).clip(0, 1)

    rank_columns = [
        "mutation",
        "position",
        "wt_aa",
        "mut_aa",
        "severity_class",
        "mutation_severity_score",
        "predicted_mutation_severity_score",
        "predicted_high_severity_probability",
        "destabilization_proxy_score",
        "predicted_destabilization_probability",
        "chemical_impact_score",
        "evolutionary_importance_score",
        "structural_sensitivity_score",
        "grantham",
        "blosum62",
        "delta_hydropathy",
        "delta_volume",
        "conservation_wt_freq",
        "position_entropy",
        "af2_plddt",
        "af2_neighbor_count_8A",
        "af2_is_buried_proxy",
        "is_globular_domain_region",
    ]
    ranked = ranked.sort_values(
        ["predicted_mutation_severity_score", "predicted_high_severity_probability"],
        ascending=False,
    )
    ranked[rank_columns].to_csv(output_dir / "mutation_impact_ranked.csv", index=False)
    df.to_csv(output_dir / "mutation_impact_features_with_proxy_targets.csv", index=False)

    save_feature_importance(
        clf,
        X_test,
        y_clf.iloc[test_idx],
        output_dir / "feature_importance_high_severity.csv",
    )
    save_feature_importance(
        reg,
        X_test,
        y_reg.iloc[test_idx],
        output_dir / "feature_importance_severity_score.csv",
    )

    joblib.dump(clf, output_dir / "random_forest_high_severity.joblib")
    joblib.dump(reg, output_dir / "random_forest_severity_score.joblib")

    with open(output_dir / "model_metrics.json", "w", encoding="utf-8") as handle:
        json.dump(metrics, handle, indent=2)

    print("Rows:", len(df))
    print("Numeric features:", X.shape[1])
    print("Train rows:", len(train_idx))
    print("Test rows:", len(test_idx))
    print("Classifier ROC-AUC:", round(metrics["high_severity_classifier"]["roc_auc"], 4))
    print("Regressor R2:", round(metrics["severity_regressor"]["r2"], 4))
    print("Saved outputs to:", output_dir)


if __name__ == "__main__":
    main()
