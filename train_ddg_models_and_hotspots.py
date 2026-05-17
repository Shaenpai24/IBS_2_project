import argparse
import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier, HistGradientBoostingRegressor, RandomForestClassifier, RandomForestRegressor
from sklearn.inspection import permutation_importance
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    f1_score,
    mean_absolute_error,
    mean_squared_error,
    precision_score,
    r2_score,
    recall_score,
    roc_auc_score,
)
from sklearn.model_selection import GroupKFold, GroupShuffleSplit


DEFAULT_INPUT = "prnp_mutation_features_ddg.csv"
DEFAULT_OUTPUT_DIR = "results_ddg_models"
RANDOM_STATE = 42


METADATA_COLUMNS = {
    "mutation",
    "wt_aa",
    "mut_aa",
    "local_window_wt",
    "af2_ref_aa",
    "af2_struct_aa",
    "foldx_mutation",
}

TARGET_COLUMNS = {
    "ddg",
    "ddg_label_available",
    "destabilizing_ddg_label",
}

PROXY_TARGET_COLUMNS = {
    "chemical_impact_score",
    "evolutionary_importance_score",
    "structural_sensitivity_score",
    "mutation_severity_score",
    "destabilization_proxy_score",
    "high_severity_label",
    "destabilizing_label",
    "predicted_high_severity_probability",
    "predicted_mutation_severity_score",
    "predicted_destabilization_probability",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Train biologically grounded DDG prediction models and generate residue-level "
            "hotspot, heatmap, and explanation outputs."
        )
    )
    parser.add_argument("--input_csv", default=DEFAULT_INPUT)
    parser.add_argument("--output_dir", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--ddg_destabilizing_threshold", type=float, default=1.0)
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--n_estimators", type=int, default=600)
    return parser.parse_args()


def numeric_feature_matrix(df):
    numeric = df.select_dtypes(include=[np.number]).copy()
    leakage_cols = set()
    leakage_cols.update(TARGET_COLUMNS)
    leakage_cols.update(PROXY_TARGET_COLUMNS)
    leakage_cols.update(col for col in numeric.columns if col.startswith("foldx_"))
    leakage_cols.update(col for col in numeric.columns if col.lower().startswith("ddg"))
    X = numeric.drop(columns=[col for col in leakage_cols if col in numeric.columns])
    X = X.replace([np.inf, -np.inf], np.nan)
    X = X.fillna(X.median(numeric_only=True)).fillna(0)
    return X


def load_labeled_data(path, threshold):
    df = pd.read_csv(path)
    if "ddg" not in df.columns:
        raise ValueError("Input dataset does not contain required `ddg` column.")
    df["ddg"] = pd.to_numeric(df["ddg"], errors="coerce")
    labeled = df.loc[df["ddg"].notna()].copy()
    if labeled.empty:
        raise ValueError("No DDG labels are available yet. Run FoldX and merge DDG values first.")
    labeled["destabilizing_ddg_label"] = (labeled["ddg"] >= threshold).astype(int)
    return df, labeled


def grouped_holdout(X, y, groups):
    splitter = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=RANDOM_STATE)
    return next(splitter.split(X, y, groups=groups))


def cv_metrics(X, y_reg, y_class, groups, folds, n_estimators):
    n_splits = min(folds, int(pd.Series(groups).nunique()))
    splitter = GroupKFold(n_splits=n_splits)
    rows = []
    for fold, (train_idx, test_idx) in enumerate(splitter.split(X, y_reg, groups=groups), start=1):
        X_train = X.iloc[train_idx]
        X_test = X.iloc[test_idx]
        y_reg_train = y_reg.iloc[train_idx]
        y_reg_test = y_reg.iloc[test_idx]
        y_class_train = y_class.iloc[train_idx]
        y_class_test = y_class.iloc[test_idx]

        reg = RandomForestRegressor(
            n_estimators=n_estimators,
            random_state=RANDOM_STATE + fold,
            n_jobs=-1,
            min_samples_leaf=2,
        )
        clf = RandomForestClassifier(
            n_estimators=n_estimators,
            random_state=RANDOM_STATE + fold,
            n_jobs=-1,
            min_samples_leaf=2,
            class_weight="balanced",
        )
        reg.fit(X_train, y_reg_train)
        clf.fit(X_train, y_class_train)
        reg_pred = reg.predict(X_test)
        class_pred = clf.predict(X_test)
        class_prob = clf.predict_proba(X_test)[:, 1]
        rows.append(
            {
                "fold": fold,
                "reg_mae": mean_absolute_error(y_reg_test, reg_pred),
                "reg_rmse": float(np.sqrt(mean_squared_error(y_reg_test, reg_pred))),
                "reg_r2": r2_score(y_reg_test, reg_pred),
                "clf_accuracy": accuracy_score(y_class_test, class_pred),
                "clf_precision": precision_score(y_class_test, class_pred, zero_division=0),
                "clf_recall": recall_score(y_class_test, class_pred, zero_division=0),
                "clf_f1": f1_score(y_class_test, class_pred, zero_division=0),
                "clf_roc_auc": roc_auc_score(y_class_test, class_prob),
                "clf_average_precision": average_precision_score(y_class_test, class_prob),
            }
        )
    return pd.DataFrame(rows)


def summarize_cv(cv_df):
    rows = []
    for col in cv_df.columns:
        if col == "fold":
            continue
        rows.append({"metric": col, "mean": cv_df[col].mean(), "std": cv_df[col].std(ddof=0)})
    return pd.DataFrame(rows)


def train_final_models(X_train, y_reg_train, y_class_train, n_estimators):
    models = {
        "random_forest_regressor": RandomForestRegressor(
            n_estimators=n_estimators,
            random_state=RANDOM_STATE,
            n_jobs=-1,
            min_samples_leaf=2,
        ),
        "hist_gradient_boosting_regressor": HistGradientBoostingRegressor(
            random_state=RANDOM_STATE,
            max_iter=250,
            learning_rate=0.05,
            l2_regularization=0.01,
        ),
        "random_forest_classifier": RandomForestClassifier(
            n_estimators=n_estimators,
            random_state=RANDOM_STATE,
            n_jobs=-1,
            min_samples_leaf=2,
            class_weight="balanced",
        ),
        "hist_gradient_boosting_classifier": HistGradientBoostingClassifier(
            random_state=RANDOM_STATE,
            max_iter=250,
            learning_rate=0.05,
            l2_regularization=0.01,
        ),
    }
    models["random_forest_regressor"].fit(X_train, y_reg_train)
    models["hist_gradient_boosting_regressor"].fit(X_train, y_reg_train)
    models["random_forest_classifier"].fit(X_train, y_class_train)
    models["hist_gradient_boosting_classifier"].fit(X_train, y_class_train)
    return models


def holdout_metrics(models, X_test, y_reg_test, y_class_test):
    rows = []
    for name, model in models.items():
        if "regressor" in name:
            pred = model.predict(X_test)
            rows.append(
                {
                    "model": name,
                    "target": "ddg_regression",
                    "mae": mean_absolute_error(y_reg_test, pred),
                    "rmse": float(np.sqrt(mean_squared_error(y_reg_test, pred))),
                    "r2": r2_score(y_reg_test, pred),
                }
            )
        else:
            pred = model.predict(X_test)
            prob = model.predict_proba(X_test)[:, 1]
            rows.append(
                {
                    "model": name,
                    "target": "destabilizing_classification",
                    "accuracy": accuracy_score(y_class_test, pred),
                    "precision": precision_score(y_class_test, pred, zero_division=0),
                    "recall": recall_score(y_class_test, pred, zero_division=0),
                    "f1": f1_score(y_class_test, pred, zero_division=0),
                    "roc_auc": roc_auc_score(y_class_test, prob),
                    "average_precision": average_precision_score(y_class_test, prob),
                }
            )
    return pd.DataFrame(rows)


def feature_importance(model, X_test, y_test, output_csv):
    result = permutation_importance(
        model,
        X_test,
        y_test,
        n_repeats=10,
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


def explanation_for_row(row):
    reasons = []
    if pd.to_numeric(row.get("conservation_wt_freq"), errors="coerce") >= 0.8:
        reasons.append("high wild-type conservation")
    if pd.to_numeric(row.get("af2_is_buried_proxy"), errors="coerce") == 1:
        reasons.append("buried residue environment")
    if pd.to_numeric(row.get("af2_neighbor_count_8A"), errors="coerce") >= 12:
        reasons.append("dense local packing")
    if abs(pd.to_numeric(row.get("delta_volume"), errors="coerce")) >= 80:
        reasons.append("large side-chain volume change")
    if abs(pd.to_numeric(row.get("delta_hydropathy"), errors="coerce")) >= 3:
        reasons.append("large hydropathy shift")
    if pd.to_numeric(row.get("grantham"), errors="coerce") >= 120:
        reasons.append("severe biochemical substitution")
    if pd.to_numeric(row.get("af2_plddt"), errors="coerce") < 70:
        reasons.append("low-confidence AlphaFold region, interpret structurally with caution")
    if not reasons:
        reasons.append("moderate feature changes relative to other PRNP mutations")
    return "Predicted destabilizing due to " + ", ".join(reasons[:4]) + "."


def add_predictions(df, labeled, X_all, X_labeled, models):
    out = df.copy()
    shared = X_all.reindex(columns=X_labeled.columns, fill_value=0)
    out["predicted_ddg_rf"] = models["random_forest_regressor"].predict(shared)
    out["predicted_ddg_hgb"] = models["hist_gradient_boosting_regressor"].predict(shared)
    out["predicted_destabilizing_probability_rf"] = models["random_forest_classifier"].predict_proba(shared)[:, 1]
    out["predicted_destabilizing_probability_hgb"] = models["hist_gradient_boosting_classifier"].predict_proba(shared)[:, 1]
    out["ddg_explanation"] = out.apply(explanation_for_row, axis=1)
    return out


def hotspot_outputs(pred_df, output_dir):
    agg = pred_df.groupby("position", as_index=False).agg(
        wt_aa=("wt_aa", "first"),
        observed_ddg_mean=("ddg", "mean"),
        observed_ddg_max=("ddg", "max"),
        observed_destabilizing_fraction=("destabilizing_ddg_label", "mean"),
        predicted_ddg_mean=("predicted_ddg_rf", "mean"),
        predicted_ddg_max=("predicted_ddg_rf", "max"),
        predicted_destabilizing_probability_mean=("predicted_destabilizing_probability_rf", "mean"),
        af2_plddt=("af2_plddt", "first"),
        af2_neighbor_count_8A=("af2_neighbor_count_8A", "first"),
        af2_is_buried_proxy=("af2_is_buried_proxy", "first"),
        is_globular_domain_region=("is_globular_domain_region", "first"),
    )
    agg["hotspot_score"] = (
        agg["predicted_ddg_mean"].rank(pct=True)
        + agg["predicted_destabilizing_probability_mean"].rank(pct=True)
        + agg["af2_neighbor_count_8A"].rank(pct=True)
    ) / 3.0
    agg = agg.sort_values("hotspot_score", ascending=False)
    agg.to_csv(output_dir / "residue_hotspot_rankings.csv", index=False)

    heatmap_ddg = pred_df.pivot_table(index="position", columns="mut_aa", values="ddg", aggfunc="mean")
    heatmap_pred = pred_df.pivot_table(index="position", columns="mut_aa", values="predicted_ddg_rf", aggfunc="mean")
    heatmap_prob = pred_df.pivot_table(index="position", columns="mut_aa", values="predicted_destabilizing_probability_rf", aggfunc="mean")
    heatmap_ddg.to_csv(output_dir / "heatmap_observed_ddg.csv")
    heatmap_pred.to_csv(output_dir / "heatmap_predicted_ddg.csv")
    heatmap_prob.to_csv(output_dir / "heatmap_destabilizing_probability.csv")

    structural_map = agg[
        [
            "position",
            "wt_aa",
            "hotspot_score",
            "predicted_ddg_mean",
            "predicted_destabilizing_probability_mean",
            "af2_plddt",
            "af2_neighbor_count_8A",
            "af2_is_buried_proxy",
            "is_globular_domain_region",
        ]
    ].copy()
    structural_map.to_csv(output_dir / "structural_sensitivity_map.csv", index=False)
    return agg


def maybe_plot_heatmap(csv_path, png_path, title):
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return False
    data = pd.read_csv(csv_path, index_col=0)
    plt.figure(figsize=(12, 8))
    plt.imshow(data.values, aspect="auto", interpolation="nearest")
    plt.colorbar(label=title)
    plt.xlabel("Mutant amino acid")
    plt.ylabel("PRNP position")
    plt.xticks(range(len(data.columns)), data.columns)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(png_path, dpi=220)
    plt.close()
    return True


def main():
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    full_df, labeled = load_labeled_data(args.input_csv, args.ddg_destabilizing_threshold)
    X_labeled = numeric_feature_matrix(labeled)
    X_all = numeric_feature_matrix(full_df)
    y_reg = labeled["ddg"]
    y_class = labeled["destabilizing_ddg_label"]
    groups = labeled["position"]

    train_idx, test_idx = grouped_holdout(X_labeled, y_reg, groups)
    models = train_final_models(
        X_labeled.iloc[train_idx],
        y_reg.iloc[train_idx],
        y_class.iloc[train_idx],
        args.n_estimators,
    )
    cv_df = cv_metrics(X_labeled, y_reg, y_class, groups, args.folds, max(150, args.n_estimators // 2))
    holdout_df = holdout_metrics(models, X_labeled.iloc[test_idx], y_reg.iloc[test_idx], y_class.iloc[test_idx])

    cv_df.to_csv(output_dir / "ddg_cross_validation_fold_metrics.csv", index=False)
    summarize_cv(cv_df).to_csv(output_dir / "ddg_cross_validation_summary.csv", index=False)
    holdout_df.to_csv(output_dir / "ddg_holdout_metrics.csv", index=False)
    feature_importance(
        models["random_forest_regressor"],
        X_labeled.iloc[test_idx],
        y_reg.iloc[test_idx],
        output_dir / "feature_importance_ddg_regression.csv",
    )
    feature_importance(
        models["random_forest_classifier"],
        X_labeled.iloc[test_idx],
        y_class.iloc[test_idx],
        output_dir / "feature_importance_destabilizing_classification.csv",
    )

    pred_df = add_predictions(full_df, labeled, X_all, X_labeled, models)
    pred_df.to_csv(output_dir / "mutation_ddg_predictions_with_explanations.csv", index=False)
    hotspot_outputs(pred_df, output_dir)
    maybe_plot_heatmap(output_dir / "heatmap_predicted_ddg.csv", output_dir / "heatmap_predicted_ddg.png", "Predicted DDG")
    maybe_plot_heatmap(output_dir / "heatmap_destabilizing_probability.csv", output_dir / "heatmap_destabilizing_probability.png", "Destabilizing Probability")

    for name, model in models.items():
        joblib.dump(model, output_dir / f"{name}.joblib")

    manifest = {
        "input_csv": args.input_csv,
        "labeled_rows": int(len(labeled)),
        "features": int(X_labeled.shape[1]),
        "grouped_by": "position",
        "destabilizing_threshold_ddg": args.ddg_destabilizing_threshold,
        "excluded_feature_prefixes": ["foldx_", "ddg"],
        "note": "FoldX DDG labels are computational labels generated on an AlphaFold structure, not experimental measurements.",
    }
    (output_dir / "ddg_model_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    print("Labeled rows:", len(labeled))
    print("Features:", X_labeled.shape[1])
    print("Outputs:", output_dir)


if __name__ == "__main__":
    main()
