"""Independent re-implementation of the tier-B depression models in scikit-learn.

Reads the design matrices exported by scripts/05_mh_ml.R (identical train/test split),
tunes each learner with stratified 10-fold CV on the training set only, and evaluates once on
the hold-out set. Purpose: check that the R results do not depend on one software stack.
Run from the project root:  .venv/bin/python python/ml_crosscheck.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.inspection import permutation_importance
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import brier_score_loss, roc_auc_score
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

SEED = 20210601
ROOT = Path(__file__).resolve().parents[1]
train = pd.read_csv(ROOT / "data/processed/ml_train_tierB.csv")
test = pd.read_csv(ROOT / "data/processed/ml_test_tierB.csv")
Xtr, ytr = train.drop(columns="dep"), train["dep"]
Xte, yte = test.drop(columns="dep"), test["dep"]
assert list(Xtr.columns) == list(Xte.columns) and "depression" not in Xtr.columns

cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=SEED)
models = {
    "Logistic regression (L2)": GridSearchCV(
        make_pipeline(StandardScaler(), LogisticRegression(max_iter=5000)),
        {"logisticregression__C": [0.01, 0.03, 0.1, 0.3, 1, 3, 10]}, cv=cv, scoring="roc_auc"),
    "Random forest": GridSearchCV(
        RandomForestClassifier(n_estimators=500, random_state=SEED, n_jobs=-1),
        {"max_features": [2, 3, 5, 8], "min_samples_leaf": [5, 15, 30]}, cv=cv, scoring="roc_auc"),
    "Hist. gradient boosting": GridSearchCV(
        HistGradientBoostingClassifier(random_state=SEED, early_stopping=False),
        {"max_depth": [2, 3, 4], "learning_rate": [0.03, 0.1], "max_iter": [100, 300]}, cv=cv, scoring="roc_auc"),
}


def bootstrap_auc_ci(y, p, n_boot=2000, seed=SEED):
    rng = np.random.default_rng(seed)
    y, p = np.asarray(y), np.asarray(p)
    stats = []
    for _ in range(n_boot):
        i = rng.integers(0, len(y), len(y))
        if y[i].min() != y[i].max():
            stats.append(roc_auc_score(y[i], p[i]))
    return np.percentile(stats, [2.5, 97.5])


rows, imps = [], []
for name, gs in models.items():
    gs.fit(Xtr, ytr)
    p = gs.predict_proba(Xte)[:, 1]
    lo, hi = bootstrap_auc_ci(yte, p)
    rows.append({"learner": name, "cv_auc": gs.best_score_, "test_auc": roc_auc_score(yte, p),
                 "lo": lo, "hi": hi, "brier": brier_score_loss(yte, p), "params": json.dumps(gs.best_params_)})
    pi = permutation_importance(gs.best_estimator_, Xte, yte, scoring="roc_auc", n_repeats=30, random_state=SEED)
    imps.append(pd.DataFrame({"learner": name, "feature": Xte.columns, "imp": pi.importances_mean}))
    print(f"{name:28s} CV AUC {gs.best_score_:.3f}  test AUC {rows[-1]['test_auc']:.3f} [{lo:.3f}, {hi:.3f}]")

out = ROOT / "output"
pd.DataFrame(rows).to_csv(out / "ml_results_python.csv", index=False)
pd.concat(imps).to_csv(out / "ml_importance_python.csv", index=False)
