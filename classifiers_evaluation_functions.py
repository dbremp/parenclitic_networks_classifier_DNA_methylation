"""
This script includes the function to evaluate classification performance.
This functions is used in classifiers_comparison_table.py
"""


import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn import metrics


def classifiers_evaluation(X_train: pd.DataFrame, X_test: pd.DataFrame,
                           y_train: pd.DataFrame, y_test: pd.DataFrame, names: list[str],
                           classifiers: list):
    # perform classification using different classifiers and store their scores in a pd DataFrame
    # names is a list of the classification methods names

    classifiers_performance = {}
    for name, clf in zip(names, classifiers):
        clf = make_pipeline(StandardScaler(), clf)
        clf.fit(X_train, y_train)

        eval = {}
        eval["roc_auc"] = metrics.roc_auc_score(clf.predict(X_test), pd.get_dummies(y_test)["yes"])
        eval["accuracy"] = metrics.accuracy_score(clf.predict(X_test), y_test)
        eval["balanced_accuracy"] = metrics.balanced_accuracy_score(clf.predict(X_test), y_test)

        classifiers_performance[name] = eval

    classifiers_performance = pd.DataFrame.from_dict(classifiers_performance)
    return classifiers_performance
