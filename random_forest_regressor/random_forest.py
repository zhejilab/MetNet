import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import itertools
import scipy.stats
import os
import sys

from sklearn.tree import DecisionTreeClassifier
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import accuracy_score, confusion_matrix, classification_report, roc_curve, auc, mean_squared_error
from sklearn.model_selection import cross_val_score

out_path = sys.argv[3]

X = pd.read_csv(sys.argv[1], index_col = 0)
y = pd.read_csv(sys.argv[2], index_col = 0)

rand_forest = RandomForestRegressor(random_state=12)
hp_grid_rf = {'n_estimators': list(range(1,10))+list(range(10,501,5)),
              'max_depth': [2, 3, 4]+list(range(5,51,5))+[None],
              'max_features': [None, 'sqrt', 'log2']}

clf_rf = GridSearchCV(rand_forest, hp_grid_rf, cv=3, scoring='neg_mean_squared_error', n_jobs = 5)

clf_rf.fit(X, y)
cv_rf_results = clf_rf.cv_results_
df_rf_results = pd.DataFrame(cv_rf_results)
df_rf_results.to_csv(out_path+'/gs_hm_df.csv')

# Extract the best hyperparameter values
best_hyperparam_rf = clf_rf.best_params_
pd.DataFrame(best_hyperparam_rf, index=[0]).to_csv(out_path+'/best_parameters.csv')

best_n_est_rf = best_hyperparam_rf['n_estimators']
best_max_depth_rf = best_hyperparam_rf['max_depth']
best_max_features = best_hyperparam_rf['max_features']

rf_best = RandomForestRegressor(n_estimators=best_n_est_rf, max_depth=best_max_depth_rf, max_features = best_max_features, random_state=12)

# Train the final models
rf_best.fit(X, y)
y_pred_rf = rf_best.predict(X)

best_model_acc_cm = pd.DataFrame([mean_squared_error(y, y_pred_rf), scipy.stats.pearsonr(y, y_pred_rf)[0], scipy.stats.pearsonr(y, y_pred_rf)[1], scipy.stats.spearmanr(y, y_pred_rf)[0], scipy.stats.spearmanr(y, y_pred_rf)[1]])
best_model_acc_cm.index = ['Mean_squared_error', "Pearson_r", "Pearson_p", "Spearman_r", "Spearman_p"]
best_model_acc_cm.to_csv(out_path+'/best_model_results.csv')

gini_temp_df = pd.DataFrame(rf_best.feature_importances_)
gini_temp_df.index = X.columns
gini_temp_df.to_csv(out_path+'/gini_index.csv')

pd.DataFrame([X.index, y, y_pred_rf]).T.to_csv(out_path+'/predicted_values.csv')

