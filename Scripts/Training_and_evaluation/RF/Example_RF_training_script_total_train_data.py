#  Original script from https://github.com/cimm-kzn/CIMtools.git
#  DOI: 10.1080/1062936X.2021.188310
#  Modified in the current work to allow for description of non-aqueous binary solvents and ionic strength

import numpy as np
from scipy.stats import sem
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import KFold
import pickle
from CIMtools.datasets.reactions import _load
import pandas as pd
from CIMtools.model_selection import TransformationOut
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler, FunctionTransformer
from CIMtools.preprocessing import Fragmentor, CGR, EquationTransformer
from CIMtools.preprocessing.solvent_mod import SolventVectorizer
from CIMtools.preprocessing.conditions_container_mod import DictToConditions, ConditionsToDataFrame

train, train_Y = _load('Total_train', 'load_sn2', return_X_y=True)
total_test, total_test_Y = _load('Total_test', 'load_sn2', return_X_y=True)



def extract_meta(x):
    return [y[0].meta for y in x]


features = ColumnTransformer([('temp', EquationTransformer('1/x'), ['temperature']),
                              ('solv_1', SolventVectorizer(), ['solvent.1']),
                              ('solv_2', SolventVectorizer(), ['solvent.2']),
                              ('amount_1', 'passthrough', ['solvent_amount.1']),
                              ('amount_2', 'passthrough', ['solvent_amount.2']),
                              ('ionic_strength', 'passthrough', ['ionic_strength'])])

conditions = Pipeline([('meta', FunctionTransformer(extract_meta)),
                       ('cond', DictToConditions(solvents=('additive.1', 'additive.2'),
                                                 temperature='temperature',
                                                 amounts=('amount.1', 'amount.2'),
                                                 ionic_strength='ionic_strength')),
                       ('desc', ConditionsToDataFrame()),
                       ('final', features),
                       ('scaler', StandardScaler())])

graph = Pipeline([('CGR', CGR()),
                  ('frg', Fragmentor(fragment_type=3, max_length=4, useformalcharge=True, version='2017')),
                  ('scaler', StandardScaler())])  # All descriptors were normalized to zero mean and unit variance.

pp = ColumnTransformer([('cond', conditions, [0]),
                        ('graph', graph, [0])])


def grouper(cgrs, params):
    '''
    Helper function for transformation-out and solvent-out CV

    Parameters
    ----------
    cgrs: list
        The list dataset.
    params: list
        What we want to sort by.
        If these are solvents, then we transfer the general
        name of the solvents in the meta data in square brackets.
        For example: ['additive.1']
    Returns
    -------
    groups: tuple
        Tuple of parameters
    '''
    groups = []
    for cgr in cgrs:
        group = tuple(cgr.meta[param] for param in params)
        groups.append(group)
    return groups


groups = grouper(train, ['additive.1', 'additive.2'])

transform_out_cv = TransformationOut(n_splits=5, n_repeats=1, random_state=1, shuffle=True)

regular_cv = KFold(n_splits=5, random_state=1, shuffle=True)

Y_pred_train, Y_pred_val, Y_pred_total_test = [], [], []
Y_true_train, Y_true_val, Y_true_total_test = [], [], []
x_trains, x_vals, x_total_tests = [], [], []
y_trains, y_vals = [], []

train_RMSEs = []
val_RMSEs = []
total_test_RMSEs = []


for train_index, val_index in transform_out_cv.split(X=train, groups=groups):
    reactions_train = train[train_index]
    reactions_val = train[val_index]
    x_train = pp.fit_transform(reactions_train.reshape(-1, 1))
    x_val = pp.transform(reactions_val.reshape(-1, 1))
    x_total_test = pp.transform(total_test.reshape(-1, 1))


    y_train = train_Y[train_index]
    y_val = train_Y[val_index]

    x_trains.append(x_train)
    x_vals.append(x_val)
    y_trains.append(y_train)
    y_vals.append(y_val)
    x_total_tests.append(x_total_test)


    est = GridSearchCV(RandomForestRegressor(random_state=1, n_estimators=500),
                       {'max_features': [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 'auto', 'log2', None]},
                       cv=regular_cv, verbose=1, scoring='neg_mean_squared_error', n_jobs=-1).fit(x_train, y_train)
    pickle.dump(est, open(
        f'RF_model_total_train_val_index_{val_index[0]}.sav','wb'))

    y_pred_train, y_pred_val, y_pred_total_test = est.predict(x_train), est.predict(x_val), est.predict(x_total_test)

    train_RMSE = np.sqrt(mean_squared_error(y_pred_train, y_train))
    val_RMSE = np.sqrt(mean_squared_error(y_pred_val, y_val))
    total_test_RMSE = np.sqrt(mean_squared_error(y_pred_total_test, total_test_Y))

    train_RMSEs.append(train_RMSE)
    val_RMSEs.append(val_RMSE)
    total_test_RMSEs.append(total_test_RMSE)


av_train_rmse = round(np.mean(train_RMSEs), 1)
train_rmse_sem = round(sem(train_RMSEs), 1)
av_val_rmse = round(np.mean(val_RMSEs), 1)
val_rmse_sem = round(sem(val_RMSEs), 1)
av_total_test_rmse = round(np.mean(total_test_RMSEs), 1)
total_test_rmse_sem = round(sem(total_test_RMSEs), 1)


results_df = pd.DataFrame(
    {'Train': [av_train_rmse, train_rmse_sem], 'Val': [av_val_rmse, val_rmse_sem],
     'Total test': [av_total_test_rmse, total_test_rmse_sem]})
results_df.to_excel(
    r'Total_train_results_df.xlsx')

with pd.ExcelWriter(
        r'Total_train_cv_splits.xlsx') as writer:
    for i in range(5):
        df_train_x = pd.DataFrame(x_trains[i])
        df_val_x = pd.DataFrame(x_vals[i])
        df_total_test_x = pd.DataFrame(x_total_tests[i])


        df_train_y = pd.DataFrame(y_trains[i])
        df_val_y = pd.DataFrame(y_vals[i])

        df_train_x.to_excel(writer, sheet_name=f'Train X fold {i}')
        df_val_x.to_excel(writer, sheet_name=f'Val X fold {i}')
        df_total_test_x.to_excel(writer, sheet_name=f'Total test X fold {i}')


        df_train_y.to_excel(writer, sheet_name=f'Train Y fold {i}')
        df_val_y.to_excel(writer, sheet_name=f'Val Y fold {i}')