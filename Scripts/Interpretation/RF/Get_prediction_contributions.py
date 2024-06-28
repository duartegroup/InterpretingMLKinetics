from treeinterpreter import treeinterpreter as ti
import pandas as pd
import pickle

model_ids = [15, 0, 1221, 4096, 21]
cv_split_data = pd.read_excel(r'Total_train_cv_splits.xlsx', sheet_name = None)
train_data_with_feature_names = pd.read_excel(r'Total_train_cv_splits_with_feature_names.xlsx', sheet_name = None)

with pd.ExcelWriter(r'Importances_for_each_test_sample.xlsx') as writer:
    for i in range(len(model_ids)):
        feature_names_sheet = train_data_with_feature_names[f'Fold {i}'].iloc[:,1:]
        fold_feature_names = feature_names_sheet.keys()
        est = pickle.load(open(rf'RF_model_total_train_val_index_{model_ids[i]}', 'rb'))
        ext_data_fold = cv_split_data[f'Total test X fold {i}'].iloc[:,1:]
        best_estimator = est.best_estimator_
        prediction, bias, contributions = ti.predict(best_estimator, ext_data_fold)
        for sample in range(len(ext_data_fold)):
            df = pd.DataFrame({'Feature': fold_feature_names, 'Contribution': contributions[sample]}).sort_values(by='Contribution', ascending=False)
            df.to_excel(writer, sheet_name = f'Total test index {sample} fold {i}')