# Example for RF
import pandas as pd

mean_importance_data = pd.read_excel(
    r'Mean_importances_for_each_test_sample.xlsx',
    sheet_name=None)
percentile_data = pd.read_excel(
    r'Accurate_and_inaccurate_predictions.xlsx',
    sheet_name=None)

all_percentile_indices = []
for percentile in ['Accurate predictions', 'Inaccurate predictions']:
    sheet = percentile_data[f'{percentile}']
    percentile_indices = sheet['Total test index'].to_list()
    all_percentile_indices.extend(percentile_indices)

thresh = 0.001 # greater than 75% of assigned importances

with pd.ExcelWriter(
        r'High_impact_features.xlsx') as writer:
    for sheet_name in mean_importance_data.keys():
        sample = int(sheet_name.split(' ')[-1])
        if sample in all_percentile_indices:
            sheet = mean_importance_data[sheet_name]
            high_impact = [row for i, row in sheet.iterrows() if abs(row['Mean importance']) >= thresh]
            high_impact_df = pd.DataFrame(high_impact)
            high_impact_df.to_excel(writer, sheet_name)
