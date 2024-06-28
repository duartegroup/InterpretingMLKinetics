import pandas as pd

common_feats_data = pd.read_excel(
    r'High_impact_feature_counts.xlsx',
    sheet_name=None)
percentile_data = pd.read_excel(
    r'Accurate_and_inaccurate_predictions.xlsx',
    sheet_name=None)


def get_percentages(percentile):
    percentile_indices = percentile_data[percentile]['Total test index'].to_list()
    num_percentile_indices = len(percentile_indices)
    percentile_common_feats_data = common_feats_data[percentile]
    percentile_common_feats = percentile_common_feats_data['Feature'].to_list()
    times_in_high_impact = percentile_common_feats_data['Times high impact'].to_list()
    common_feat_dict = {i: j for i, j in zip(percentile_common_feats, times_in_high_impact)}
    temp_times_in_high_impact = common_feat_dict['Reciprocal temp \K-1']
    perc_in_high_impact = (temp_times_in_high_impact / num_percentile_indices) * 100
    print(perc_in_high_impact)


get_percentages('Accurate predictions')
get_percentages('Inaccurate predictions')