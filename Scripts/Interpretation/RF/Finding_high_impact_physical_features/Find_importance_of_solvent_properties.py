import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

common_feats_data = pd.read_excel(
    r'High_impact_feature_counts.xlsx',
    sheet_name=None)
percentile_data = pd.read_excel(
    r'Accurate_and_inaccurate_predictions.xlsx',
    sheet_name=None)
label_data = pd.read_excel(
    r'Labels_for_all_continuous_features.xlsx',
    sheet_name=None)

label_dict = {}
for sheet_name in label_data.keys():
    sheet = label_data[sheet_name]
    properties = sheet['Property'].to_list()
    labels = sheet['Label'].to_list()
    for prop, lab in zip(properties, labels):
        label_dict[prop] = lab


def get_percentages(percentile):
    percentile_indices = percentile_data[percentile]['Total test index'].to_list()
    num_percentile_indices = len(percentile_indices)
    percentile_common_feats_data = common_feats_data[percentile]
    percentile_common_feats = percentile_common_feats_data['Feature'].to_list()
    times_in_high_impact = percentile_common_feats_data['Times high impact'].to_list()
    common_feat_dict = {i: j for i, j in zip(percentile_common_feats, times_in_high_impact)}
    solvent_feats = [item for item in percentile_common_feats if 'Solvent' in item and 'ratio' not in item]
    solvent_feats_times_in_high_impact = [common_feat_dict[i] for i in solvent_feats]
    solvent_feats_perc_in_high_impact = [(i / num_percentile_indices) * 100 for i in solvent_feats_times_in_high_impact]

    print(solvent_feats_perc_in_high_impact)

    fig, ax = plt.subplots(figsize=(14, 6))
    xpos = range(len(solvent_feats))

    colors = []
    axis_labels = []

    for item in solvent_feats:
        if 'Solvent 1' in item:
            colors.append('orange')
            axis_labels.append(label_dict[item].split('Solvent 1 ')[-1])
        else:
            colors.append('yellow')
            axis_labels.append(label_dict[item].split('Solvent 2 ')[-1])

    plt.bar(xpos, solvent_feats_perc_in_high_impact, color=colors, alpha=0.6, width=0.5, edgecolor='black')

    plt.xlabel('Property', fontsize=20, labelpad=20)
    plt.ylabel('% Property is high impact', fontsize=20, labelpad=10)

    print(len(axis_labels))

    plt.xticks(xpos, axis_labels, fontsize=20, rotation=90)
    plt.yticks(fontsize=20)
    plt.tight_layout()
    plt.ylim(0, 110)

    if percentile == 'Accurate predictions':
        custom_lines = [Patch(facecolor='orange', edgecolor='black', alpha=0.6,
                              label='Solvent 1'),
                        Patch(facecolor='yellow', edgecolor='black', alpha=0.6,
                              label='Solvent 2')]

        ax.legend(custom_lines, ['Solvent 1', 'Solvent 2'], frameon=False, fontsize=20, loc='upper right')
    plt.savefig(f'RF_{percentile}_solvent_importances.pdf',dpi=1000, bbox_inches='tight')


get_percentages('Accurate predictions')
print('\n')
get_percentages('Inaccurate predictions')