# Example for accurate predictions

import pandas as pd
import matplotlib.pyplot as plt
from pywaffle import Waffle

important_feats_data = pd.read_excel(
    r'High_impact_features.xlsx',
    sheet_name=None)
substruct_match_data = pd.read_excel(
    r'Substructure_matches_to_rxn_center_allylic.xlsx',
    sheet_name=None)
percentile_data = pd.read_excel(
    r'Accurate_and_inaccurate_predictions.xlsx',
    sheet_name=None)

frags = ['C-C=C', 'C-C#C', 'C-C:C']
isida_dict = {'C-C=C': 'C-C=C', 'C-C#C': 'C+C-C', 'C-C:C': 'C*C-C'}

adjacent_o = [74]


def _get_prediction_contribution_sign(sample, isida_frag):
    sample_important_feats_data = important_feats_data[f'Total test index {sample}']
    for i, row in sample_important_feats_data.iterrows():
        if row['Feature'] == isida_frag:
            importance = row['Mean importance']
            sem = row['Sem']

            if importance > 0 and (importance - sem) > 0:
                res = '+'
            elif importance < 0 and (importance + sem) < 0:
                res = '-'
            else:
                res = 'o'
    return res


def get_prediction_contribution_sign(fragment_status_samples, frag):
    positive = []
    negative = []
    both = []
    isida_frag = isida_dict[frag]
    for sample in fragment_status_samples:
        sample_important_feats_data = important_feats_data[f'Total test index {sample}']
        sample_important_feats = sample_important_feats_data['Feature'].to_list()
        if isida_frag in sample_important_feats:
            sign = _get_prediction_contribution_sign(sample, isida_frag)
        else:
            sign = 'o'

        if sign == '+':
            positive.append(sample)
        elif sign == '-':
            negative.append(sample)
        elif sign == 'o':
            both.append(sample)

    return positive, negative, both


def generate_plot(percentile):
    percentile_samples = percentile_data[f'{percentile}']['Total test index'].to_list()

    num_pos_each_frag = []
    num_neg_each_frag = []
    num_both_each_frag = []
    padding_each_frag = []

    for frag in frags:
        electrophilic_match = []
        for sample in percentile_samples:
            sample_substruct_match_data = substruct_match_data[f'Total test index {sample}']
            sample_substruct_match_columns = [item for item in sample_substruct_match_data.columns if
                                              'count' not in item]
            frag_substruct_match_columns = [item for item in sample_substruct_match_columns if frag in item.split(' ')]
            rxn_center_matches = [item for item in frag_substruct_match_columns if
                                  sample_substruct_match_data[item].to_list()[0] == 'Rxn center match']
            if rxn_center_matches == [f'C atom in {frag}', f'C atom (major prod) in {frag}']:
                electrophilic_match.append(sample)

        if frag == 'C-C:C':
            electrophilic_match = [item for item in electrophilic_match if item not in adjacent_o]

        prediction_contribution_signs = get_prediction_contribution_sign(electrophilic_match, frag)

        num_pos = len(prediction_contribution_signs[0])
        num_neg = len(prediction_contribution_signs[1])
        num_both = len(prediction_contribution_signs[2])
        total = num_pos + num_neg + num_both
        padding = 33 - total

        num_pos_each_frag.append(num_pos)
        num_neg_each_frag.append(num_neg)
        num_both_each_frag.append(num_both)
        padding_each_frag.append(padding)

        print(frag)
        if total > 0:
            print(f'% pos: {round((num_pos / total) * 100, 0)}')
            print(f'% neg: {round((num_neg / total) * 100, 0)}')
            print(f'% LI: {round((num_both / total) * 100, 0)}')

        print('\n')

    plt.figure(
        figsize=[10, 1.75],
        FigureClass=Waffle,
        plots={
            311: {'values': [num_pos_each_frag[0], num_neg_each_frag[0], num_both_each_frag[0], padding_each_frag[0]]},
            312: {'values': [num_pos_each_frag[1], num_neg_each_frag[1], num_both_each_frag[1], padding_each_frag[1]]},
            313: {'values': [num_pos_each_frag[2], num_neg_each_frag[2], num_both_each_frag[2], padding_each_frag[2]]}},
        rows=1,
        colors=['crimson', 'lightblue', 'grey', 'white'],
        characters='⬤',
        font_size=20,

    )

    plt.savefig(
        rf'RF_{percentile}_allylic_effects.pdf',
        dpi=1000, bbox_inches='tight')

generate_plot('Accurate predictions')