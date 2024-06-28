# Example for accurate predictions/reactants

import pandas as pd
import matplotlib.pyplot as plt
from pywaffle import Waffle

important_feats_data = pd.read_excel(
    r'High_impact_features.xlsx',
    sheet_name=None)
important_feats_data = important_feats_data
rxn_center_data = pd.read_excel(
    r'Rxn_center_matches_to_halides.xlsx',
    sheet_name=None)
rxn_center_index_data = pd.read_excel(
    r'Rxn_center_token_indices.xlsx',
    sheet_name=None)

halide_tokens = ['I', 'Br', 'Cl', 'F']
halide_dict = {'I': '[I-]', 'Br': '[Br-]', 'Cl': '[Cl-]', 'F': '[F-]'}


def _get_prediction_contribution_sign(sample, token_index):
    sample_important_feats_data = important_feats_data[f'Total test index {sample}']

    for i, row in sample_important_feats_data.iterrows():
        if row['Token index'] == token_index:
            importance = row['Mean IG']
            sem = row['Sem']

    if importance > 0 and (importance - sem) > 0:
        res = '+'
    elif importance < 0 and (importance + sem) < 0:
        res = '-'
    else:
        res = 'o'

    return res


def get_prediction_contribution_sign(token_status_samples, rxn_center_index_data):
    positive = []
    negative = []
    both = []

    for sample in token_status_samples:
        sample_important_feats_data = important_feats_data[f'Total test index {sample}']
        sample_important_feats = sample_important_feats_data['Token index'].to_list()
      
        for i, row in rxn_center_index_data.iterrows():
            if row['Total test index'] == sample:
                token_index = row['LG atom token indices']

        if token_index in sample_important_feats:
            sign = _get_prediction_contribution_sign(sample, token_index)
        else:
            sign = 'o'

        if sign == '+':
            positive.append(sample)
        if sign == '-':
            negative.append(sample)
        if sign == 'o':
            both.append(sample)

    return positive, negative, both


def generate_percentages(percentile):
    percentile_rxn_center_data = rxn_center_data[f'{percentile}']
    percentile_rxn_center_index_data = rxn_center_index_data[f'{percentile}']

    num_pos_each_frag = []
    num_neg_each_frag = []
    num_both_each_frag = []
    padding_each_frag = []

    for token in halide_tokens:
        lg_match = []

        for i, row in percentile_rxn_center_data.iterrows():
            if row[halide_dict[token]] == 'LG match':
                lg_match.append(row['Total test index'])

        prediction_contribution_signs_lg = get_prediction_contribution_sign(lg_match, 
                                                                            percentile_rxn_center_index_data)
        num_pos = len(prediction_contribution_signs_lg[0])
        num_neg = len(prediction_contribution_signs_lg[1])
        num_both = len(prediction_contribution_signs_lg[2])
        total = num_pos + num_neg + num_both
        padding = 33 - total

        num_pos_each_frag.append(num_pos)
        num_neg_each_frag.append(num_neg)
        num_both_each_frag.append(num_both)
        padding_each_frag.append(padding)

        print(token)
        if total > 0:
            print(f'% pos: {round((num_pos / total) * 100, 0)}')
            print(f'% neg: {round((num_neg / total) * 100, 0)}')
            print(f'% LI: {round((num_both / total) * 100, 0)}')

        print('\n')

    plt.figure(
        figsize=[10, 2],
        FigureClass=Waffle,
        plots={
            411: {'values': [num_pos_each_frag[0], num_neg_each_frag[0], num_both_each_frag[0], padding_each_frag[0]]},
            412: {'values': [num_pos_each_frag[1], num_neg_each_frag[1], num_both_each_frag[1], padding_each_frag[1]]},
            413: {'values': [num_pos_each_frag[2], num_neg_each_frag[2], num_both_each_frag[2], padding_each_frag[2]]},
            414: {'values': [num_pos_each_frag[3], num_neg_each_frag[3], num_both_each_frag[3], padding_each_frag[3]]}},
        rows=1,
        colors=['crimson', 'lightblue', 'grey', 'white'],
        characters='⬤',
        font_size=20,
    )

    plt.savefig(
        f'BERT_{percentile}_reacts_LG_effects.pdf',
        dpi=1000, bbox_inches='tight')


generate_percentages('Accurate predictions')