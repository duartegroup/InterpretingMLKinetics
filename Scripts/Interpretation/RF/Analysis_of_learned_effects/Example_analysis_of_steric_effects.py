# Example for accurate predictions

import pandas as pd
import matplotlib.pyplot as plt
from pywaffle import Waffle

important_feats_data = pd.read_excel(
    r'High_impact_features.xlsx',
    sheet_name=None)

substruct_match_data = pd.read_excel(
    r'Substructure_matches_to_rxn_center_sterics.xlsx',
    sheet_name=None)
center_category_data = pd.read_excel(
    r'Unsubstituted_and_C_substituted_substrates.xlsx',
    sheet_name=None)

frags = ['C-C', 'C-C-C', 'C-C-C-C']


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
    isida_frag = frag
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
    percentile_center_category_data = center_category_data[f'{percentile}']

    c_substituted = percentile_center_category_data['Total test index (C substituted)'].to_list()
    c_substituted = [int(item) for item in c_substituted if pd.isna(item) == False]

    unsubstituted = percentile_center_category_data['Total test index (unsubstituted)'].to_list()
    unsubstituted = [int(item) for item in unsubstituted if pd.isna(item) == False]

    num_pos_each_frag_c_substituted = []
    num_neg_each_frag_c_substituted = []
    num_both_each_frag_c_substituted = []
    padding_each_frag_c_substituted = []

    num_pos_each_frag_unsubstituted = []
    num_neg_each_frag_unsubstituted = []
    num_both_each_frag_unsubstituted = []
    padding_each_frag_unsubstituted = []

    for frag in frags:
        electrophilic_match = []
        for sample in c_substituted:
            sample_substruct_match_data = substruct_match_data[f'Total test index {sample}']
            sample_substruct_match_columns = [item for item in sample_substruct_match_data.columns if
                                              'count' not in item]
            frag_substruct_match_columns = [item for item in sample_substruct_match_columns if frag in item.split(' ')]
            rxn_center_matches = [item for item in frag_substruct_match_columns if
                                  sample_substruct_match_data[item].to_list()[0] == 'Rxn center match']
            if rxn_center_matches == [f'C atom in {frag}', f'C atom (major prod) in {frag}']:
                electrophilic_match.append(sample)

        unsubstituted_no_nu_match = []
        for sample in unsubstituted:
            sample_substruct_match_data = substruct_match_data[f'Total test index {sample}']
            sample_substruct_match_columns = [item for item in sample_substruct_match_data.columns if
                                              'count' not in item]
            frag_substruct_match_columns = [item for item in sample_substruct_match_columns if frag in item.split(' ')]
            rxn_center_matches = [item for item in frag_substruct_match_columns if
                                  sample_substruct_match_data[item].to_list()[0] == 'Rxn center match']
            if rxn_center_matches == []:
                unsubstituted_no_nu_match.append(sample)

        prediction_contribution_signs_c_substituted = get_prediction_contribution_sign(electrophilic_match, frag)
        prediction_contribution_signs_unsubstituted = get_prediction_contribution_sign(unsubstituted_no_nu_match, frag)

        num_pos_c_substituted = len(prediction_contribution_signs_c_substituted[0])
        num_neg_c_substituted = len(prediction_contribution_signs_c_substituted[1])
        num_both_c_substituted = len(prediction_contribution_signs_c_substituted[2])
        total_c_substituted = num_pos_c_substituted + num_neg_c_substituted + num_both_c_substituted
        padding_c_substituted = 33 - total_c_substituted

        num_pos_each_frag_c_substituted.append(num_pos_c_substituted)
        num_neg_each_frag_c_substituted.append(num_neg_c_substituted)
        num_both_each_frag_c_substituted.append(num_both_c_substituted)
        padding_each_frag_c_substituted.append(padding_c_substituted)

        num_pos_unsubstituted = len(prediction_contribution_signs_unsubstituted[0])
        num_neg_unsubstituted = len(prediction_contribution_signs_unsubstituted[1])
        num_both_unsubstituted = len(prediction_contribution_signs_unsubstituted[2])

        total_unsubstituted = num_pos_unsubstituted + num_neg_unsubstituted + num_both_unsubstituted
        padding_unsubstituted = 33 - total_unsubstituted

        num_pos_each_frag_unsubstituted.append(num_pos_unsubstituted)
        num_neg_each_frag_unsubstituted.append(num_neg_unsubstituted)
        num_both_each_frag_unsubstituted.append(num_both_unsubstituted)
        padding_each_frag_unsubstituted.append(padding_unsubstituted)

        print(frag, '\n')

        print('substituted', '\n')

        if total_c_substituted > 0:
            print(f'% pos: {round((num_pos_c_substituted / total_c_substituted) * 100, 0)}')
            print(f'% neg: {round((num_neg_c_substituted / total_c_substituted) * 100, 0)}')
            print(f'% LI: {round((num_both_c_substituted / total_c_substituted) * 100, 0)}')

        print('\n')

        print('unsubstituted')

        if total_unsubstituted > 0:
            print(f'% pos: {round((num_pos_unsubstituted / total_unsubstituted) * 100, 0)}')
            print(f'% neg: {round((num_neg_unsubstituted / total_unsubstituted) * 100, 0)}')
            print(f'% LI: {round((num_both_unsubstituted / total_unsubstituted) * 100, 0)}')

        print('\n')

    plt.figure(figsize=[10, 3],
               FigureClass=Waffle,
               plots={611: {'values': [num_pos_each_frag_c_substituted[0], num_neg_each_frag_c_substituted[0],
                                       num_both_each_frag_c_substituted[0], padding_each_frag_c_substituted[0]]},
                      614: {'values': [num_pos_each_frag_unsubstituted[0], num_neg_each_frag_unsubstituted[0],
                                       num_both_each_frag_unsubstituted[0], padding_each_frag_unsubstituted[0]]},
                      612: {'values': [num_pos_each_frag_c_substituted[1], num_neg_each_frag_c_substituted[1],
                                       num_both_each_frag_c_substituted[1], padding_each_frag_c_substituted[1]]},
                      615: {'values': [num_pos_each_frag_unsubstituted[1], num_neg_each_frag_unsubstituted[1],
                                       num_both_each_frag_unsubstituted[1], padding_each_frag_unsubstituted[1]]},
                      613: {'values': [num_pos_each_frag_c_substituted[2], num_neg_each_frag_c_substituted[2],
                                       num_both_each_frag_c_substituted[2], padding_each_frag_c_substituted[2]]},
                      616: {'values': [num_pos_each_frag_unsubstituted[2], num_neg_each_frag_unsubstituted[2],
                                       num_both_each_frag_unsubstituted[2], padding_each_frag_unsubstituted[2]]}},
               rows=1,
               colors=['crimson', 'lightblue', 'grey', 'white'],
               characters='⬤',
               font_size=20,
               )

    plt.savefig(
        rf'RF_{percentile}_steric_effects.pdf',
        dpi=1000, bbox_inches='tight')


generate_plot('Accurate predictions')