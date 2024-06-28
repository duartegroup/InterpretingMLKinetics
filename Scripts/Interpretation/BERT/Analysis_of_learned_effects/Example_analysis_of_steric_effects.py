# Example for accurate predictions/reacts

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pywaffle import Waffle

important_feats_data = pd.read_excel(
    r'High_impact_features.xlsx',
    sheet_name=None)
center_category_data = pd.read_excel(
    fr'Electrophilic_center_types.xlsx',
    sheet_name=None)
reaction_center_atom_data = pd.read_excel(
    r'Rxn_center_token_indices.xlsx',
    sheet_name=None)
percentile_data = pd.read_excel(
    r'Accurate_and_inaccurate_predictions.xlsx',
    sheet_name=None)

c_substituted_rc_token_index_data = pd.read_excel(
    r'C_substituted_rxn_center_token_indices.xlsx',
    sheet_name=None)


def _get_prediction_contribution_sign(sample, token_index):
    sample_important_feats_data = important_feats_data[f'Total test index {sample}']

    if type(token_index) == int:
        for i, row in sample_important_feats_data.iterrows():
            if row['Token index'] == token_index:
                importance = row['Mean IG']
                sem = row['Sem']

    if type(token_index) == list:
        importances = []
        sems = []
        for index in token_index:
            for i, row in sample_important_feats_data.iterrows():
                if row['Token index'] == index:
                    importances.append(row['Mean IG'])
                    sems.append(row['Sem'])
        importance = sum(importances)
        sem = np.sqrt(sum([x ** 2 for x in sems]))

    if importance > 0 and (importance - sem) > 0:
        res = '+'
    elif importance < 0 and (importance + sem) < 0:
        res = '-'
    else:
        res = 'o'

    return res


def _get_prediction_contribution_sign_center(sample, index_data, sample_important_feats):
    for i, row in index_data.iterrows():
        if row['Total test index'] == sample:
            token_index = row['C atom token indices']
            if token_index in sample_important_feats:
                sign = _get_prediction_contribution_sign(sample, token_index)
            else:
                sign = 'o'
    return sign


def _get_prediction_contribution_sign_non_center(sample, index_data, sample_important_feats):
    sample_index_data = index_data[f'Total test index {sample}']
    token_index = sample_index_data['Token indices'].to_list()
    ind = int(len(token_index) / 2)
    token_index = token_index[:ind]

    token_index = [item for item in token_index if item in sample_important_feats]
    if len(token_index) > 0:
        sign = _get_prediction_contribution_sign(sample, token_index)
    else:
        sign = 'o'

    return sign


def get_prediction_contribution_sign(token_status_samples, token_status, index_data):
    positive = []
    negative = []
    both = []

    for sample in token_status_samples:
        sample_important_feats_data = important_feats_data[f'Total test index {sample}']
        sample_important_feats = sample_important_feats_data['Token index'].to_list()

        if token_status == 'center':
            sign = _get_prediction_contribution_sign_center(sample, index_data, sample_important_feats)
        else:
            sign = _get_prediction_contribution_sign_non_center(sample, index_data, sample_important_feats)

        if sign == '+':
            positive.append(sample)
        if sign == '-':
            negative.append(sample)
        if sign == 'o':
            both.append(sample)

    return positive, negative, both


def generate_percentages(percentile):
    percentile_rxn_center_atom_data = reaction_center_atom_data[f'{percentile}']
    percentile_center_category_data = center_category_data[f'{percentile}']

    unsubstituted = percentile_center_category_data['Total test index (unsubstituted)'].to_list()
    unsubstituted = [int(item) for item in unsubstituted if pd.isna(item) == False]
    c_substituted = percentile_center_category_data['Total test index (C substituted)'].to_list()
    c_substituted = [int(item) for item in c_substituted if pd.isna(item) == False]

    prediction_contribution_signs_unsubstituted = get_prediction_contribution_sign(unsubstituted, 'center',
                                                                                   percentile_rxn_center_atom_data)
    prediction_contribution_signs_c_substituted = get_prediction_contribution_sign(c_substituted, 'center',
                                                                                   percentile_rxn_center_atom_data)
    prediction_contribution_signs_substituting = get_prediction_contribution_sign(c_substituted, 'non_center',
                                                                                  c_substituted_rc_token_index_data)

    padding_unsubstituted = 33 - sum([len(item) for item in prediction_contribution_signs_unsubstituted])
    padding_c_substituted = 33 - sum([len(item) for item in prediction_contribution_signs_c_substituted])
    padding_substituting = 33 - sum([len(item) for item in prediction_contribution_signs_substituting])

    total_c_substituted = len(c_substituted)
    total_unsubstituted = len(unsubstituted)

    print('substituted', '\n')

    if total_c_substituted > 0:
        print(f'% pos: {round((len(prediction_contribution_signs_c_substituted[0]) / total_c_substituted) * 100, 0)}')
        print(f'% neg: {round((len(prediction_contribution_signs_c_substituted[1]) / total_c_substituted) * 100, 0)}')
        print(f'% LI: {round((len(prediction_contribution_signs_c_substituted[2]) / total_c_substituted) * 100, 0)}')

    print('\n')

    print('substituting', '\n')

    if total_c_substituted > 0:
        print(f'% pos: {round((len(prediction_contribution_signs_substituting[0]) / total_c_substituted) * 100, 0)}')
        print(f'% neg: {round((len(prediction_contribution_signs_substituting[1]) / total_c_substituted) * 100, 0)}')
        print(f'% LI: {round((len(prediction_contribution_signs_substituting[2]) / total_c_substituted) * 100, 0)}')

    print('\n')

    print('unsubstituted')

    if total_unsubstituted > 0:
        print(f'% pos: {round((len(prediction_contribution_signs_unsubstituted[0]) / total_unsubstituted) * 100, 0)}')
        print(f'% neg: {round((len(prediction_contribution_signs_unsubstituted[1]) / total_unsubstituted) * 100, 0)}')
        print(f'% LI: {round((len(prediction_contribution_signs_unsubstituted[2]) / total_unsubstituted) * 100, 0)}')

    print('\n')

    plt.figure(
        figsize=[10, 1.5],
        FigureClass=Waffle,
        plots={311: {'values': [len(prediction_contribution_signs_c_substituted[0]),
                                len(prediction_contribution_signs_c_substituted[1]),
                                len(prediction_contribution_signs_c_substituted[2]), padding_c_substituted]},
               313: {'values': [len(prediction_contribution_signs_unsubstituted[0]),
                                len(prediction_contribution_signs_unsubstituted[1]),
                                len(prediction_contribution_signs_unsubstituted[2]), padding_unsubstituted]},
               312: {'values': [len(prediction_contribution_signs_substituting[0]),
                                len(prediction_contribution_signs_substituting[1]),
                                len(prediction_contribution_signs_substituting[2]), padding_substituting]}},
        rows=1,
        colors=['crimson', 'lightblue', 'grey', 'white'],
        characters='⬤',
        font_size=20,
    )

    plt.savefig(
        rf'BERT_{percentile}_reacts_steric_effects.pdf',
        dpi=1000, bbox_inches='tight')


generate_percentages('Accurate predictions')