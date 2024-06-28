# Example for accurate predictions

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pywaffle import Waffle
from rdkit import Chem

important_feats_data = pd.read_excel(
    r'High_impact_features.xlsx',
    sheet_name=None)
solvent_category_data = pd.read_excel(
    fr'Solvent_types.xlsx',
    sheet_name=None)
percentile_data = pd.read_excel(
    r'Accurate_and_inaccurate_predictions.xlsx',
    sheet_name=None)
ext_data = pd.read_excel(
    r'Total_test_processed.xlsx')


def _get_prediction_contribution_sign(sample, token_index):
    sample_important_feats_data = important_feats_data[f'Total test index {sample}']

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


def get_prediction_contribution_sign(token_status_samples, index_data):
    positive = []
    negative = []
    both = []
    for sample in token_status_samples:
        sample_important_feats_data = important_feats_data[f'Total test index {sample}']
        sample_important_feats = sample_important_feats_data['Token index'].to_list()

        sample_index_data = index_data[f'Total test index {sample}']
        token_index = sample_index_data['Token indices'].to_list()
        token_index = [item for item in token_index if item in sample_important_feats]

        if len(token_index) > 0:
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


def generate_heat_map(percentile, nu_samples, nu_type):
    percentile_indices = percentile_data[f'{percentile}']['Total test index'].to_list()

    solvent_categories = ['Polar_protic', 'Polar_aprotic', 'Non_polar']

    num_pos_each_frag = []
    num_neg_each_frag = []
    num_both_each_frag = []
    padding_each_frag = []

    for category in solvent_categories:
        index_data = pd.read_excel(
            fr'{category}_solvent_token_indices.xlsx',
            sheet_name=None)

        category_samples = []
        for sheet_name in index_data:
            sample = int(sheet_name.split(' ')[-1])
            if sample in percentile_indices and sample in nu_samples:
                category_samples.append(sample)

        prediction_contribution_signs = get_prediction_contribution_sign(category_samples, index_data)

        num_pos = len(prediction_contribution_signs[0])
        num_neg = len(prediction_contribution_signs[1])
        num_both = len(prediction_contribution_signs[2])
        total = num_pos + num_neg + num_both
        padding = 33 - total

        num_pos_each_frag.append(num_pos)
        num_neg_each_frag.append(num_neg)
        num_both_each_frag.append(num_both)
        padding_each_frag.append(padding)

        print(category)
        if total > 0:
            print(f'% pos: {round((len(prediction_contribution_signs[0]) / total) * 100, 0)}')
            print(f'% neg: {round((len(prediction_contribution_signs[1]) / total) * 100, 0)}')
            print(f'% LI: {round((len(prediction_contribution_signs[2]) / total) * 100, 0)}')

    print('\n')

    plt.figure(
        figsize=[10, 1.5],
        FigureClass=Waffle,
        plots={
            311: {'values': [num_pos_each_frag[0], num_neg_each_frag[0], num_both_each_frag[0], padding_each_frag[0]]},
            312: {'values': [num_pos_each_frag[1], num_neg_each_frag[1], num_both_each_frag[1], padding_each_frag[1]]},
            313: {
                'values': [num_pos_each_frag[2], num_neg_each_frag[2], num_both_each_frag[2], padding_each_frag[2]]}, },
        rows=1,
        colors=['crimson', 'lightblue', 'grey', 'white'],
        characters='⬤',
        font_size=20,
    )

    plt.savefig(
        fr'BERT_{percentile}_solvent_effects_{nu_type}.pdf',
        dpi=1000, bbox_inches='tight')


nus = ext_data['Nucleophile'].to_list()

neutral = []
anionic = []
for i in range(len(nus)):
    nu = nus[i]
    nu_mol = Chem.MolFromSmiles(nu)
    nu_charge = Chem.GetFormalCharge(nu_mol)
    if nu_charge == 0:
        neutral.append(i)
    else:
        anionic.append(i)

generate_heat_map('Accurate predictions', neutral, 'neutral')
print('\n')
generate_heat_map('Accurate predictions', anionic, 'anionic')
