# To generate input data from train isida rxn center counts:
# Drop columns where feature value = 0 for every reaction
# Append temperature and solvent feature values (from regular RF input data)
# Standardize all feature values

import pandas as pd
from rdkit import Chem


train_data = pd.read_excel(
    r"Total_train_processed.xlsx")
index_data = pd.read_excel(
    r'Train_rxn_center_indices_resonance_correction.xlsx')
smarts_and_labels_data = pd.read_excel(
    r'Labels_and_smarts_for_all_features.xlsx')


def find_rxn_center_frag_matches(smi, index, patt):
    smarts = Chem.MolFromSmarts(patt)
    mol = Chem.MolFromSmiles(smi)
    matches = mol.GetSubstructMatches(smarts)
    matches_at_rxn_center = [item for item in matches if index in item]
    return matches_at_rxn_center


label_dict = {}
smarts_dict = {}
isida_fragments = smarts_and_labels_data['ISIDA fragment'].to_list()
labels = smarts_and_labels_data['Label'].to_list()
smarts = smarts_and_labels_data['SMARTS'].to_list()
for fragment, label in zip(isida_fragments, labels):
    label_dict[fragment] = label
for fragment, smart in zip(isida_fragments, smarts):
    smarts_dict[fragment] = smart

all_res = []
for sample, row in train_data.iterrows():
    sample_res = []
    index_row = index_data.iloc[sample, :]
    for isida_frag, patt in zip(isida_fragments, smarts):
        nu_matches = find_rxn_center_frag_matches(row['Nucleophile'], int(index_row['Nu index']), patt)
        lg_matches = find_rxn_center_frag_matches(row['Substrate'], int(index_row['LG index']), patt)
        c_matches = find_rxn_center_frag_matches(row['Substrate'], int(index_row['C index']), patt)
        c_matches = [item for item in c_matches if item not in lg_matches]
        nu_mp_matches = find_rxn_center_frag_matches(row['Major product'], int(index_row['Nu index (major prod)']),
                                                     patt)
        c_mp_matches = find_rxn_center_frag_matches(row['Major product'], int(index_row['C index (major prod)']), patt)
        c_mp_matches = [item for item in c_mp_matches if item not in nu_mp_matches]
        lg_nuc_matches = find_rxn_center_frag_matches(row['Nucleofuge'], int(index_row['LG index (nucleofuge)']), patt)
        rxn_center_matches = nu_matches + lg_matches + c_matches + nu_mp_matches + c_mp_matches + lg_nuc_matches
        rxn_center_count = len(rxn_center_matches)
        sample_res.append(rxn_center_count)
    all_res.append(sample_res)

df = pd.DataFrame(all_res, columns=isida_fragments)

df.to_excel(
    r"Train_isida_rxn_center_counts.xlsx")