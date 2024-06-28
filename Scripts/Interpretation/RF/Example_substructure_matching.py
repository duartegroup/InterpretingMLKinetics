import pandas as pd
from rdkit import Chem


ext_data = pd.read_excel(
    r'Total_test_processed.xlsx')
index_data = pd.read_excel(
    r'Rxn_center_indices_resonance_correction.xlsx')
smarts_and_labels_data = pd.read_excel(
    r'Labels_and_smarts_for_all_features.xlsx')
high_impact_data = pd.read_excel(
    r'High_impact_features.xlsx',
    sheet_name=None)


def find_rxn_center_frag_matches(smi, index, patt):
    substruct = Chem.MolFromSmarts(patt)
    mol = Chem.MolFromSmiles(smi)
    matches = mol.GetSubstructMatches(substruct)
    index_in_match = []
    for match in matches:
        if index in match:
            index_in_match.append('matched')
    if len(matches) == 0:
        res = 'No matches in molecule'
    elif len(index_in_match) > 0:
        res = 'Rxn center match'
    else:
        res = 'None'
    count = len(index_in_match)
    return res, count


label_dict = {}
smarts_dict = {}
isida_fragments = smarts_and_labels_data['ISIDA fragment'].to_list()
labels = smarts_and_labels_data['Label'].to_list()
smarts = smarts_and_labels_data['SMARTS'].to_list()
for fragment, label in zip(isida_fragments, labels):
    label_dict[fragment] = label
for fragment, smart in zip(isida_fragments, smarts):
    smarts_dict[fragment] = smart

rxn_center_atoms = ['Nu atom', 'LG atom', 'C atom', 'Nu atom (major prod)', 'C atom (major prod)',
                    'LG atom (nucleofuge)']

with pd.ExcelWriter(
        r'Substructure_matches_to_rxn_center.xlsx') as writer:
    for sheet_name in high_impact_data.keys():
        sample = int(sheet_name.split(' ')[-1])
        sample_isida_fragments = high_impact_data[sheet_name]["Feature"].to_list()
        sample_isida_fragments = [item for item in sample_isida_fragments if
                                  "Solvent" not in item and "temp" not in item and "Ionic" not in item]
        sample_patts = [smarts_dict[item] for item in sample_isida_fragments]
        sample_labels = [label_dict[item] for item in sample_isida_fragments]
        sample_res = []
        for i, row in ext_data.iterrows():
            if i == sample:
                index_row = index_data.iloc[i, :]
                for patt in sample_patts:
                    nu_res, nu_count = find_rxn_center_frag_matches(row['Nucleophile'], int(index_row['Nu index']),
                                                                    patt)
                    lg_res, lg_count = find_rxn_center_frag_matches(row['Substrate'], int(index_row['LG index']), patt)
                    c_res, c_count = find_rxn_center_frag_matches(row['Substrate'], int(index_row['C index']), patt)
                    nu_mp_res, nu_mp_count = find_rxn_center_frag_matches(row['Major product'],
                                                                          int(index_row['Nu index (major prod)']), patt)
                    c_mp_res, c_mp_count = find_rxn_center_frag_matches(row['Major product'],
                                                                        int(index_row['C index (major prod)']), patt)
                    lg_nuc_res, lg_nuc_count = find_rxn_center_frag_matches(row['Nucleofuge'],
                                                                            int(index_row['LG index (nucleofuge)']),
                                                                            patt)
                    patt_res = [nu_res, nu_count, lg_res, lg_count, c_res, c_count, nu_mp_res, nu_mp_count, c_mp_res,
                                c_mp_count, lg_nuc_res, lg_nuc_count]
                    sample_res.extend(patt_res)

        cols = []
        for lab in sample_labels:
            for atom in rxn_center_atoms:
                cols.append(f'{atom} in {lab}')
                cols.append(f'{atom} in {lab} count')

        df = pd.DataFrame(sample_res).transpose()
        df.columns = cols
        df.to_excel(writer, sheet_name=f'Total test index {sample}')