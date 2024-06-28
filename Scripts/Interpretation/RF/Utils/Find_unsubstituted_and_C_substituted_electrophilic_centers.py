import pandas as pd
from rdkit import Chem

percentile_data = pd.read_excel(
    r'Accurate_and_inaccurate_predictions.xlsx',
    sheet_name=None)
ext_data = pd.read_excel(
    r'Total_test_processed.xlsx')
index_data = pd.read_excel(
    r'Rxn_center_indices_resonance_correction.xlsx')

substrates = ext_data['Substrate'].to_list()
c_indices = index_data['C index'].to_list()
lg_indices = index_data['LG index'].to_list()

allyl_patt = '[C,c;+0]-[C,c;+0]=[C,c;+0]'
triple_patt = '[C,c;+0]-[C,c;+0]#[C,c;+0]'
aryl_patt = '[c+0]:[c+0]-[C,c;+0]'
carbonyl_patt = '[C,c;+0]-[C,c;+0]=[O+0]'


def _substruct_match(mol, patt, c_index):
    smarts = Chem.MolFromSmarts(patt)
    matches = mol.GetSubstructMatches(smarts)
    if any(c_index in match for match in matches):
        return True
    else:
        return False


def categorize_substrate(percentile):
    percentile_samples = percentile_data[percentile]['Total test index'].to_list()

    unsubstituted = []
    c_substituted = []

    for sample in percentile_samples:
        substrate = substrates[sample]
        substrate_mol = Chem.AddHs(Chem.MolFromSmiles(substrate))
        c_index = c_indices[sample]
        lg_index = lg_indices[sample]
        c_atom = substrate_mol.GetAtomWithIdx(c_index)
        c_neighbor_atoms = c_atom.GetNeighbors()
        c_neighbor_indices = [atom.GetIdx() for atom in c_neighbor_atoms]
        c_neighbor_atoms = [substrate_mol.GetAtomWithIdx(x) for x in c_neighbor_indices if x != lg_index]
        c_neighbor_symbols = [atom.GetSymbol() for atom in c_neighbor_atoms]

        if all(item == 'H' for item in c_neighbor_symbols):
            unsubstituted.append(sample)

        else:
            allyl_match = _substruct_match(substrate_mol, allyl_patt, c_index)
            triple_match = _substruct_match(substrate_mol, triple_patt, c_index)
            aryl_match = _substruct_match(substrate_mol, aryl_patt, c_index)
            carbonyl_match = _substruct_match(substrate_mol, carbonyl_patt, c_index)

            all_res = [allyl_match, aryl_match, carbonyl_match, triple_match]
            if all(item == False for item in all_res):
                if all(item == 'H' or item == 'C' for item in c_neighbor_symbols):
                    c_substituted.append(sample)

    other = [item for item in percentile_samples if item not in unsubstituted + c_substituted]

    unsubstituted_subs = []

    for item in unsubstituted:
        unsubstituted_subs.append(substrates[item])

    for item in list(set(unsubstituted_subs)):
        print(item)
    print('\n')

    c_substituted_subs = []
    for item in c_substituted:
        c_substituted_subs.append(substrates[item])

    for item in list(set(c_substituted_subs)):
        print(item)
    print('\n')

    other_subs = []
    for item in other:
        other_subs.append(substrates[item])

    for item in list(set(other_subs)):
        print(item)
    print('\n')

    df = pd.DataFrame({'Total test index (unsubstituted)': pd.Series(unsubstituted), 'Total test index (C substituted)': pd.Series(c_substituted)})
    return df


print('Accurate', '\n')
top_df = categorize_substrate('Accurate predictions')
print('Inaccurate', '\n')
bottom_df = categorize_substrate('Inaccurate predictions')

with pd.ExcelWriter(
        r'Unsubstituted_and_C_substituted_substrates.xlsx') as writer:
    top_df.to_excel(writer, 'Accurate predictions')
    bottom_df.to_excel(writer, 'Inaccurate predictions')