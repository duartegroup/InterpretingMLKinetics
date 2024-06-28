import pandas as pd
from rdkit import Chem

external_data = pd.read_excel(r'BERT_input_data.xlsx',sheet_name='Total test fold 0')
rxn_center_data = pd.read_excel(r'Rxn_center_token_indices.xlsx',sheet_name=None)

text = external_data['text'].to_list()
reacts = [item.split('>>')[0] for item in text]
substrates = [item.split('.')[1] for item in reacts]

substrate_mols = [Chem.AddHs(Chem.MolFromSmiles(item)) for item in substrates]

allyl_patt = '[C,c;+0]-[C,c;+0]=[C,c;+0]'
triple_patt = '[C,c;+0]-[C,c;+0]#[C,c;+0]'
aryl_patt = '[c+0]:[c+0]-[C,c;+0]'
carbonyl_patt = '[C,c;+0]-[C,c;+0]=[O+0]'

adjacent_o = [74]

def _substruct_match(mol, patt, c_index):
    smarts = Chem.MolFromSmarts(patt)
    matches = mol.GetSubstructMatches(smarts)
    if any(c_index in match for match in matches):
        return True
    else:
        return False


def find_centers(percentile):
    percentile_rxn_center_data = rxn_center_data[f'{percentile}']
    percentile_indices = percentile_rxn_center_data['Total test index'].to_list()
    percentile_c_indices = percentile_rxn_center_data['C atomic index'].to_list()
    percentile_lg_indices = percentile_rxn_center_data['LG atomic index'].to_list()
    c_index_dict = {i: j for i, j in zip(percentile_indices, percentile_c_indices)}
    lg_index_dict = {i: j for i, j in zip(percentile_indices, percentile_lg_indices)}

    unsubstituted = []
    c_substituted = []
    allylic = []
    triple = []
    aryl = []
    carbonyl = []

    for i in range(len(substrate_mols)):
        if i in percentile_indices:
            substrate_mol = substrate_mols[i]
            c_index = c_index_dict[i]
            lg_index = lg_index_dict[i]
            c_atom = substrate_mol.GetAtomWithIdx(c_index)
            c_neighbor_atoms = c_atom.GetNeighbors()
            c_neighbor_indices = [atom.GetIdx() for atom in c_neighbor_atoms]
            c_neighbor_atoms = [substrate_mol.GetAtomWithIdx(x) for x in c_neighbor_indices if x != lg_index]
            c_neighbor_symbols = [atom.GetSymbol() for atom in c_neighbor_atoms]

            if all(item == 'H' for item in c_neighbor_symbols):
                unsubstituted.append(i)
            else:
                allyl_match = _substruct_match(substrate_mol, allyl_patt, c_index)
                triple_match = _substruct_match(substrate_mol, triple_patt, c_index)
                aryl_match = _substruct_match(substrate_mol, aryl_patt, c_index)
                carbonyl_match = _substruct_match(substrate_mol, carbonyl_patt, c_index)
                if allyl_match == True:
                    allylic.append(i)
                if aryl_match == True and i not in adjacent_o:
                    aryl.append(i)
                if carbonyl_match == True:
                    carbonyl.append(i)
                if triple_match == True:
                    triple.append(i)

                all_res = [allyl_match, aryl_match, carbonyl_match, triple_match]
                if all(item == False for item in all_res):
                    if all(item == 'H' or item == 'C' for item in c_neighbor_symbols):
                        c_substituted.append(i)

    accounted_for = unsubstituted + c_substituted + allylic + triple + aryl + carbonyl
    other = [item for item in percentile_indices if item not in accounted_for]

    print('Unsubstituted')
    unsubstituted_subs = []
    for item in unsubstituted:
        unsubstituted_subs.append(substrates[item])

    for item in list(set(unsubstituted_subs)):
        print(item)
    print('\n')

    print('Alkyl')
    c_substituted_subs = []
    for item in c_substituted:
        c_substituted_subs.append(substrates[item])

    for item in list(set(c_substituted_subs)):
        print(item)
    print('\n')

    print('Allylic')
    allylic_subs = []
    for item in allylic:
        allylic_subs.append(substrates[item])

    for item in list(set(allylic_subs)):
        print(item)
    print('\n')

    print('Triple')
    triple_subs = []
    for item in triple:
        triple_subs.append(substrates[item])

    for item in list(set(triple_subs)):
        print(item)
    print('\n')

    print('Aryl')
    aryl_subs = []
    for item in aryl:
        aryl_subs.append(substrates[item])

    for item in list(set(aryl_subs)):
        print(item)
    print('\n')

    print('Carbonyl')
    carbonyl_subs = []
    for item in carbonyl:
        carbonyl_subs.append(substrates[item])
    for item in list(set(carbonyl_subs)):
        print(item)
    print('\n')

    print('Excluded')
    print(substrates[74])
    print('\n')

    print('Other')
    other_subs = []
    for item in other:
        other_subs.append(substrates[item])

    for item in list(set(other_subs)):
        print(item)
    print('\n')

    df = pd.DataFrame({'Total test index (unsubstituted)': pd.Series(unsubstituted), 'Total test index (C substituted)': pd.Series(c_substituted),
                       'Total test index (allylic)': pd.Series(allylic), 'Total test index (triple)': pd.Series(triple), 'Total test index (aryl)': pd.Series(aryl),
                       'Total test index (carbonyl)': pd.Series(carbonyl), 'Total test index (other)': pd.Series(other)})
    return df


top_df = find_centers('Accurate predictions')
bottom_df = find_centers('Inaccurate predictions')

with pd.ExcelWriter(
        rf'Electrophilic_center_types.xlsx') as writer:
    top_df.to_excel(writer, sheet_name='Accurate predictions')
    bottom_df.to_excel(writer, sheet_name='Inaccurate predictions')

