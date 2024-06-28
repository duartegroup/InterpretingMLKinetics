import pandas as pd
from rdkit import Chem

external_data = pd.read_excel(
    r'BERT_input_data.xlsx',
    sheet_name='Total test fold 0')
rxn_center_atom_data = pd.read_excel(
    r'Rxn_center_token_indices.xlsx',
    sheet_name=None)
ig_data = pd.read_excel(
    r'IGs_fold_0.xlsx',
    sheet_name=None)
center_category_data = pd.read_excel(
    r'Electrophilic_center_types.xlsx',
    sheet_name=None)

text = external_data['text'].to_list()
reacts = [item.split('>>')[0] for item in text]
substrates = [item.split('.')[1] for item in reacts]
prods = [item.split('>>')[1].split('[RecipTemp]')[0] for item in text]
major_prods = [item.split('.')[0] for item in prods]

substrate_mols = [Chem.AddHs(Chem.MolFromSmiles(item)) for item in substrates]
major_prod_mols = [Chem.AddHs(Chem.MolFromSmiles(item)) for item in major_prods]


def _split_into_molecules(rxn_side_tokens, rxn_side_token_indices):
    delimiter_indices = [i for i, v in enumerate(rxn_side_tokens, 1) if v == '.']
    delimiter_indices = delimiter_indices + [(len(rxn_side_tokens) + 1)]
    molecule_split_tokens = [rxn_side_tokens[i:j - 1] for i, j in zip([0] + delimiter_indices, delimiter_indices)]
    molecule_split_token_indices = [rxn_side_token_indices[i:j - 1] for i, j in
                                    zip([0] + delimiter_indices, delimiter_indices)]

    return molecule_split_tokens, molecule_split_token_indices


def _get_atom_tokens_only(molecule_tokens, molecule_token_indices):
    molecule_tokens_atoms_only = []
    molecule_token_indices_atoms_only = []
    for token, index in zip(molecule_tokens, molecule_token_indices):
        if (token.isalpha() == True and token != 'DOUBLEBOND') or ('[' in token):
            molecule_tokens_atoms_only.append(token)
            molecule_token_indices_atoms_only.append(index)

    return molecule_tokens_atoms_only, molecule_token_indices_atoms_only


def _get_rxn_side_tokens(sample_tokens, side):
    rxn_arrow_token_index = sample_tokens.index('>>')
    recip_temp_index = sample_tokens.index('[RecipTemp]')
    if side == 'reacts':
        rxn_side_token_indices = [x for x in range(len(sample_tokens)) if x < rxn_arrow_token_index and x != 0]
    else:
        rxn_side_token_indices = [x for x in range(len(sample_tokens)) if
                                  x > rxn_arrow_token_index and x < recip_temp_index]

    rxn_side_tokens = [sample_tokens[x] for x in rxn_side_token_indices]

    return rxn_side_tokens, rxn_side_token_indices


def _replace_double_bond(tokens):
    for x in range(len(tokens)):
        if tokens[x] == 0:
            tokens[x] = 'DOUBLEBOND'

    return tokens


def _find_substituting_atom_token_indices(electrophilic_c_index, nu_index_mp, molecule, center_category, species_tokens,
                                          species_token_indices, reacts_or_prods):
    electrophilic_c_atom = molecule.GetAtomWithIdx(electrophilic_c_index)
    substituting_atoms = electrophilic_c_atom.GetNeighbors()
    substituting_c_atoms = [atom for atom in substituting_atoms if atom.GetSymbol() == 'C']

    if center_category == 'C_substituted':
        substituting_atomic_indices = [atom.GetIdx() for atom in substituting_c_atoms]
    if center_category == 'Aryl':
        substituting_aromatic_atoms = [atom for atom in substituting_c_atoms if atom.GetIsAromatic() == True]
        substituting_atomic_indices = [atom.GetIdx() for atom in substituting_aromatic_atoms]

    substituting_atom_token_indices = []
    substituting_atom_tokens = []

    if reacts_or_prods == 'prods':
        for x in substituting_atomic_indices:
            if x != nu_index_mp:
                substituting_atom_token_indices.append(species_token_indices[x])
                substituting_atom_tokens.append(species_tokens[x])

    else:
        for x in substituting_atomic_indices:
            substituting_atom_token_indices.append(species_token_indices[x])
            substituting_atom_tokens.append(species_tokens[x])

    return substituting_atom_tokens, substituting_atom_token_indices


def _find_adjacent_bond_token_indices(center_category, species_tokens, species_token_indices):
    bond_tokens = []
    bond_token_indices = []

    if center_category == 'Allylic':
        toke = 'DOUBLEBOND'
    if center_category == 'Triple':
        toke = '#'

    for x in range(len(species_tokens)):
        if species_tokens[x] == toke:
            bond_tokens.append(species_tokens[x])
            bond_token_indices.append(species_token_indices[x])

    return bond_tokens, bond_token_indices


def find_indices(center_category):
    with pd.ExcelWriter(f'{center_category}_rxn_center_token_indices.xlsx') as writer:
        for percentile in ['Accurate predictions', 'Inaccurate predictions']:
            percentile_center_category_data = center_category_data[f'{percentile}']
            percentile_rxn_center_atom_data = rxn_center_atom_data[f'{percentile}']
            percentile_samples = percentile_rxn_center_atom_data['Total test index'].to_list()

            percentile_electrophilic_c_indices_reacts = percentile_rxn_center_atom_data['C atomic index'].to_list()
            c_index_dict_reacts = {i: j for i, j in zip(percentile_samples, percentile_electrophilic_c_indices_reacts)}

            percentile_electrophilic_c_indices_prods = percentile_rxn_center_atom_data[
                'C atomic index (major prod)'].to_list()
            c_index_dict_prods = {i: j for i, j in zip(percentile_samples, percentile_electrophilic_c_indices_prods)}

            percentile_nu_indices_mp = percentile_rxn_center_atom_data['Nu atomic index (major prod)'].to_list()
            nu_index_mp_dict = {i: j for i, j in zip(percentile_samples, percentile_nu_indices_mp)}

            if center_category == 'C_substituted':
                id = center_category.replace('_',' ')
            else:
                id = center_category.lower()
            center_category_samples = percentile_center_category_data[f'Total test index ({id})'].to_list()
            center_category_samples = [int(item) for item in center_category_samples if pd.isna(item) == False]

            for sample in center_category_samples:
                sample_ig_data = ig_data[f'Total test index {sample}']
                sample_tokens = sample_ig_data['Token'].to_list()
                sample_tokens = _replace_double_bond(sample_tokens)
                reactant_tokens, reactant_token_indices = _get_rxn_side_tokens(sample_tokens, 'reacts')
                product_tokens, product_token_indices = _get_rxn_side_tokens(sample_tokens, 'prods')

                reactant_molecule_tokens, reactant_molecule_token_indices = _split_into_molecules(reactant_tokens,
                                                                                                  reactant_token_indices)
                substrate_tokens, substrate_token_indices = reactant_molecule_tokens[1], \
                                                            reactant_molecule_token_indices[1]
                substrate_tokens_atoms_only, substrate_token_indices_atoms_only = _get_atom_tokens_only(
                    substrate_tokens, substrate_token_indices)

                product_molecule_tokens, product_molecule_token_indices = _split_into_molecules(product_tokens,
                                                                                                product_token_indices)
                major_prod_tokens, major_prod_token_indices = product_molecule_tokens[0], \
                                                              product_molecule_token_indices[0]
                major_prod_tokens_atoms_only, major_prod_token_indices_atoms_only = _get_atom_tokens_only(
                    major_prod_tokens, major_prod_token_indices)

                if center_category == 'C_substituted' or center_category == 'Aryl':
                    electrophilic_c_index_reacts = c_index_dict_reacts[sample]
                    electrophilic_c_index_prods = c_index_dict_prods[sample]
                    nu_index_mp = nu_index_mp_dict[sample]
                    substrate = substrate_mols[sample]
                    major_prod = major_prod_mols[sample]
                    tokens_reacts, token_indices_reacts = _find_substituting_atom_token_indices(
                        electrophilic_c_index_reacts, nu_index_mp, substrate, center_category,
                        substrate_tokens_atoms_only, substrate_token_indices_atoms_only, 'reacts')
                    tokens_prods, token_indices_prods = _find_substituting_atom_token_indices(
                        electrophilic_c_index_prods, nu_index_mp, major_prod, center_category,
                        major_prod_tokens_atoms_only, major_prod_token_indices_atoms_only, 'prods')
                    tokens = tokens_reacts + tokens_prods
                    token_indices = token_indices_reacts + token_indices_prods

                elif center_category == 'Allylic' or center_category == 'Triple':
                    tokens_reacts, token_indices_reacts = _find_adjacent_bond_token_indices(center_category,
                                                                                            substrate_tokens,
                                                                                            substrate_token_indices)
                    tokens_prods, token_indices_prods = _find_adjacent_bond_token_indices(center_category,
                                                                                          major_prod_tokens,
                                                                                          major_prod_token_indices)
                    tokens = tokens_reacts + tokens_prods
                    token_indices = token_indices_reacts + token_indices_prods
                else:
                    tokens, token_indices = ['O', 'DOUBLEBOND', 'C', 'C', '(', 'DOUBLEBOND', 'O', ')'], [15, 16, 17, 49,
                                                                                                         50, 51, 52, 53]

                df = pd.DataFrame({'Tokens': tokens, 'Token indices': token_indices})
                df.to_excel(writer, sheet_name=f'Total test index {sample}')


for center_category in ['C_substituted', 'Allylic', 'Triple', 'Aryl', 'Carbonyl']:
    find_indices(center_category)