import pandas as pd

ig_data = pd.read_excel('IGs_fold_0.xlsx',sheet_name=None)
solvent_category_data = pd.read_excel(
    r'Solvent_types.xlsx',
    sheet_name=None)


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


def find_indices(solvent_category):
    with pd.ExcelWriter(
            f'{solvent_category}_solvent_token_indices.xlsx') as writer:
        for percentile in ['Accurate predictions', 'Inaccurate predictions']:
            percentile_solvent_category_data = solvent_category_data[f'{percentile}']

            solvent_category_samples = percentile_solvent_category_data[f'Total test index ({solvent_category.replace("_"," ").lower()})'].to_list()
            solvent_category_samples = [int(item) for item in solvent_category_samples if pd.isna(item) == False]

            for sample in solvent_category_samples:
                sample_ig_data = ig_data[f'Total test index {sample}']
                sample_tokens = sample_ig_data['Token'].to_list()
                sample_tokens = _replace_double_bond(sample_tokens)
                reactant_tokens, reactant_token_indices = _get_rxn_side_tokens(sample_tokens, 'reacts')
                product_tokens, product_token_indices = _get_rxn_side_tokens(sample_tokens, 'prods')

                reactant_molecule_tokens, reactant_molecule_token_indices = _split_into_molecules(reactant_tokens,
                                                                                                  reactant_token_indices)
                solvent_tokens_reacts, solvent_token_indices_reacts = reactant_molecule_tokens[2], \
                                                                      reactant_molecule_token_indices[2]

                product_molecule_tokens, product_molecule_token_indices = _split_into_molecules(product_tokens,
                                                                                                product_token_indices)
                solvent_tokens_prods, solvent_token_indices_prods = product_molecule_tokens[2], \
                                                                    product_molecule_token_indices[2]

                tokens = solvent_tokens_reacts + solvent_tokens_prods
                token_indices = solvent_token_indices_reacts + solvent_token_indices_prods

                df = pd.DataFrame({'Tokens': tokens, 'Token indices': token_indices})
                df.to_excel(writer, sheet_name=f'Total test index {sample}')


for solvent_category in 'Polar_protic', 'Polar_aprotic', 'Non_polar':
    find_indices(solvent_category)