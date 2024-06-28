import pandas as pd

rxn_center_index_data = pd.read_excel(
    r'Rxn_center_indices_resonance_correction.xlsx')
percentile_data = pd.read_excel(
    'Accurate_and_inaccurate_predictions.xlsx',
    sheet_name=None)
ig_data = pd.read_excel(
    r'IGs_fold_0.xlsx',
    sheet_name=None)


# only need fold 0 to map token to index - same in every fold


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


def find_token_index(percentile):
    percentile_indices = percentile_data[f'{percentile}']['Total test index'].to_list()

    percentile_nu_atom_indices = []
    percentile_nu_atom_tokens = []
    percentile_nu_atom_token_indices = []

    percentile_c_atom_indices = []
    percentile_c_atom_tokens = []
    percentile_c_atom_token_indices = []

    percentile_lg_atom_indices = []
    percentile_lg_atom_tokens = []
    percentile_lg_atom_token_indices = []

    percentile_lg_atom_indices_nucleofuge = []
    percentile_lg_atom_tokens_nucleofuge = []
    percentile_lg_atom_token_indices_nucleofuge = []

    percentile_c_atom_indices_mp = []
    percentile_c_atom_tokens_mp = []
    percentile_c_atom_token_indices_mp = []

    percentile_nu_atom_indices_mp = []
    percentile_nu_atom_tokens_mp = []
    percentile_nu_atom_token_indices_mp = []

    for sample in percentile_indices:
        sample_ig_data = ig_data[f'Total test index {sample}']
        sample_tokens = sample_ig_data['Token'].to_list()
        sample_tokens = _replace_double_bond(sample_tokens)

        reactant_tokens, reactant_token_indices = _get_rxn_side_tokens(sample_tokens, 'reacts')
        product_tokens, product_token_indices = _get_rxn_side_tokens(sample_tokens, 'prods')

        reactant_molecule_tokens, reactant_molecule_token_indices = _split_into_molecules(reactant_tokens,
                                                                                          reactant_token_indices)
        product_molecule_tokens, product_molecule_token_indices = _split_into_molecules(product_tokens,
                                                                                        product_token_indices)

        nu_tokens, nu_token_indices = reactant_molecule_tokens[0], reactant_molecule_token_indices[0]
        substrate_tokens, substrate_token_indices = reactant_molecule_tokens[1], reactant_molecule_token_indices[1]
        major_prod_tokens, major_prod_token_indices = product_molecule_tokens[0], product_molecule_token_indices[0]
        nucleofuge_tokens, nucleofuge_token_indices = product_molecule_tokens[1], product_molecule_token_indices[1]

        nu_tokens_atoms_only, nu_token_indices_atoms_only = _get_atom_tokens_only(nu_tokens, nu_token_indices)
        substrate_tokens_atoms_only, substrate_token_indices_atoms_only = _get_atom_tokens_only(substrate_tokens,
                                                                                                substrate_token_indices)
        major_prod_tokens_atoms_only, major_prod_token_indices_atoms_only = _get_atom_tokens_only(major_prod_tokens,
                                                                                                  major_prod_token_indices)
        nucleofuge_tokens_atoms_only, nucleofuge_token_indices_atoms_only = _get_atom_tokens_only(nucleofuge_tokens,
                                                                                                  nucleofuge_token_indices)

        sample_rxn_center_index_row = rxn_center_index_data.iloc[sample, :]

        nu_atom_index = sample_rxn_center_index_row['Nu index']
        percentile_nu_atom_indices.append(nu_atom_index)
        percentile_nu_atom_tokens.append(nu_tokens_atoms_only[nu_atom_index])
        percentile_nu_atom_token_indices.append(nu_token_indices_atoms_only[nu_atom_index])

        c_atom_index = sample_rxn_center_index_row['C index']
        percentile_c_atom_indices.append(c_atom_index)
        percentile_c_atom_tokens.append(substrate_tokens_atoms_only[c_atom_index])
        percentile_c_atom_token_indices.append(substrate_token_indices_atoms_only[c_atom_index])

        lg_atom_index = sample_rxn_center_index_row['LG index']
        percentile_lg_atom_indices.append(lg_atom_index)
        percentile_lg_atom_tokens.append(substrate_tokens_atoms_only[lg_atom_index])
        percentile_lg_atom_token_indices.append(substrate_token_indices_atoms_only[lg_atom_index])

        lg_atom_index_nucleofuge = sample_rxn_center_index_row['LG index (nucleofuge)']
        percentile_lg_atom_indices_nucleofuge.append(lg_atom_index_nucleofuge)
        percentile_lg_atom_tokens_nucleofuge.append(nucleofuge_tokens_atoms_only[lg_atom_index_nucleofuge])
        percentile_lg_atom_token_indices_nucleofuge.append(
            nucleofuge_token_indices_atoms_only[lg_atom_index_nucleofuge])

        c_atom_index_mp = sample_rxn_center_index_row['C index (major prod)']
        percentile_c_atom_indices_mp.append(c_atom_index_mp)
        percentile_c_atom_tokens_mp.append(major_prod_tokens_atoms_only[c_atom_index_mp])
        percentile_c_atom_token_indices_mp.append(major_prod_token_indices_atoms_only[c_atom_index_mp])

        nu_atom_index_mp = sample_rxn_center_index_row['Nu index (major prod)']
        percentile_nu_atom_indices_mp.append(nu_atom_index_mp)
        percentile_nu_atom_tokens_mp.append(major_prod_tokens_atoms_only[nu_atom_index_mp])
        percentile_nu_atom_token_indices_mp.append(major_prod_token_indices_atoms_only[nu_atom_index_mp])

    df = pd.DataFrame({'Total test index': percentile_indices, 'Nu atomic index': percentile_nu_atom_indices,
                       'Nu atom token': percentile_nu_atom_tokens,
                       'Nu atom token indices': percentile_nu_atom_token_indices,
                       'C atomic index': percentile_c_atom_indices, 'C atom token': percentile_c_atom_tokens,
                       'C atom token indices': percentile_c_atom_token_indices,
                       'LG atomic index': percentile_lg_atom_indices, 'LG atom token': percentile_lg_atom_tokens,
                       'LG atom token indices': percentile_lg_atom_token_indices,
                       'LG atomic index (nucleofuge)': percentile_lg_atom_indices_nucleofuge,
                       'LG atom token (nucleofuge)': percentile_lg_atom_tokens_nucleofuge,
                       'LG atom token indices (nucleofuge)': percentile_lg_atom_token_indices_nucleofuge,
                       'C atomic index (major prod)': percentile_c_atom_indices_mp,
                       'C atom token (major prod)': percentile_c_atom_tokens_mp,
                       'C atom token indices (major prod)': percentile_c_atom_token_indices_mp,
                       'Nu atomic index (major prod)': percentile_nu_atom_indices_mp,
                       'Nu atom token (major prod)': percentile_nu_atom_tokens_mp,
                       'Nu atom token indices (major prod)': percentile_nu_atom_token_indices_mp})
    return df


top_df = find_token_index('Accurate predictions')
bottom_df = find_token_index('Inaccurate predictions')

with pd.ExcelWriter('Rxn_center_token_indices.xlsx') as writer:
    top_df.to_excel(writer, 'Accurate predictions')
    bottom_df.to_excel(writer, 'Inaccurate predictions')