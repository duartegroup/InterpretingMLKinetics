import pandas as pd

ig_data = pd.read_excel(
    r"IGs_fold_0.xlsx",
    sheet_name=None)
set_of_feats_data = pd.read_excel(
    r"High_impact_features.xlsx",
    sheet_name=None)
token_data = ig_data
percentile_data = pd.read_excel(
    "Accurate_and_inaccurate_predictions.xlsx",
    sheet_name=None)
solvent_data = pd.read_excel(
    r"Solvent_smiles.xlsx",
    sheet_name=None)

solvents = solvent_data["Total test"]["Solvent"].to_list()


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
        if (token.isalpha() == True and token != "DOUBLEBOND") or ("[" in token):
            molecule_tokens_atoms_only.append(token)
            molecule_token_indices_atoms_only.append(index)

    return molecule_tokens_atoms_only, molecule_token_indices_atoms_only


def _get_rxn_side_tokens(sample_tokens, side):
    rxn_arrow_token_index = sample_tokens.index(">>")
    recip_temp_index = sample_tokens.index("[RecipTemp]")
    if side == "reacts":
        rxn_side_token_indices = [x for x in range(len(sample_tokens)) if x < rxn_arrow_token_index and x != 0]
    else:
        rxn_side_token_indices = [x for x in range(len(sample_tokens)) if
                                  x > rxn_arrow_token_index and x < recip_temp_index]

    rxn_side_tokens = [sample_tokens[x] for x in rxn_side_token_indices]

    return rxn_side_tokens, rxn_side_token_indices


def _replace_double_bond(tokens):
    for x in range(len(tokens)):
        if tokens[x] == 0:
            tokens[x] = "DOUBLEBOND"

    return tokens


with pd.ExcelWriter(
        f"High_impact_solvent_tokens.xlsx") as writer:
    for percentile in ["Accurate predictions", "Inaccurate predictions"]:
        percentile_samples = percentile_data[percentile]["Total test index"].to_list()
        for sample in percentile_samples:
            if solvents[sample] != "solvent is nu":
                sample_ig_data = ig_data[f"Total test sample {sample}"]
                sample_tokens = sample_ig_data["Token"].to_list()
                sample_tokens = _replace_double_bond(sample_tokens)
                reactant_tokens, reactant_token_indices = _get_rxn_side_tokens(sample_tokens, "reacts")

                reactant_molecule_tokens, reactant_molecule_token_indices = _split_into_molecules(reactant_tokens,
                                                                                                  reactant_token_indices)
                solvent_tokens_reacts, solvent_token_indices_reacts = reactant_molecule_tokens[2], \
                                                                      reactant_molecule_token_indices[2]

                tokens = solvent_tokens_reacts
                token_indices = solvent_token_indices_reacts

                sample_set_of_feats_data = set_of_feats_data[f"Total test sample {sample}"]
                sample_set_of_feats = sample_set_of_feats_data["Token index"].to_list()

                high_impact_tokens = []
                high_impact_token_indices = []
                for i in range(len(token_indices)):
                    if token_indices[i] in sample_set_of_feats:
                        high_impact_token_indices.append(token_indices[i])
                        high_impact_tokens.append(tokens[i])

                df = pd.DataFrame({"Token": high_impact_tokens, "Token indices": high_impact_token_indices})
                df.to_excel(writer, sheet_name=f"Total test sample {sample}")

