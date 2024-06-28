from autode.reactions.reaction import Reaction
from autode.bond_rearrangement import get_bond_rearrangs

entry_num = "0"
forwards_rxn = "Cn1c(N)nc2ccccc21.CCI>>CC[n+]1c(N)n(C)c2ccccc21.[I-]"
reacts = forwards_rxn.split(">>")[0]
prods = forwards_rxn.split(">>")[1]
backwards_rxn = ">>".join([prods, reacts])


def atom_mapping_forwards(reaction, entry_number):
    br = get_bond_rearrangs(reactant=reaction.reactant,
                            product=reaction.product,
                            name='bond_rearrange',
                            save=True)[0]

    fbond, bbond = br.fbonds[0], br.bbonds[0]

    c_atom = set(fbond).intersection(bbond)
    nu_atom = [element for element in fbond if element not in bbond][0]
    lg_atom = [element for element in bbond if element not in fbond][0]
    with open(f"atom_mapping_out_forwards_{entry_number}.txt", "w") as myfile:
        myfile.write(
            f"fbond = {fbond}\nbbond = {bbond}\nc atom (reacts) = {c_atom}\nnu atom (reacts) = {nu_atom}\nlg atom (reacts) = {lg_atom}")


def atom_mapping_backwards(reaction, entry_number):
    br = get_bond_rearrangs(reactant=reaction.reactant,
                            product=reaction.product,
                            name='bond_rearrange',
                            save=True)[0]

    fbond, bbond = br.fbonds[0], br.bbonds[0]

    c_atom = set(fbond).intersection(bbond)
    lg_atom = [element for element in fbond if element not in bbond][
        0]  # nu and lg are reversed because we've reversed prods and reacts
    nu_atom = [element for element in bbond if element not in fbond][0]
    with open(f"atom_mapping_out_backwards_{entry_number}.txt", "w") as myfile:
        myfile.write(
            f"fbond = {fbond}\nbbond = {bbond}\nc atom in major prod = {c_atom}\nnu atom in major prod  = {nu_atom}\nlg atom in nucleofuge = {lg_atom}")


forwards_rxn = Reaction(smiles=forwards_rxn)
atom_mapping_forwards(forwards_rxn, entry_num)

backwards_rxn = Reaction(smiles=backwards_rxn)
atom_mapping_backwards(backwards_rxn, entry_num)