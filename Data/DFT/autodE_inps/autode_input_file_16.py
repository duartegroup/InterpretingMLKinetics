import autode
from autode.reactions.reaction import Reaction
from autode.wrappers.keywords import SinglePointKeywords


autode.Config.ORCA.keywords.sp = SinglePointKeywords(['SP', 'CCSD(T)', 'def2-TZVP'])
autode.utils.log_time()

def reaction_profile(reaction):
    reaction.calculate_reaction_profile(free_energy=True, with_complexes=True)


entry_number = "343"
smiles = "CCCBr.[Cl-]>>CCCCl.[Br-]"
solvent = "dmf"
temp = 298.15

sn2 = Reaction(name=f"dft_benchmark_data_point_{entry_number}_outputs", smiles=smiles, solvent_name=solvent, temp=temp)
reaction_profile(sn2)
