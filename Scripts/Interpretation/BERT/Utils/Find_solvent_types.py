import pandas as pd

percentile_data = pd.read_excel(r'Accurate_and_inaccurate_predictions.xlsx',sheet_name=None)
solvent_data = pd.read_excel(r'Solvent_smiles.xlsx', sheet_name = None)

ext_solvent_data = solvent_data['Total test']

polar_protic = ['CO', 'O']
polar_aprotic = ['CS(C)=O', 'CC#N', 'CC(C)=O', 'CN(C)C=O']
non_polar = ['CCCCCCCCO']

def find_solvents(percentile):
    rxns_with_polar_protic = []
    rxns_with_polar_aprotic = []
    rxns_with_non_polar = []
    percentile_indices = percentile_data[f'{percentile}']['Total test index'].to_list()
    for i,row in ext_solvent_data.iterrows():
        if i in percentile_indices:
            if row['Solvent'] in polar_protic:
                rxns_with_polar_protic.append(i)
            if row['Solvent'] in polar_aprotic:
                rxns_with_polar_aprotic.append(i)
            if row['Solvent'] in non_polar:
                rxns_with_non_polar.append(i)
    df = pd.DataFrame({'Total test index (polar protic)': pd.Series(rxns_with_polar_protic), 'Total test index (polar aprotic)': pd.Series(rxns_with_polar_aprotic), 'Total test index (non polar)': pd.Series(rxns_with_non_polar)})
    return df

top_df = find_solvents('Accurate predictions')
bottom_df = find_solvents('Inaccurate predictions')
with pd.ExcelWriter(r'Solvent_types.xlsx') as writer:
    top_df.to_excel(writer, sheet_name = 'Accurate predictions')
    bottom_df.to_excel(writer, sheet_name = 'Inaccurate predictions')