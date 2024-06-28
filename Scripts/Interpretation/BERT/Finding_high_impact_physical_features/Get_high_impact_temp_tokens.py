import pandas as pd

set_of_feats_data = pd.read_excel(
    r"High_impact_features.xlsx",
    sheet_name=None)
token_data = pd.read_excel(
    r"IGs_fold_0.xlsx",
    sheet_name=None)
percentile_data = pd.read_excel(
    "Accurate_and_inaccurate_predictions.xlsx",
    sheet_name=None)

with pd.ExcelWriter(
        r"High_impact_temp_tokens.xlsx") as writer:
    for percentile in ["Accurate predictions", "Inaccurate predictions"]:
        percentile_samples = percentile_data[percentile]["Total test index"].to_list()
        for sample in percentile_samples:
            sample_set_of_feats_data = set_of_feats_data[f"Total test sample {sample}"]
            sample_set_of_feats = sample_set_of_feats_data["Token index"].to_list()
            sample_token_data = token_data[f"Total test sample {sample}"]
            sample_tokens = sample_token_data["Token"].to_list()
            recip_temp_index = sample_tokens.index("[RecipTemp]")
            ion_str_index = sample_tokens.index("[IonStr]")
            temp_token_indices = [i for i in range(recip_temp_index, ion_str_index)]
            sig_temp_token_indices = [temp_token_indices[0]] + temp_token_indices[5:]
            high_impact_temp_token_indices = [item for item in sample_set_of_feats if item in sig_temp_token_indices]
            high_impact_temp_tokens = [sample_tokens[i] for i in high_impact_temp_token_indices]
            sample_df = pd.DataFrame({"Token": high_impact_temp_tokens, "Token index": high_impact_temp_token_indices})
            sample_df.to_excel(writer, sheet_name=f"Total test sample {sample}")
