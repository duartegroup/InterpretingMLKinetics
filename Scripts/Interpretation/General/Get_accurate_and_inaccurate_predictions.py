# Example for RF
import pandas as pd

predictions_data = pd.read_excel(r'Average_predictions_over_folds_sorted_by_error.xlsx')
twenty_five_percent = int(round(len(predictions_data) * 0.25,0))
top_25 = predictions_data.iloc[:twenty_five_percent,:]
bottom_25 = predictions_data.iloc[-twenty_five_percent:,:]

with pd.ExcelWriter('Accurate_and_inaccurate_predictions.xlsx') as writer:
    top_25.to_excel(writer, sheet_name = 'Accurate predictions')
    bottom_25.to_excel(writer, sheet_name = 'Inaccurate predictions')