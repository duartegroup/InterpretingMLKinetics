import pandas as pd
import random
import numpy as np
from scipy.stats import sem
from scipy.stats import spearmanr

importance_data = pd.read_excel(
    rf"Mean_IGs_for_each_test_sample.xlsx",
    sheet_name=None)

solvent_categories = ["Polar_protic", "Polar_aprotic", "Non_polar"]


def get_cov():
    rs = []
    ps = []

    for iteration in range(50):
        token_A = []
        token_B = []
        for category in solvent_categories:
            token_index_data = pd.read_excel(
                rf"{category}_solvent_token_indices.xlsx",
                sheet_name=None)
            for sheet_name in token_index_data.keys():
                sheet = token_index_data[sheet_name]
                token_indices = sheet["Token indices"].to_list()  # summed over reacts and prods, no split
                token_indices = random.sample(token_indices, 2)
                sample = sheet_name.split(" ")[-1]
                sample_importance_data = importance_data[f"Total test sample {sample}"]
                sems = [row["Sem"] for i, row in sample_importance_data.iterrows() if
                        row["Token index"] in token_indices]
                token_A.append(sems[0])
                token_B.append(sems[1])

        r = spearmanr(token_A, token_B)[0]
        p = spearmanr(token_A, token_B)[1]
        rs.append(r)
        ps.append(p)

    print(round(np.mean(rs), 2), "+-", round(sem(rs), 2), len(token_A))
    print(round(np.mean(ps), 2), "+-", round(sem(ps), 2))


get_cov()