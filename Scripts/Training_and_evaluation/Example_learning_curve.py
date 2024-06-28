import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

full = pd.read_excel(r"Total_train_results_df.xlsx")


def decay(x, a, b, c):
    y = a * np.exp(-b * x) + c
    return y


rmse = [full["Total test"].to_list()[0]]
sem = [full["Total test"].to_list()[1]]

for percentage in [75, 50, 25, 10, 5, 1]:
    percentage_res = pd.read_excel(
        fr"{percentage}_per_cent_results_df.xlsx")
    percentage_rmse = percentage_res["Total test"].to_list()[0]
    percentage_sem = percentage_res["Total test"].to_list()[1]
    rmse.append(percentage_rmse)
    sem.append(percentage_sem)

fig, ax = plt.subplots(figsize=[5, 5])
x = [100, 75, 50, 25, 10, 5, 1][::-1]
plt.errorbar(x, rmse[::-1], yerr=sem[::-1], fmt="o", ecolor="black", markeredgecolor="black", color="black")

print(rmse[::-1])
print(sem[::-1])
parameters, covariance = curve_fit(decay, x, rmse[::-1])
param_a = parameters[0]
param_b = parameters[1]
param_c = parameters[2]

fit_y = decay(range(1, 201), param_a, param_b, param_c)
g_y = np.gradient(fit_y)

index_zeros = np.where(abs(g_y) < 0.01)
print(index_zeros[0][0] + 1)  # indexing  - we start at 1%, not 0% of the data

plt.plot(range(1, 201), fit_y, '-', color="black")

plt.ylim(0.4, 2)
plt.xlabel("% training data", fontsize=20)
plt.ylabel("RMSE /log$k$", fontsize=20)
plt.xticks(range(0, 225, 50), fontsize=20)
plt.yticks(fontsize=20)