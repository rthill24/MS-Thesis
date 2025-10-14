import matplotlib.pyplot as plt
import math
import numpy as np
import sys
sys.path.insert(0, 'C:/Users/rthill/Documents/MS-Thesis')
import post_process

dbname_list = []

#cycle through static FOS cases
for i in range(1,9):
    name = "C:/Users/rthill/Documents/MS-Thesis/Midship Section Optimization/First Yield Moment/Static FOS/Case " + str(i) + "/SD_Midship_Section_Optimizer_Output_Case_" + str(i)
    dbname_list.append(name)

#cycle through reliability based cases
for i in range(1,9):
    name = "C:/Users/rthill/Documents/MS-Thesis/Midship Section Optimization/First Yield Moment/Reliability-Based/Case " + str(i) + "/SD_Midship_Section_Optimizer_Output_Case_" + str(i)
    dbname_list.append(name)

#cycle through return to optimality & reliability-based cases
for i in range(1,9):
    name = "C:/Users/rthill/Documents/MS-Thesis/Midship Section Optimization/First Yield Moment/RB & Return to Optimality/Case " + str(i) + "/SD_Midship_Section_Optimizer_Output_Case_" + str(i)
    dbname_list.append(name)

comb_list = []
# Define three plot styles for the three groups of cases
group_styles = ['o', 's', '^']
group_labels = [
    'Static FOS',
    'Reliability-Based',
    'RB & Return to Optimality'
]

cases_per_group = 8

# Define three plot styles and colors for the three groups of cases
group_styles = ['o', 's', '^']
group_colors = ['tab:green', 'tab:red', 'tab:blue']

for idx, i in enumerate(dbname_list):
    testobj = post_process.NSGA_PostProcess(i)
    data = testobj.frontPlotFormat(200, 0, [1, 2])
    # Remove rows with NaN values
    data = np.array([row for row in data if not (math.isnan(row[0]) or math.isnan(row[1]))])
    comb_list.append(data)
    x = [item[0] for item in data]
    y = [item[1] for item in data]
    group_idx = idx // cases_per_group  # Use integer division
    marker = group_styles[group_idx % len(group_styles)]
    color = group_colors[group_idx % len(group_colors)]
    # Only add label once per group
    label = group_labels[group_idx] if idx % cases_per_group == 0 else None
    plt.scatter(x, y, marker=marker, color=color, label=label)

plt.xlabel('Cost [$]')
plt.ylabel('Weight [kg]')
plt.legend()
plt.title('Pareto Fronts for All FYM Cases')
plt.savefig('All_FYM_Cases_Pareto_Fronts.png', dpi=144)
plt.show()




