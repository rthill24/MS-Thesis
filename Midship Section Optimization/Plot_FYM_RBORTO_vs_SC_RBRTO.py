import matplotlib.pyplot as plt
import math
import numpy as np
import sys
sys.path.insert(0, 'C:/Users/rthill/Documents/MS-Thesis')
import post_process
from matplotlib.lines import Line2D

dbname_list = []

#cycle through FYM return to optimality & reliability-based cases
for i in range(1,9):
    name = "C:/Users/rthill/Documents/MS-Thesis/Midship Section Optimization/First Yield Moment/RB & Return to Optimality/Case " + str(i) + "/SD_Midship_Section_Optimizer_Output_Case_" + str(i)
    dbname_list.append(name)

#cycle through SC return to optimality & reliability-based cases
for i in range(1,9):
    name = "C:/Users/rthill/Documents/MS-Thesis/Midship Section Optimization/Smith Collapse/RB & Return to Optimality/Case " + str(i) + "/SD_Midship_Section_Optimizer_Output_Case_" + str(i)
    dbname_list.append(name)

comb_list = []
# Define two plot styles for the two groups of cases
group_styles = ['o', 's']
group_labels = [
    'R-B',
    'P-B'
]

cases_per_group = 8

# Customize this to change marker size (in points)
# This size is used for Line2D legend markers as "markersize" (points)
# and for scatter points as area = marker_size**2 (points^2)
marker_size = 5

# Define two plot colors for the two groups and distinct markers per case
group_colors = ['tab:green', 'tab:red']
case_markers = ['o', 's', '^', 'v', '<', '>', 'D', 'P']

for idx, i in enumerate(dbname_list):
    testobj = post_process.NSGA_PostProcess(i)
    data = testobj.frontPlotFormat(200, 0, [1, 2])
    # Remove rows with NaN values
    data = np.array([row for row in data if not (math.isnan(row[0]) or math.isnan(row[1]))])
    comb_list.append(data)
    x = [item[0] for item in data]
    y = [item[1] for item in data]
    group_idx = idx // cases_per_group  # Use integer division
    case_idx = idx % cases_per_group
    marker = case_markers[case_idx]
    color = group_colors[group_idx]
    # add a thin black outline to each marker; scatter 's' is area in points^2
    plt.scatter(x, y, marker=marker, facecolors=color, edgecolors='k', linewidths=0.05, s=marker_size**2)

plt.xlabel('Cost [$]')
plt.ylabel('Weight [kg]')
plt.title('Static FOS vs RB & RTO Pareto Fronts')

# Build legend handles for groups (colors) with black outline
group_handles = [
    Line2D([0], [0], marker='o', color='w',
           markerfacecolor=group_colors[g], markeredgecolor='k',
           markeredgewidth=0.6, markersize=marker_size, label=group_labels[g])
    for g in range(len(group_labels))
]

# Build legend handles for cases (markers) with black outline and white fill
case_handles = [
    Line2D([0], [0], marker=case_markers[c], color='w', linestyle='None',
           markerfacecolor='white', markeredgecolor='k', markeredgewidth=0.6,
           markersize=marker_size, label=f'Case {c+1}')
    for c in range(cases_per_group)
]

ax = plt.gca()
legend1 = ax.legend(handles=group_handles, loc='upper right', title='Group Type')
ax.add_artist(legend1)
# place case legend to the right of the plot
case_legend = ax.legend(handles=case_handles, loc='center left', bbox_to_anchor=(1.02, 0.5), title='Case Marker')

# make room on the right for the external legend
plt.subplots_adjust(right=0.78)

plt.savefig('SC_RBRTO_vs_FYM_RBORTO_Paper.png', dpi=144, bbox_inches='tight')
plt.show()




