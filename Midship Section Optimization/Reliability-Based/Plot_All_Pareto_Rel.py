import matplotlib.pyplot as plt
import math
import numpy as np
import sys
sys.path.insert(0, 'C:/Users/rthill/Documents/MS-Thesis')
import post_process

dbname_list = []
for i in range(1,5):
    name = "C:/Users/rthill/Documents/MS-Thesis/Midship Section Optimization/Reliability-Based/Case " + str(i) + "/SD_Midship_Section_Optimizer_Output_Case_" + str(i)
    dbname_list.append(name)

comb_list = []
styles = ['o-', 's--', '^-', 'd:', 'x-', '*--', 'p-', 'h--']  # Unique styles for up to 8 cases

for idx, i in enumerate(dbname_list):
    testobj = post_process.NSGA_PostProcess(i)
    data = testobj.frontPlotFormat(200, 0, [1, 2])
    for i in range(len(data)):
            if math.isnan(data[i][0]) or math.isnan(data[i][1]):
                data = np.delete(data, i, 0)
    comb_list.append(data)
    x = [item[0] for item in data]
    y = [item[1] for item in data]
    plt.scatter(x, y, marker=styles[idx % len(styles)][0], label=f'Case {idx+1}')

plt.xlabel('Cost [$]')
plt.ylabel('Weight [kg]')
plt.legend()
plt.title('Pareto Fronts for All Reliability-Based Cases')
plt.savefig('All_Reliability_Based_Cases_Pareto_Fronts.png', dpi=144)
plt.show()




