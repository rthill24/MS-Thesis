import post_process
import matplotlib.pyplot as plt
from operator import itemgetter

dbname = "C:/Users/Administrator/Documents/Michigan Documents/First Term/Master's Thesis/Code/Working/test_sql"
testobj = post_process.NSGA_PostProcess(dbname)
retlist = testobj.getFront(100,0,[1,2])
inputs = testobj.getIndVariables(retlist)

plot_data = testobj.frontPlotFormat(100,1,[1,2])

xvector = list(map(itemgetter(0), plot_data))
yvector = list(map(itemgetter(1), plot_data))

plt.plot(xvector, yvector, 'ro')
plt.xlabel("Production Cost")
plt.ylabel("Volume")
plt.title("Pareto Front")
plt.axis([min(xvector)/1.15, max(xvector)*1.15, min(yvector)/1.15, max(yvector)*1.15])
plt.show()

""" plot_data = testobj.GenPlot2D(100, [1,2], [1,1], "Pareto_Test_Plot")
plt.show(plot_data) """