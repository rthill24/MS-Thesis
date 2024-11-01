import post_process
import matplotlib.pyplot as plt

dbname = "C:/Users/rthill/Documents/MS-Thesis/Optimizer_Output"
testobj = post_process.NSGA_PostProcess(dbname)
optvar_front = testobj.getFront(100,0,[1,2])
optvars = testobj.getIndVariables(optvar_front)

#generate just Pareto front for single generation
plot_pareto = testobj.SingleFront(dbname, 10,1, [1,2], 1)
plt.show()

#generate and show all fronts for single generation
plot_data = testobj.GenPlot2D(10, [1,2], [1,1], "SingleGen_All_Fronts_Plot")
plt.show()

movie = testobj.ObjMovie(1,10,[1,2],[1,1],"All_Fronts_Movie")