import post_process_rewrite
import matplotlib.pyplot as plt

dbname = "C:/Users/rthill/Documents/MS-Thesis/Optimizer_Output"
testobj = post_process_rewrite.NSGA_PostProcess(dbname)
optvar_front = testobj.getFront(100,50,[1,2])
optvars = testobj.getIndVariables(optvar_front)

#generate just Pareto front for single generation
plot_pareto = testobj.SingleFront(10, 0, [1,2], 1.2, "SingleGen_Pareto_Plot")
plt.show()

#generate and show all fronts for single generation
plot_data = testobj.GenPlot2D(1, [1,2], [1.75,1.75], "SingleGen_All_Fronts_Plot", True)
plt.show()

#generate a gif to show front progression over generations
movie = testobj.ObjMovie(1,10,[1,2],[1.3,1.3],"All_Fronts_Movie")