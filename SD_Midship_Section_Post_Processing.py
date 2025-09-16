import post_process
import matplotlib.pyplot as plt

dbname = "C:/Users/rthill/Documents/MS-Thesis/SD_Midship_Section_Optimizer_Output"
testobj = post_process.NSGA_PostProcess(dbname)
optvar_front = testobj.getFront(100,0,[1,2])
optvars = testobj.getIndVariables(optvar_front)
print (optvars)

#generate just Pareto front for single generation
plot_pareto = testobj.SingleFront(100, 50, [1,2], 1.2, "SD_Midship_Section_SingleGen_Pareto_Plot")
plt.show()

#generate and show all fronts for single generation
plot_data = testobj.GenPlot2D(10, [1,2], [1.75,1.75], "SD_Midship_Section_SingleGen_All_Fronts_Plot", True)
plt.show()

#generate a gif to show front progression over generations
movie = testobj.ObjMovie(1,10,[1,2],[1.3,1.3],"SD_Midship_Section_All_Fronts_Movie", option="SD")

