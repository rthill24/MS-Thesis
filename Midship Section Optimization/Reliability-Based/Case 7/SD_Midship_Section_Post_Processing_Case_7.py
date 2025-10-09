import sys
sys.path.insert(0, 'C:/Users/rthill/Documents/MS-Thesis')
import post_process
import matplotlib.pyplot as plt


dbname = "C:/Users/rthill/Documents/MS-Thesis/Midship Section Optimization/Reliability-Based/Case 7/SD_Midship_Section_Optimizer_Output_Case_7"
testobj = post_process.NSGA_PostProcess(dbname)
#optvar_front = testobj.getFront(200,0,[1,2])
#optvars = testobj.getIndVariables(optvar_front)
#print (optvars)

#generate just Pareto front for single generation
plot_pareto = testobj.SingleFront(200, 0, [1,2], 1.2, "SD_Midship_Section_SingleGen_Pareto_Plot")
#plt.show()

#generate and show all fronts for single generation
plot_data = testobj.GenPlot2D(200, [1,2], [1.75,1.75], "SD_Midship_Section_SingleGen_All_Fronts_Plot", 0, False)
#plt.show()

#generate a gif to show front progression over generations
movie = testobj.ObjMovie(190,200,[1,2], 1.0,"SD_Midship_Section_All_Fronts_Movie", "C:/Users/rthill/Documents/MS-Thesis/Midship Section Optimization/Reliability-Based/Case 7")

#feas_stats = testobj.Feasibility_Stats(list(range(120, 123)))
#print(feas_stats)