using Plots
backend(:plotly)
using Main.kneedle_algo
x = 0:0.001:1
y = 1 ./ (x .+ 0.001)
plot(x,y)
idx, xk,yk = kneedle_algo.kneedle(x,y,curve=:elbow)