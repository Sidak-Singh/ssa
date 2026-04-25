using CSV
using DataFrames
using Plots
using Main.SSA_module
using BenchmarkTools
using Peaks

data= CSV.read("./sunspots/sunspot.csv", DataFrame, header=false, delim=';')
rename!(data, :Column5 => :series)

ts = data.series
L = 32000
K = length(ts) - L + 1

X = SSA_module.hankel(ts, L)
U, Σ, V = SSA_module.optim_svd_a(X, 20)
elementary_recons = SSA_module.elementary_reconstructions(U, Σ, V, 20)
plot(log.(Σ))
W_matrix = SSA_module.w_heatmap(elementary_recons, L, K, 20)

groups = [1, 2:3]
grouped_recons = SSA_module.grouped_series(elementary_recons, L, K, groups)

plot(grouped_recons[1])
plot(grouped_recons[2])
plot(ts)
title!("Smoothed Trend and Oscillation component of Daily Sunspots")