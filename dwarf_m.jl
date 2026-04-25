using CSV
using DataFrames
using Plots
using Main.SSA_module
gr()

mat = CSV.read("Dwarfst.dat", DataFrame; header = false)
rename!(mat, :Column1 => :timeseries)

time_series = mat.timeseries
plot(time_series)

X = SSA_module.hankel(time_series, 200)
L = 200
K = length(time_series) - L + 1

U, Σ, V = SSA_module.svd_a(X)

plot((Σ))
elementary_recons = SSA_module.elementary_reconstructions(U, Σ, V, 30)

w_heatmap = SSA_module.w_heatmap(elementary_recons, L, K, 30)

groups = [1:2, 3:5, 1:11]
grouped_series = SSA_module.grouped_series(elementary_recons, L, K, groups)

plot(grouped_series[1])
plot(grouped_series[end])
plot!( title = "Reconstructed Series",
    xlabel = "Time",
    ylabel = "Signal",
    titlefont = font(8),
    guidefont = font(8),)
plot!(time_series)

relative_error(time_series, grouped_series[end])