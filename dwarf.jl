using CSV
using DataFrames
using Plots
using Main.SSA


mat = CSV.read("Dwarfst.dat", DataFrame; header = false)
rename!(mat, :Column1 => :timeseries)

time_series = mat.timeseries
plot(time_series)

X = hankel(time_series, 200)

(U, Σ, V) = svd_and_analysis(X)

# Analysis
n=length(Σ)
singular_value_plot = plot(1:n, Σ, xlims=(1,30))
title!("Singular Values of the Trajectory Matrix")

# eigenvalues analysis
plot(title="Right Singular Eigenvectos")
plot!(1:n, V[:,4])

plot(V[:,1],V[:,2] )
plot(V[:,2],V[:,4] )

grouped_matrices = grouped_svd(U, Σ, V, [1:3, 3:6, 7:8, 9:10])
# (trend,_)  = hankelization(grouped_matrices[1])
# (reconstructed,_) = hankelization(grouped_matrices[4])
# plot(trend)
# plot(reconstructed)
# plot!(time_series)

for (i, x) in enumerate(grouped_matrices)
    (data, _) = hankelization(x)
    p = plot(data)
    plot!(time_series)
    display(p)
    
end

plot(hankelization(grouped_matrices[3])[1])
W = abs.(w_matrix(X, 30))

df = DataFrame(W, :auto) # ':auto' creates column names x1, x2, ...
vscodedisplay(df)
a = hankelization(grouped_matrices[1])
typeof(a)

heatmap(
    W,
   c = cgrad(:grays, rev=true),  # 🔥 reverse it
    clim = (0, 1),
    xlabel="i",
    ylabel="j",
    title="Grayscale Heatmap"
)