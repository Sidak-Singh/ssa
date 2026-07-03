using Plots
using Main.SSA
sample = 100

t = 0:(1/sample):(10-1/sample)

x_signal = sin.(2π.*t) .+ 2sin.(6π.*t) .+ sin.(2π.*(t ./ 2))

function gaussian_noise(σ,signal)
    signal_dim=length(signal)
    noise = σ .* randn(signal_dim)
    return noise
end

x = x_signal .+ gaussian_noise(2, x_signal)

plot(x_signal)
plot!(x)

X = hankel(x,400)

U, Σ, V = svd_and_analysis(X)
s=plot(Σ, xlims=(1,20))


W = abs.(w_matrix(X,20))
heatmap(
    W,
   c = cgrad(:grays, rev=true),  # 🔥 reverse it
    clim = (0, 1),
    xlabel="i",
    ylabel="j",
    title="Grayscale Heatmap"
)

grouped_matrices = grouped_svd(U, Σ, V, [1:2,3:4,5:6, 1:6, 9:12]);

for (i, z) in enumerate(grouped_matrices)
    (data, _) = hankelization(z)
    p = plot(data)
    display(p) 
end

smoothing, _ = hankelization(grouped_matrices[4])
m = plot(smoothing, label="Smoothed Signal")
plot!(x_signal, label="Original Signal")
# plot!(x, label="Noisy Signal")

plot!(size=(1400, 800)) 

gui(m)