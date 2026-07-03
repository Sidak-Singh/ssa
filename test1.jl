using Random
using Plots
using Main.SSA_module
f(x) = 0.3*x + (1+0.1*x)*sin((2π/20)*x)

function gaussian_noise(σ,signal)
    signal_dim=length(signal)
    noise = σ .* randn(signal_dim)
    return noise
end

i = 1:400
x_true=f.(i)
x = f.(i) + gaussian_noise(5,i)

plot(x)

X = hankel(x, 180)
L = 180
K = length(x) - L + 1
U, Σ, V = svd_a(X)





plot(log.((Σ)))
elementary_recons = elementary_reconstructions(U, Σ, V, 30)

ρ = w_heatmap(elementary_recons, L, K, 30)

using Clustering, StatsPlots

D = 1 .- ρ
D = (D + D') / 2   # optional safety

hc = hclust(D, linkage=:average)

plot(hc, labels=1:size(ρ,1))
clusters = cutree(hc, h=0.2)
function clusters_to_groups(clusters::AbstractVector{<:Integer})
    k = maximum(clusters)
    groups = [Int[] for _ in 1:k]

    for (i, c) in enumerate(clusters)
        push!(groups[c], i)
    end

    return filter(!isempty, groups)
end
gr=clusters_to_groups(clusters)
groups = []

reconstructed_series = grouped_series(elementary_recons,L, K, groups)

plot([x for x in reconstructed_series])
plot(reconstructed_series[2]+reconstructed_series[4])
plot!([(1+0.1*x)*sin((2π/20)*x) for x in 1:400])
plot(reconstructed_series[1] + reconstructed_series[2] + reconstructed_series[3] + reconstructed_series[4])
plot!(x_true)

plot(V[:,1])


groups_components = [[1,4],vcat(2:3, 5:6)]
reconstructed_series = grouped_series(elementary_recons, L, K, groups_components)
trend = reconstructed_series[1]
oscillation = reconstructed_series[2]

true_trend = [0.3*x for x in 1:400]
true_oscillation = [(1+0.1*x)*sin((2π/20)*x) for x in 1:400]

relative_error(true_trend, trend)
relative_error(true_oscillation, oscillation)