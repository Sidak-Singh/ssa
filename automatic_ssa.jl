module automatic_ssa
include("ssa_module.jl")
include("./kneedle_jl/kneedle_algo.jl")
export auto_ssa

function auto_ssa(x, y, L)
N = length(x)
K = N - L + 1

X = SSA_module.hankel(x, L)
U, Σ, V = SSA_module.svd_a(X)

idx, xk, yk = kneedle_algo.kneedle(1:length(Σ), Σ)

elementary_recons = SSA_module.elementary_reconstructions(U, Σ, V, length(Σ))

groups = [1:floor(idx), floor(idx)+1:length(Σ)]
reconstructed_series = SSA_module.grouped_series(elementary_recons, L, K, groups)

return smooth_series, residual = reconstructed_series[1], reconstructed_series[2]
end

end