module SSA_module
export hankel, optim_svd_a, svd_a, elementary_reconstructions,  w_heatmap, grouped_series, relative_error
using Plots 
backend(:plotly)
using LinearAlgebra
using LowRankApprox

    function hankel(data, embedding_window)
        N = length(data)
        L = embedding_window
        K = N - L + 1

        X = zeros(K, L)

        for (i,j) in Iterators.product(1:K,1:L)
            X[i,j] = data[(i + j - 1)]
        end
    # X = 1/floor(sqrt(length(data))) * X
        return X
    end

    function optim_svd_a(A, r::Int)
        U, S, V = psvd(A, rank=r)
        return U, S, V
    end

    function svd_a(A)
        U, Σ, V = svd(A)

        return U, Σ, V
    end 

    function weights(L, K)
        N = K + L - 1
        weight = zeros(N)

        weight = [min(n, L, N - n + 1) for n in 1:N]
    end 

    function hankelized_series(u, σ, v)
        K = size(u, 1)
        L = size(v, 1)
        N = K + L - 1
        a_diag_sum = zeros(N)
        counts = weights(L, K)

        for j in 1:L
            v_j = v[j]
            for i in 1:K
                n = i + j - 1
                a_diag_sum[n] += σ * u[i] * v_j
            end
        end

        return a_diag_sum ./ counts
    end 
# It might be better to just do the hankelized reconstruction in the elementary reconstruction as it isnt called anymore after that point. 
    function elementary_reconstructions(U, Σ, V, k::Int)
        K = size(U, 1)
        L = size(V, 1)

        N = K + L - 1
        
        elementary_recons = [zeros(N) for i in 1:k]

        for i in 1:k
            elementary_recons[i] = hankelized_series(@view(U[:,i]),Σ[i],@view(V[:,i]))
        end
        
        return elementary_recons
    end 

    function w_heatmap(elementary_recons::AbstractVector{<:AbstractVector}, L, K,  k::Int)
        ρ = zeros(k,k)
        weight = weights(L, K)
        N = L + K - 1
        recon_matrix = zeros(k,N)
        weighted_matrix = zeros(k,N)

        for i in 1:k
            recon_matrix[i,:] = elementary_recons[i]
            weighted_matrix[i,:] = elementary_recons[i] .* weight
        end 

        W = recon_matrix * weighted_matrix'
        norms = sqrt.(diag(W))

        ρ = W ./ (norms * norms')

        ρ = abs.(ρ)

        display(heatmap(
            ρ,
            c = cgrad(:grays, rev=true),  # 🔥 reverse it
            clim = (0, 1),
            xlabel="i",
            ylabel="j",
            title="Grayscale Heatmap"
        ))

        return ρ
    end 

    function grouped_series(elementary_recons::AbstractVector{<:AbstractVector}, L, K, groups)
        n = length(groups)
        N = K + L - 1
        grouped_recons = [zeros(N) for i in 1:n]

        for i in 1:n 
            for j in groups[i]
                grouped_recons[i] += elementary_recons[j]
            end 
        end 

        return grouped_recons
    end 
    function relative_error(x::AbstractVector, y::AbstractVector)
    s = zero(eltype(x))
    @inbounds @simd for i in eachindex(x, y)
        d = x[i] - y[i]
        s += d*d
    end

    rmse = sqrt(s / length(x))
    range = maximum(x) - minimum(x)

    return rmse / range
end
end