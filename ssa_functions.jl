module SSA
    export hankel, svd_and_analysis, grouped_svd, hankelization, w_matrix
    using Plots
    backend(:plotly)
    using LinearAlgebra




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

    function svd_and_analysis(A)
        (U, Σ, V) = svd(A)
        return (U, Σ, V)
    end


    function grouped_svd(U, Σ, V, groups)
    K = size(U, 1)
    L = size(V, 1)

    n=length(groups)
    grouped_matrices = [zeros(K,L) for i in 1:n]

    for i in 1:n 
        for j in groups[i]
            grouped_matrices[i] += Σ[j] * U[:,j] * V[:,j]'
        end 
    end
    return grouped_matrices
    end 

    function hankelization(A)
    K,L= size(A)
    N = L + K - 1

    time_series = zeros(N)
    weighted_series = zeros(N)

    anti_diagonal_sum = zeros(N)
    counts = zeros(N)
    for i in 1:K, j in 1:L
        n = i + j - 1
        anti_diagonal_sum[n] += A[i,j]

        counts[n] += 1
    end 

    time_series = anti_diagonal_sum ./ counts
    weighted_series = anti_diagonal_sum

    return (time_series, weighted_series)
    end 

    function w_matrix(X, rcs)
    K, L= size(X)
    N = L + K - 1
    
    w = zeros(rcs, N)
    xt = zeros(rcs, N)
    (U, Σ, V) = svd(X)
    
    for i in 1:rcs
        x_rcs = Σ[i] * U[:,i] * V[:,i]'
       (x_rcs_series, x_w_rcs_series) = hankelization(x_rcs)
        
       w[i,:] = x_w_rcs_series
       xt[i,:] = x_rcs_series
    end

    W= w*xt'

    norms = sqrt.(diag(W))
    return Wcorr = W ./ (norms * norms')
    end

end