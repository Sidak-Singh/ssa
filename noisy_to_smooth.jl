using Plots
plotlyjs
using LinearAlgebra
using Statistics

using Main.SSA
function smooth_function(x)
    y=0.2*x + sin.(2x)
end


function gaussian_noise(σ,signal)
    signal_dim=length(signal)
    noise = σ .* randn(signal_dim)
    return noise
end

function data_mat(N, y)
    y_matrix = zeros()
end

x = range(-20,20,100);
y = smooth_function(x)
plot(x, y)

σ = 1
noise = gaussian_noise(σ,x)

noisy_y = y + noise
plot!(x,noisy_y)


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


X = hankel(noisy_y, 20)
function svd_vals(X)
    (U, Σ, V) = svd(X)
end




(U, Σ, V) = svd(X)



plot(1:20, Σ, xlims=(1,5))
plot(1:20, V[:,1])
plot!(1:20, V[:,2])
plot!(1:20, V[:,3])
plot!(1:20, V[:,4])

function grouped_svd(U, Σ, V, groups)
L = size(U, 1)
K = size(V, 1)

n=length(groups)
grouped_matrices = [zeros(L,K) for i in 1:n]

for i in 1:n 
    for j in groups[i]
        grouped_matrices[i] += Σ[j] * U[:,j] * V[:,j]'
    end 
end
return grouped_matrices
end 

function hankelization(A)
    L, K = size(A)
    N = L + K - 1

    time_series = zeros(N)

    anti_diagonal_sum = zeros(N)
    counts = zeros(N)
    for i in 1:L, j in 1:K
        n = i + j - 1
        anti_diagonal_sum[n] += A[i,j]

        counts[n] += 1
    end 

    time_series = anti_diagonal_sum ./ counts
end 

grouped_matrices = grouped_svd(U, Σ, V, [1,2:3,1:4, 1:2])

trend = hankelization(grouped_matrices[1])
periods = hankelization(grouped_matrices[2])
recon1 = hankelization(grouped_matrices[3])
recon = hankelization(grouped_matrices[4])
plot(1:100, trend)
plot(1:100, periods, xlims=(0,4π))
plot(1:100, recon1)
plot(1:100, recon)

W=w_matrix(X,10)