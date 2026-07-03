using Main.automatic_ssa

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

smooth, residual = automatic_ssa.auto_ssa(x,y,200)
plot!(smooth)
plot(residual)