module kneedle_algo

export kneedle
function normalize_curve(x,y,x_min,x_max,y_min,y_max)
    x_n = (x .- x_min) ./ (x_max - x_min)
    y_n = (y .- y_min) ./ (y_max - y_min) 
    return (x_n,y_n)
end

function difference_curve(x_n, y_n, curve)

    if curve == :knee
        y_d = (y_n) .- x_n
    elseif curve== :elbow
        y_d = (1 .- y_n) .- x_n
    end
    return(y_d)
end

function find_maxima(y_d)
    maximas = Int[]

    for i in 2:length(y_d)-1
        if y_d[i] > y_d[i-1] && y_d[i] > y_d[i+1]
            push!(maximas, i)
        end
    end

    if isempty(maximas)
        maximas = [argmax(y_d)]
    end
    return maximas
end

function sensitivity_filter(maximas, y_d, N, S)
    x_avg = 1 / (N - 1)

    if length(maximas) == 1
        return maximas[1]
    end

    for i in 1:length(maximas)-1
        m = maximas[i]
        next_m = maximas[i+1]

        threshold = y_d[m] - S * x_avg

        for j in m+1:next_m
            if y_d[j] < threshold
                return m
            end
        end
    end

    return maximas[argmax(y_d[maximas])]
end

function kneedle(x, y; S=1.0, curve=:elbow)
    x_min, x_max = minimum(x), maximum(x)
    y_min, y_max = minimum(y), maximum(y)

    N = length(x)

    x_n, y_n = normalize_curve(x, y, x_min, x_max, y_min, y_max)
    y_d = difference_curve(x_n, y_n, curve)

    maximas_idx = find_maxima(y_d)
    knee_idx = sensitivity_filter(maximas_idx, y_d, N, S)

    return knee_idx, x[knee_idx], y[knee_idx]
end

end