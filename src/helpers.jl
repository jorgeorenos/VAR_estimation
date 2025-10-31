function hodrick_prescott_filter(y, lambda)
    T = size(y, 1)
    matrix = zeros(T, T)

    matrix[1, 1:3] = [1 + lambda, -2 * lambda, lambda]
    matrix[2, 1:4] = [-2 * lambda, 1 + 5 * lambda, -4 * lambda, lambda]

    for i = 3 : T - 2
        matrix[i, i-2 : i+2] = [lambda, -4*lambda, 1 + 6 * lambda, -4 * lambda, lambda]
    end

    matrix[T-1, T-3:T] = [lambda, -4 * lambda, 1 + 5 * lambda, -2 * lambda]
    matrix[T, T-2:T] = [lambda, -2 * lambda, 1 + lambda]

    trend = matrix \ y
    cycle = y - trend

    return cycle
end

function remove_outliers(data; lw = 0.10, up = 0.90)
    q1 = quantile(data, lw)
    q2 = quantile(data, up)
    lower_bound = q1
    upper_bound = q2
    return filter(x -> lower_bound < x < upper_bound, data)
end