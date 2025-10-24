"""
VAR(data, p, intercept)

Is a function to estimate the reduced form VAR(p)

Following Killian and Lütkepohl the LS estimator is:
A = [v, A1, A2,..., Ap] = YZ'(ZZ')^-1

Where Y = [y1, y2, ..., yT] is KxT
Z = [Z0,...,ZT-1] is (Kp+1)xTs
Zt-1 = (1, y't-1, ... , y't-p)

"""
function VAR(data::AbstractMatrix, p::Int=1, intercept::Bool=false)
    T, K = size(data)
    Teff = T - p # Number of effective periods

    # Define the Y matrix
    Y = mapreduce(hcat, (p+1):T) do x
        data[x, :]
    end

    Z_array = Array{Float64}(undef, Teff, K * p)
    for i in 1:p
        # columna bloque para el desfase i: y_{t-i}
        Z_array[:, (K*(i-1)+1):(K*i)] = data[(p+1-i):(T-i), :]
    end

    # Define the Z matrix
    if intercept
        Z = hcat(hones(Teff), Z_array)'
    else
        Z = Z_array'
    end

    A = (Z' \ Y')'

    Σu = (Y - A * Z) * (Y - A * Z)' / (T - K * p - 1)

    #return (A = A, Y=Y, Z=Z, Σu = Σu)

    return Dict(
        "A" => A,
        "Y" => Y,
        "Z" => Z,
        "Σu" => Σu
    )
end


"""
get_params()

Funtion that provides a matrix with the estimate coefficients
"""
function get_params(VAR, K::Int, p::Int)
    T = K * p

    get = (VAR, x::Int, y::Int) -> VAR[:, x:y]
    y = K:K:T
    x = y .- (K - 1)

    mats =
        map(x, y) do x, y
            get(VAR, x, y)
        end

    return mats

end

