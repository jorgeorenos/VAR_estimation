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
        "Σu" => Σu,
        "p" => p,
        "K" => K
    )
end


"""
get_params()

Funtion that provides a matrix with the estimate coefficients
"""
function get_params(VAR_est)
    K = VAR_est["K"]
    p = VAR_est["p"]

    T = K * p

    get = (VAR_est, x::Int, y::Int) -> VAR_est[:, x:y]
    y = K:K:T
    x = y .- (K - 1)

    mats =
        map(x, y) do x, y
            get(VAR_est["A"], x, y)
        end

    return mats

end

"""
    BQ function identifys the VAR by the Blanchard an Qua method
    
    we need the Theta matrix: chol(A(1)^1ΣuA(1)^1')
    where A(1) = (I_k - A1 - A2 - ... - Ap)

"""

function BQ_VAR(VAR_est)
    K = VAR_est["K"]
    
    A = get_params(VAR_est)

    A_1 = Matrix(I, K, K)

    for i in 1:length(A)
        A_1 -= A[i]
    end

    A_1_inv = inv(A_1)

    Θ_1 = cholesky(Hermitian(A_1_inv * VAR_est["Σu"] * A_1_inv')).L

    B_0_inv = A_1 * Θ_1
    B_0 = inv(B_0_inv)

    return Dict("B_0_inv" => B_0_inv, "B_0" => B_0)

end

"""
    The IRF function provides the Impulse Response Function following Killian and Lütkepohl

    The function needs the A matrices and the J matrix.

    the argument is the VAR object, a boolean parameter that indicates if we wat the reduced or
    structural IRFs under the Blanchard and Qua method, and the number of periods for the IRFs. 
    The VAR object contains the K variables and p number of lags

    example:

    IRF(VAR, structural = false)
"""

function IRF(VAR_est, periods::Int=20, structural::Bool=false)

    K = VAR_est["K"]
    
    # J Matrix
    J = hcat(Matrix(I, K, K), zeros(K, K*K))
    
    # A Matrix
    A = vcat(
        VAR_est["A"],
        hcat(Matrix(I, K^2, K^2), zeros(K^2, K))
    )

    IRF = Array{Float64}(undef, K, K, periods)

    if structural

        for i in 1:periods
            IRF[:, :, i] = (J * A^i * J') * BQ_VAR(VAR_est)["B_0_inv"]
        end

    else

        for i in 1:periods
            IRF[:, :, i] = (J * A^i * J') 
        end

    end

    return IRF

end