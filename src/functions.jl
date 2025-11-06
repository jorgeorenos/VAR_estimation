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
        Z = hcat(ones(Teff), Z_array)'
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

    # Obtain the structural matrices


    return Dict("B_0_inv" => B_0_inv, "B_0" => B_0, "A_1" => A_1)

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

function IRF(VAR_est, periods::Int=20, constant::Bool=false, structural::Bool=false)

    if constant == false
        K = VAR_est["K"]
        p = VAR_est["p"]

        # J Matrix
        J = hcat(Matrix(I, K, K), zeros(K, K*(p-1)))
        
        # A Matrix
        A = vcat(
            VAR_est["A"],
            hcat(Matrix(I, K*(p-1), K*(p-1)), zeros(K*(p-1), K))
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

    else
        K = VAR_est["K"]
        p = VAR_est["p"]

        # J Matrix
        J = hcat(Matrix(I, K, K), zeros(K, K*(p-1)))
        
        # A Matrix
        A = vcat(
            VAR_est["A"][:,2:end],
            hcat(Matrix(I, K*(p-1), K*(p-1)), zeros(K*(p-1), K))
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
    end

    return IRF

end

function Forecast_pre_data(VAR_est, periods, data_proj, constant::Bool=false)
    
    K = VAR_est["K"]

    if constant == false
    
        # Projections without constant
        
        A = VAR_est["A"]

        proj_matrix = Array{Float64}(undef, periods, K)

        projection = Float64[]

        for i in 1:periods
            x_proj = A*data_proj

            # reorganize the projection vector

            data_proj = vcat(x_proj, data_proj[1:end-K])
            
            push!(projection, x_proj...)

        end

        for j in 1:K
            proj_matrix[:,j] = projection[j:K:end]
        end

    else

        # Projections with constant

        A = VAR_est["A"][:,2:end]
        C = VAR_est["A"][:,1]
        
        proj_matrix = Array{Float64}(undef, periods, K)

        projection = Float64[]

        for i in 1:periods
            x_proj = A*data_proj + C

            # reorganize the projection vector

            data_proj = vcat(x_proj, data_proj[1:end-K])
            
            push!(projection, x_proj...)

        end

        for j in 1:K
            proj_matrix[:,j] = projection[j:K:end]
        end

    end

    return proj_matrix

end

"""
    theta_sequence(VAR_est, H) -> Thetas

Devuelve un Vector{Matrix{Float64}} con Θ₀,…,Θ_{H-1}, donde
Θᵢ = Φᵢ * B0_inv y Φᵢ = J * A_c^i * J'.
No usa observaciones pasadas.
"""
function theta_sequence(VAR_est, H::Int)
    K = VAR_est["K"]; p = VAR_est["p"]

    # J matrix
    J = hcat(Matrix(I, K, K), zeros(K, K*(p-1)))

    # Companion matrix
    A_top = VAR_est["A"]                                
    A_bottom = hcat(Matrix(I, K*(p-1), K*(p-1)), zeros(K*(p-1), K))
    A_c = vcat(A_top, A_bottom)

    # structural contemporaneous impact matrix
    B0_inv = BQ_VAR(VAR_est)["B_0_inv"]                         

    Thetas = Vector{Matrix{Float64}}(undef, H)

    # Initial P = A_c^0
    P = Matrix(I, size(A_c,1), size(A_c,2))
    for i in 0:H-1
        Φi = J * A_c^i * J'                                 
        Thetas[i+1] = Φi * B0_inv                      
                                          
    end
    return Thetas
end

"""
    forecast_from_structural_shocks(VAR_est, H; shock_path=nothing,
                                    shock_idx=1, shock_size=-1.0, shock_h=1)

Pronóstico condicional "MA puro":
ŷ_{t+h|t} = ∑_{i=0}^{h-1} Θᵢ w_{t+h-i}, h=1..H.

- Si no pasas `shock_path`, aplica un único choque estructural de tamaño
  `shock_size` en el índice `shock_idx` en el paso `shock_h`.
- Si pasas `shock_path`, debe ser una matriz K×H con los choques futuros w_{t+1..t+H}.
"""
function forecast_from_structural_shocks(VAR_est, H::Int;
                                         shock_idx::Int=1, shock_size::Real=-1.0, shock_h::Int=1)

    K = VAR_est["K"]
    Θ = theta_sequence(VAR_est, H)                      # Θ₀..Θ_{H-1}

    # Camino de choques futuros (estructurales, var=I)
    W = zeros(K, H) 
    W[shock_idx, shock_h] = shock_size

    # Construye pronóstico condicional sin baseline (MA puro)
    Ycond = zeros(K, H)
    for h in 1:H
        acc = zeros(K)
        for i in 0:h-1
            acc += Θ[i+1] * W[:, h-i]
        end
        Ycond[:, h] = acc
    end
    return Ycond   # columnas: y_{t+1|t}, …, y_{t+H|t}
end

