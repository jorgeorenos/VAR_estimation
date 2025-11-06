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

    return IRF

end

function Forecast(VAR_est, periods, constant::Bool=false)
    
    K = VAR_est["K"]

    if constant == false
    
        # Projections without constant
        data_proj = VAR_est["Z"][:,end]
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

        data_proj = VAR_est["Z"][:,end][2:end]
        A = VAR_est["A"][:,2:end]
        C = VAR_est["A"][:,1]
        
        proj_matrix = Array{Float64}(undef, periods, K)

        projection = Float64[]

        for i in 1:periods
            x_proj = A*data_proj

            # reorganize the projection vector

            data_proj = vcat(x_proj, data_proj[1:end-K])
            
            push!(projection, x_proj...)

        end

        for j in 1:K
            proj_matrix[:,j] = projection[j:K:end].+C[j]
        end

    end

    return proj_matrix

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
            x_proj = A*data_proj

            # reorganize the projection vector

            data_proj = vcat(x_proj, data_proj[1:end-K])
            
            push!(projection, x_proj...)

        end

        for j in 1:K
            proj_matrix[:,j] = projection[j:K:end].+C[j]
        end

    end

    return proj_matrix

end

function ConditionalForecast(VAR_est, periods::Int=20)

    K = VAR_est["K"]
    p = VAR_est["p"]

    # J Matrix
    J = hcat(Matrix(I, K, K), zeros(K, K*(p-1)))
    
    # A Matrix
    A = vcat(
        VAR_est["A"],
        hcat(Matrix(I, K*(p-1), K*(p-1)), zeros(K*(p-1), K))
    )

    cond_forecast = Array{Float64}(undef, K, K, periods)

    for i in 1:periods

        if i == 1
                
                cond_forecast[:, :, i] = (J * A^i * J') * BQ_VAR(VAR_est)["B_0_inv"]
                
        else

                cond_forecast[:, :, i] = (J * A^i * J') * BQ_VAR(VAR_est)["B_0_inv"]

        end

    return cond_forecast


    end
end

# ===== Utilidades para Θᵢ y pronóstico condicional =====

"""
    theta_sequence(VAR_est, H) -> Thetas

Devuelve un Vector{Matrix{Float64}} con Θ₀,…,Θ_{H-1}, donde
Θᵢ = Φᵢ * B0_inv y Φᵢ = J * A_c^i * J'.
"""
function theta_sequence(VAR_est, H::Int)
    K = VAR_est["K"]; p = VAR_est["p"]
    J = hcat(Matrix(I, K, K), zeros(K, K*(p-1)))
    A_top = VAR_est["A"]                                # K × (Kp)
    A_bottom = hcat(Matrix(I, K*(p-1), K*(p-1)), zeros(K*(p-1), K))
    A_c = vcat(A_top, A_bottom)                         # Kp × Kp
    B0_inv = BQ_VAR(VAR_est)["B_0_inv"]                         # K × K

    Thetas = Vector{Matrix{Float64}}(undef, H)
    P = Matrix(I, size(A_c,1), size(A_c,2))             # A_c^0
    for i in 0:H-1
        Φi = J * P * J'                                 # K × K
        Thetas[i+1] = Φi * B0_inv                      # Θᵢ
        P = A_c * P                                     # eleva potencia
    end
    return Thetas
end

"""
    baseline_forecast(VAR_est, ylags, H) -> Ybase

Pronóstico incondicional (sin choques futuros) usando la forma compañera.
Si tus datos están centrados y c=0, este baseline puede ser cercano a 0.
Retorna matriz K×H con [y_{T+1|T} … y_{T+H|T}].
"""
function baseline_forecast(VAR_est, ylags, H::Int)
    K = VAR_est["K"]; p = VAR_est["p"]
    J = hcat(Matrix(I, K, K), zeros(K, K*(p-1)))
    A_top = VAR_est["A"]
    A_bottom = hcat(Matrix(I, K*(p-1), K*(p-1)), zeros(K*(p-1), K))
    A_c = vcat(A_top, A_bottom)
    c = haskey(VAR_est, "c") ? VAR_est["c"] : zeros(K)
    d = vcat(c, zeros(K*(p-1)))

    # estado inicial s0 = [y_T; y_{T-1}; …; y_{T-p+1}]
    @assert size(ylags,1)==K && size(ylags,2)==p
    s = vcat((ylags[:,i] for i in 1:p)...)              # Kp×1

    Ybase = Array{Float64}(undef, K, H)
    for h in 1:H
        s = A_c*s + d                                   # sin choques
        Ybase[:,h] = J*s
    end
    return Ybase
end

"""
    conditional_forecast_with_Theta(VAR_est, ylags, H; shock_path, shock_idx, shock_size, shock_h)

Construye el pronóstico condicional:
Ycond = Ybase + ∑_{i=0}^{h-1} Θᵢ * w_{h-i}.
- `shock_path` (opcional): matriz K×H con los choques estructurales futuros.
- Si no pasas `shock_path`, aplica un único choque de tamaño `shock_size`
  en el shock `shock_idx` en el paso `shock_h` (por defecto: -1 en h=1).
"""
function conditional_forecast_with_Theta(VAR_est, ylags, H::Int;
                                         shock_path::Union{Nothing,AbstractMatrix}=nothing,
                                         shock_idx::Int=1, shock_size::Real=-1.0, shock_h::Int=1)

    K = VAR_est["K"]
    # 1) Θ₀…Θ_{H-1}
    Θ = theta_sequence(VAR_est, H)

    # 2) baseline incondicional (sin choques futuros)
    Ybase = baseline_forecast(VAR_est, ylags, H)

    # 3) camino de choques futuros w_{t+1}…w_{t+H}
    W = shock_path === nothing ? zeros(K, H) : copy(shock_path)
    if shock_path === nothing && 1 <= shock_h <= H
        W[shock_idx, shock_h] = shock_size
    end

    # 4) término condicionado: para cada h, suma Θ₀ w_h + Θ₁ w_{h-1} + ...
    Yadd = zeros(K, H)
    for h in 1:H
        acc = zeros(K)
        for i in 0:h-1
            acc += Θ[i+1] * W[:, h-i]
        end
        Yadd[:,h] = acc
    end

    return Ybase + Yadd   # K×H  (columnas: T+1 … T+H)
end

