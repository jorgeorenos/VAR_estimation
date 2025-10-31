using DrWatson
@quickactivate "VAR_estimation"

using DataFrames, DataFramesMeta, CSV
using CairoMakie, UnicodePlots, Dates
using Statistics, StatsBase, Random
using TerminalPager
using LinearAlgebra
using PrettyTables

include(srcdir("functions.jl"))
include(srcdir("helpers.jl"))

# Defining some functions
d4_ln_fn = (x) -> x[5:end] - x[1:(end-4)]

# Define the length of the IRFs
l = 20

# Load the data
GT_log_data = CSV.read(
    datadir("data_log.csv"),
    DataFrame
)

dates = GT_log_data.dates[17:end]

GT_log_data.dates = Date("2001-03"):Month(3):Date("2024-12")
UnicodePlots.lineplot(hodrick_prescott_filter(GT_log_data.ln_y, 1600))
GT_log_data.y_gap = hodrick_prescott_filter(GT_log_data.ln_y, 1600)

# data for headline inflation
data_MSE = @chain GT_log_data begin
    @rsubset :dates >= Date("2003-3")
    @select :dates :y_gap :ln_cpi_sub :i
end

data_MSE_d4_ln = copy(data_MSE)

# Transformations for d4_ln data
data_MSE_d4_ln.d4_ln_cpi_sub = [fill(NaN,4); d4_ln_fn(data_MSE_d4_ln.ln_cpi_sub)]
data_MSE_d4_ln.r = [fill(NaN,4); data_MSE_d4_ln.i[5:end] - data_MSE_d4_ln.d4_ln_cpi_sub[5:end]]

# only d4_ln_cpi data
data_MSE_d4_ln_mat = (Matrix)(data_MSE_d4_ln[data_MSE_d4_ln.dates .>= Date("2005-3"), [:y_gap, :d4_ln_cpi_sub, :r]])

# plot the data
variables = [
    "Bracha del producto",
    "Inflación subyacente optima MSE",
    "Tasa de interés real"
]

fig = Figure(size = (1200, 600))

Label(
    fig[0, 1:3],
    "Datos para Guatemala",
    fontsize = 25
)

map(1:3) do var
    ax = Axis(
        fig[1,var],
        title = variables[var],
        ygridvisible = false,
        xgridvisible = false,
        xticks = (1:7:size(data_MSE_d4_ln_mat, 1), dates[1:7:end]),
        xticklabelrotation = pi/4
    )

    lines!(
        ax,
        data_MSE_d4_ln_mat[:,var],
        color = :black
    )

    hlines!(
        ax,
        0,
        color = :black
    )

end

fig
save(
    plotsdir("sub", "data_sub.png"),
    fig,
    px_per_unit = 2.0
)


##### SVAR model for headline inflation d4_ln data #########
lags = 4 # Max number of lags
models_d4_ln_sub = Dict()
IRFs_d4_ln_sub = Dict()
SR_d4_ln_sub = Array{Float64}(undef, lags) # Sacrifice ratio by a reducction in the inflation
for p in 1:lags
    VAR_est = VAR(data_MSE_d4_ln_mat, p)
    
    models_d4_ln_sub["VAR_$(p)"] = VAR_est 
    
    IRFs_d4_ln_sub["VAR_$(p)"] = IRF(VAR_est, 20, true)

    SR_d4_ln_sub[p] = sum(cumsum(IRFs_d4_ln_sub["VAR_$(p)"][1,2,:]))/sum(IRFs_d4_ln_sub["VAR_$(p)"][2,2,:])
end

SR_d4_ln_sub_data = [
    "VAR(1)" SR_d4_ln_sub[1];
    "VAR(2)" SR_d4_ln_sub[2];
    "VAR(3)" SR_d4_ln_sub[3];
    "VAR(4)" SR_d4_ln_sub[4];
]

pt = pretty_table(
    SR_d4_ln_sub_data;
    column_labels = ["Modelos", "Coeficiente de sacrificio"]
)

# Plot the IRFs #
fig = Figure(size = (1200, 1000))

Label(fig[0, 1:lags], "Impulso Respuesta para cada orden propuesto", fontsize = 27)
Label(fig[1,1:lags], "IRF de la brecha del producto", fontsize = 20)

map(1:lags) do d
    ax = Axis(
        fig[2,d],
        title = "Var$(d)",
        xgridvisible = false,
        ygridvisible = false
    )

    IRF = IRFs_d4_ln_sub["VAR_$(d)"].*-1

    lines!(
        ax,
        IRF[1,2,:]
    )

    hlines!(
        ax,
        0,
        color = :black
    )

end

Label(fig[3,1:lags], "IRF de la inflación total", fontsize = 20)

map(1:lags) do d
    ax = Axis(
        fig[4,d],
        title = "Var$(d)",
        xgridvisible = false,
        ygridvisible = false
    )

    IRF = IRFs_d4_ln_sub["VAR_$(d)"]

    lines!(
        ax,
        IRF[2,2,:].*-1 
    )

    hlines!(
        ax,
        0,
        color = :black
    )

end

Label(fig[5,1:lags], "IRF de la diferencia en la inflación total", fontsize = 20)

map(1:lags) do d
    ax = Axis(
        fig[6,d],
        title = "Var$(d)",
        xgridvisible = false,
        ygridvisible = false
    )

    IRF = IRFs_d4_ln_sub["VAR_$(d)"]

    lines!(
        ax,
        (IRF[2,2,1:end].*-1 - [0; IRF[2,2,1:end-1].*-1])
    )

    hlines!(
        ax,
        0,
        color = :black
    )

end

fig
save(
    plotsdir("sub", "IRFs_inflación_MSE.png"),
    fig,
    ps_per_unit = 2.0
)

# Plot for the VAR(4)
fig = Figure(size = (800, 600))

Label(
    fig[0, 1], 
    "Funciones Impulso Respuesta del modelo VAR estructural",
    fontsize = 25
    )

ax = Axis(
    fig[1,1],
    title = "Respuesta del producto",
    xgridvisible = false,
    ygridvisible = false
)

hidespines!(ax, :r, :t)

lines!(
    ax,
    IRFs_d4_ln_sub["VAR_4"][1,2,:].*-1,
    color = :blue,
    linewidth = 2
)

hlines!(
    ax,
    0,
    color = :black,
    linewidth = 1
)

ax = Axis(
    fig[2,1],
    title = "Respuesta de la inflación",
    xgridvisible = false,
    ygridvisible = false
)

hidespines!(ax, :r, :t)

lines!(
    ax,
    IRFs_d4_ln_sub["VAR_4"][2,2,:].*-1,
    color = :blue,
    linewidth = 2
)

hlines!(
    ax,
    0,
    color = :black,
    linewidth = 1
)

fig
save(
    plotsdir("sub", "IRFs_VAR(4).png"),
    fig,
    px_per_unit = 2.0
)

sum(cumsum(IRFs_d4_ln_sub["VAR_2"][1,3,:]))/sum(IRFs_d4_ln_sub["VAR_2"][2,3,:])


########## Simulation for VAR(4) for y-o-y data ###################
VAR_4 = VAR(data_MSE_d4_ln_mat, 4)
IRFs_4 = IRF(VAR_4, l, true)
SR_point = sum(cumsum(IRFs_4[1,2,:]))/sum(IRFs_4[2,2,:])

# simulation
Random.seed!(1234)
replications = 10000
SR_sim = Array{Float64}(undef, replications)
for j in 1:replications
    VAR_4 = VAR(data_MSE_d4_ln_mat, 4)

    est_values = VAR_4["A"] * VAR_4["Z"]
    errors = VAR_4["Y"] - est_values
    errors = errors'

    errors_boots_array = Matrix{Float64}(undef, 76, 3)
    for var in 1:3
        rand_errors = Array{Float64}(undef, 76)
        for i in 1:76
            indx = rand(1:76)
            rand_errors[i] = errors[indx, var]
        end
        errors_boots_array[:,var] = rand_errors 
    end

    boots_data = est_values' .+ errors_boots_array

    # estimating with the bootstrap data
    VAR_boots = VAR(boots_data, 4)
    IRFs_boots = IRF(VAR_boots, l, true)
    SR_sim[j] = sum(cumsum(IRFs_boots[1,2,:]))/sum(IRFs_boots[2,2,:])
end

SR_clean = remove_outliers(SR_sim, lw = 0.01, up = 0.99)
mean(SR_clean)
median(SR_clean)
std(SR_clean)
skewness(SR_clean)
kurtosis(SR_clean)

fig = Figure(size = (900, 600))

ax = Axis(
    fig[1,1],
    title = "Histograma del coeficientes de sacrifico
    9800 simulaciones",
    xgridvisible = false,
    ygridvisible =false
)

hist!(ax, SR_clean)
vlines!(ax, SR_point, color = :black, label = "Estimación para datos observados $(round(SR_point, digits = 2))")
vlines!(ax, mean(SR_clean), color = :black, linestyle = :dash, label = "simulaciones $(round(mean(SR_clean), digits = 2))")

axislegend()
fig
save(
    plotsdir("sub", "Histograma_simulaciones.png"),
    fig,
    px_per_unit = 2.0
)

quantile(SR_clean, 0.05)
quantile(SR_clean, 0.95)
sum(quantile(SR_clean, 0.05) .< SR_clean .< quantile(SR_clean, 0.95))/length(SR_clean)
