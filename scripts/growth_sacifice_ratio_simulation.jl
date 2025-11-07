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
first_difference = (x) -> x[2:end] - x[1:end-1]
d4_ln_fn = (x) -> x[5:end] - x[1:(end-4)]
dl_ln_fn = (x) -> x[2:end] - x[1:(end-1)]

# Define the length of the IRFs
l = 20

# Load the data
GT_log_data = CSV.read(
    datadir("data_log.csv"),
    DataFrame
)

GT_log_data.dates = Date("2001-03"):Month(3):Date("2025-6")
x = hodrick_prescott_filter(GT_log_data.ln_y[1:end-2], 1600)
x2 = hodrick_prescott_filter(GT_log_data.ln_y[1:end-2], 1600)
trend = x2[2]
GT_log_data.y_gap = [x[1]; fill(NaN, 2)]
trend = x[2]

# data for headline inflation
data_headline = @chain GT_log_data begin
    @rsubset :dates >= Date("2003-3")
    @select :dates :y_gap :ln_cpi :i
end

data_headline_d4_ln = copy(data_headline)

# Transformations for d4_ln data
data_headline_d4_ln.d4_ln_cpi = [fill(NaN, 4); d4_ln_fn(data_headline_d4_ln.ln_cpi)]
data_headline_d4_ln.r = [fill(NaN, 4); data_headline_d4_ln.i[5:end] - data_headline_d4_ln.d4_ln_cpi[5:end]]

# only d4_ln_cpi data
data_headline_d4_ln_mat = (Matrix)(data_headline_d4_ln[data_headline_d4_ln.dates.>=Date("2005-3"), [:y_gap, :d4_ln_cpi, :r]])
data_headline_d4_ln_mat = data_headline_d4_ln_mat[1:end-2,:]

########## Projections for VAR(4) ############
VAR_4 = VAR(data_headline_d4_ln_mat, 4, false)

# pre data
data_headline_d4_ln.y_gap[end-1:end] = x2[1][end-1:end]
data_hedline_pred_data = (Matrix)(data_headline_d4_ln[:, [:y_gap, :d4_ln_cpi, :r]])
pre_data = vcat(
    data_hedline_pred_data[end, :],
    data_hedline_pred_data[end-1, :],
    data_hedline_pred_data[end-2, :],
    data_hedline_pred_data[end-3, :],
)

forecast = Forecast_pre_data(VAR_4, 20, pre_data, false)

# simulation
Random.seed!(1234)
replications = 10000
forecast_sim = Array{Float64}(undef, 20, 4,replications)
forecast_cf_sim = Array{Float64}(undef, 20, 4,replications)
CF_sim = Array{Float64}(undef, 20, 3, replications)
for j in 1:replications
    VAR_4 = VAR(data_headline_d4_ln_mat, 4)

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
        errors_boots_array[:, var] = rand_errors
    end

    boots_data = est_values' .+ errors_boots_array

    # pre data

    pre_data_boots = vcat(
        boots_data[end, :],
        boots_data[end-1, :],
        boots_data[end-2, :],
        boots_data[end-3, :],
    )


    # estimating with the bootstrap data
    VAR_boots = VAR(boots_data, 4)

    forecast = Forecast_pre_data(VAR_boots, 20, pre_data_boots, false)
    forecast_sim[:,1:3,j] = forecast

    CF = forecast_from_structural_shocks(VAR_4, 20; shock_idx=2, shock_size= -1, shock_h=1)
    CF_sim[:, :, j] = CF'
    forecast_cf = zeros(size(forecast))
    forecast_cf = forecast + CF'
    forecast_cf_sim[:,1:3,j] = forecast_cf
end

# d4_ln_y
δ_y = trend[end] - trend[end-1]

trend_growth = [trend[end]]

for i in 1:20
    x_trend = trend_growth[end] + δ_y
    push!(trend_growth, x_trend)
end

for i in 1:replications
    trend_growth_gap = trend_growth[2:end] + forecast_cf_sim[:,1,i] 
    d4_ln_y = d4_ln_fn(vcat(GT_log_data.ln_y, trend_growth_gap))
    d4_ln_y = d4_ln_y[95:end]
    forecast_cf_sim[:,4,i] = d4_ln_y
end

############# Coeficiente de sacrificio ##############
SR_sim = Array{Float64}(undef, replications)
for j in 1:replications
    SR_sim[j] = sum(cumsum(forecast_cf_sim[:,4,j] .- 3.51)) / sum(CF_sim[:,2,j])
end

# Create an histogram of SR_sim
fig = Figure(size = (900, 600))

ax = Axis(
    fig[1, 1],
    title = "Coeficiente de sacrificio simulaciones",
    titlesize = 23,
    subtitle = "10000 simulaciones",
    xgridvisible=false,
    ygridvisible=false,
    xticksize = 20,
)

hidespines!(
    ax,
    :r,
    :t
)

hist!(ax, SR_sim)

text!(
    ax,
    6,
    2000,
    text="Intervalo de confianza (10%-90%): \n[1.2596%, 4.4674%]",
    fontsize=17
)

fig
save(
    plotsdir("headline", "SR output growth simulations.png"),
    fig,
    px_per_unit = 2.0
)

quantiles_SR = quantile(SR_sim, [0.1, 0.9])

# Plot multiple trajectories of d4_ln_y
fig = Figure(size = (900, 600))
ax = Axis(
    fig[1, 1],
    title = "Trayectorias simuladas de crecimiento del producto",
    titlesize = 23,
    xgridvisible=false,
    ygridvisible=false,
    xticksize = 20,
)
hidespines!(ax, :r, :t)
for j in 1:100
    lines!(
        ax,
        forecast_cf_sim[:,4,j],
        color = (:blue, 0.2)
    )
end


