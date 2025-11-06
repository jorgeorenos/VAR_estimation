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
IRFs_VAR_4 = IRF(VAR_4, l, false, true)
sum(cumsum(IRFs_VAR_4[1,2,:])) / sum(IRFs_VAR_4[2,2,:])

UnicodePlots.lineplot(IRFs_VAR_4[1,2,:]*-1)
UnicodePlots.lineplot(IRFs_VAR_4[2,2,:]*-1)

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

# d4_ln_y
δ_y = trend[end] - trend[end-1]

trend_growth = [trend[end]]

for i in 1:20
    x_trend = trend_growth[end] + δ_y
    push!(trend_growth, x_trend)
end

trend_growth_gap = trend_growth[2:end] + forecast[:,1] 
d4_ln_y = d4_ln_fn(vcat(GT_log_data.ln_y, trend_growth_gap))
d4_ln_y = d4_ln_y[13:end]

# Create the data frame with the forecast
forecast_var_4 = vcat(data_hedline_pred_data[9:end,:], forecast[:,:])
forecast_var_4 = hcat(forecast_var_4, d4_ln_y)
forecast_var_4 = hcat(Dates.format.(Date("2005-3"):Month(3):Date("2030-6"), "mm-yyyy"), forecast_var_4)
dates_forecast = forecast_var_4[:,1]

# conditional forecast
CF = forecast_from_structural_shocks(VAR_4, 20; shock_idx=2, shock_size= -1, shock_h=1)

UnicodePlots.lineplot(CF[1,:])
UnicodePlots.lineplot(CF[2,:])

forecast_cf = zeros(size(forecast))
forecast_cf[:,1] = forecast[:,1] + IRFs_VAR_4[1,2,:].*-1
forecast_cf[:,2] = forecast[:,2] + IRFs_VAR_4[2,2,:].*-1
forecast_cf[:,3] = forecast[:,3] + IRFs_VAR_4[3,2,:].*-1

# d4_ln_y
δ_y = trend[end] - trend[end-1]

trend_growth = [trend[end]]

for i in 1:20
    x_trend = trend_growth[end] + δ_y
    push!(trend_growth, x_trend)
end

trend_growth_gap = trend_growth[2:end] + forecast_cf[:,1] 
d4_ln_y = d4_ln_fn(vcat(GT_log_data.ln_y, trend_growth_gap))
d4_ln_y = d4_ln_y[13:end]

# Create the data frame with the forecast
forecast_var_4_cond = vcat(data_hedline_pred_data[9:end,:], forecast_cf[:,:])
forecast_var_4_cond = hcat(forecast_var_4_cond, d4_ln_y)
forecast_var_4_cond = hcat(Dates.format.(Date("2005-3"):Month(3):Date("2030-6"), "mm-yyyy"), forecast_var_4_cond)
dates_forecast_cond = forecast_var_4_cond[:,1]


# output and ouputgap
fig = Figure(size = (900, 600))

# Output gap
ax = Axis(
    fig[1, 1],
    title = "Brecha del producto",
    xgridvisible=false,
    ygridvisible=false,
    xticks = (
        (1:6:length(dates_forecast[65:end]),
        dates_forecast[65:6:end])
    ),
    xticklabelrotation=pi / 4
)

lines!(
    ax,
    forecast_var_4[65:end,2],
    label = "Pronóstico incondicional"
)

lines!(
    ax,
    forecast_var_4_cond[65:end,2],
    label = "Pronóstico condicional"
)

vlines!(
    ax,
    18,
    color = :black,
    linestyle = :dash
)

hlines!(
    ax,
    0,
    color = :black
)

axislegend()

# Output

ax = Axis(
    fig[2, 1],
    title = "Crecimiento económico",
    xgridvisible=false,
    ygridvisible=false,
    xticks = (
        (1:6:length(dates_forecast[65:end]),
        dates_forecast[65:6:end])
    ),
    xticklabelrotation=pi / 4
)


lines!(
    ax,
    forecast_var_4[65:end,5],
    label = "Pronóstico incondicional"
)

lines!(
    ax,
    forecast_var_4_cond[65:end,5],
    label = "Pronóstico condicional"
)


vlines!(
    ax,
    18,
    color = :black,
    linestyle = :dash
)

hlines!(
    ax,
    3.5,
    color=:black,
    label = "Estado Estacionario 3.5"
)

axislegend()

fig
save(
    plotsdir("headline", "conditional forecast headline, output and gap.png"),
    fig,
    px_per_unit = 2.0
)

##### Inflation conditional vs inconditional
fig = Figure(size = (900, 600))

# Output gap
ax = Axis(
    fig[1, 1],
    title = "Inflación",
    xgridvisible=false,
    ygridvisible=false,
    xticks = (
        (1:6:length(dates_forecast[65:end]),
        dates_forecast[65:6:end])
    ),
    xticklabelrotation=pi / 4
)

lines!(
    ax,
    forecast_var_4[65:end,3],
    label = "Pronóstico incondicional"
)

lines!(
    ax,
    forecast_var_4_cond[65:end,3],
    label = "Pronóstico condicional"
)

vlines!(
    ax,
    18,
    color = :black,
    linestyle = :dash
)

hlines!(
    ax,
    0,
    color = :black
)

axislegend()

fig