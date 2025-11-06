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
VAR_4 = VAR(data_headline_d4_ln_mat, 4)
structural = BQ_VAR(VAR_4)

# pre data
data_headline_d4_ln.y_gap[end-1:end] = x2[1][end-1:end]
data_hedline_pred_data = (Matrix)(data_headline_d4_ln[:, [:y_gap, :d4_ln_cpi, :r]])
pre_data = vcat(
    data_hedline_pred_data[end, :],
    data_hedline_pred_data[end-1, :],
    data_hedline_pred_data[end-2, :],
    data_hedline_pred_data[end-3, :],
)

forecast = Forecast_pre_data(VAR_4, 20, pre_data)

# d4_ln_y
δ_y = trend[end] - trend[end-1]

trend_growth = [trend[end]]

for i in 1:20
    x_trend = trend_growth[end] + δ_y
    push!(trend_growth, x_trend)
end

trend_growth_gap = trend_growth[2:end] + forecast[:,1] 

UnicodePlots.lineplot(d4_ln_fn(vcat(trend_growth_gap)))
d4_ln_y = d4_ln_fn(vcat(GT_log_data.ln_y, trend_growth_gap))
d4_ln_y = d4_ln_y[13:end]

# Create the data frame with the forecast
forecast_var_4 = vcat(data_hedline_pred_data[9:end,:], forecast[:,:])
forecast_var_4 = hcat(forecast_var_4, d4_ln_y)
forecast_var_4 = hcat(Dates.format.(Date("2005-3"):Month(3):Date("2030-6"), "mm-yyyy"), forecast_var_4)
dates_forecast = forecast_var_4[:,1] 

# Plot all variables
fig = Figure(size=(900, 1000))

Label(fig[0, 1], "Pronósticos Incondicionales VAR orden 4", fontsize=27)
Label(fig[1, 1], "Brecha del producto", fontsize=20)

map(1:1) do d

    ax = Axis(
        fig[2, d],
        xgridvisible=false,
        ygridvisible=false,
        xticks = (
            (1:6:length(dates_forecast[65:end]),
            dates_forecast[65:6:end])
        ),
        xticklabelrotation=pi / 4
    )

    y = forecast_var_4[65:end,2]

    lines!(
        ax,
        y
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
        color=:black
    )

end

Label(fig[3, 1], "Crecimiento económico", fontsize=20)

map(1:1) do d

    ax = Axis(
        fig[4, d],
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
        forecast_var_4[65:end,5]
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
        color=:black
    )

end

Label(fig[5, 1], "Inflación", fontsize=20)

map(1:1) do d
    ax = Axis(
        fig[6, d],
        xgridvisible=false,
        ygridvisible=false,
        xticks = (
            (1:6:length(dates_forecast[65:end]),
            dates_forecast[65:6:end])
        ),
        xticklabelrotation=pi / 4
    )

    y = forecast_var_4[65:end,3]

    lines!(
        ax,
        y
    )

    vlines!(
        ax,
        18,
        color = :black,
        linestyle = :dash
    )

end

Label(fig[7, 1], "Tasa de interés real", fontsize=20)

map(1:1) do d
    ax = Axis(
        fig[8, d],
        xgridvisible=false,
        ygridvisible=false,
        xticks = (
            (1:6:length(dates_forecast[65:end]),
            dates_forecast[65:6:end])
        ),
        xticklabelrotation=pi / 4
    )

    y = forecast_var_4[65:end,4]

    lines!(
        ax,
        y
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
        color=:black
    )

end

fig
save(
    plotsdir("headline", "unconditional forecast VAR4.png"),
    fig,
    px_per_unit=2.0
)


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

y = forecast_var_4[65:end,2]

lines!(
    ax,
    y
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
    forecast_var_4[65:end,5]
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
    plotsdir("headline", "unconditional forecast headline, output and gap.png"),
    fig,
    px_per_unit = 2.0
)



##################### conditional forecast ####################
pre_data_2 = hcat(
    data_hedline_pred_data[end, :],
    data_hedline_pred_data[end-1, :],
    data_hedline_pred_data[end-2, :],
    data_hedline_pred_data[end-3, :],
)
CF = conditional_forecast_with_Theta(VAR_4, pre_data_2, 20; shock_idx=2, shock_size=-2.0, shock_h=1)

# d4_ln_y reconstrucción
δ_y = trend[end] - trend[end-1]

trend_growth = [trend[end]]

for i in 1:20
    x_trend = trend_growth[end] + δ_y
    push!(trend_growth, x_trend)
end

trend_growth_gap = trend_growth[2:end] + CF[1,:] 

UnicodePlots.lineplot(d4_ln_fn(vcat(trend_growth_gap)))
d4_ln_y = d4_ln_fn(vcat(GT_log_data.ln_y, trend_growth_gap))
d4_ln_y = d4_ln_y[13:end]

# Create the data frame with the forecast
forecast_var_4_cond = vcat(data_hedline_pred_data[9:end,:], CF')
forecast_var_4_cond = hcat(forecast_var_4_cond, d4_ln_y)
forecast_var_4_cond = hcat(Dates.format.(Date("2005-3"):Month(3):Date("2030-6"), "mm-yyyy"), forecast_var_4_cond)
dates_forecast = forecast_var_4[:,1] 


# Plot all variables
fig = Figure(size=(900, 1000))

Label(fig[0, 1], "Pronósticos Condicionales VAR orden 4", fontsize=27)
Label(fig[1, 1], "Brecha del producto", fontsize=20)

map(1:1) do d

    ax = Axis(
        fig[2, d],
        xgridvisible=false,
        ygridvisible=false,
        xticks = (
            (1:6:length(dates_forecast[65:end]),
            dates_forecast[65:6:end])
        ),
        xticklabelrotation=pi / 4
    )

    y = forecast_var_4_cond[65:end,2]

    lines!(
        ax,
        y
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
        color=:black
    )

end

Label(fig[3, 1], "Crecimiento económico", fontsize=20)

map(1:1) do d

    ax = Axis(
        fig[4, d],
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
        forecast_var_4_cond[65:end,5]
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
        color=:black
    )

end

Label(fig[5, 1], "Inflación", fontsize=20)

map(1:1) do d
    ax = Axis(
        fig[6, d],
        xgridvisible=false,
        ygridvisible=false,
        xticks = (
            (1:6:length(dates_forecast[65:end]),
            dates_forecast[65:6:end])
        ),
        xticklabelrotation=pi / 4
    )

    y = forecast_var_4_cond[65:end,3]

    lines!(
        ax,
        y
    )

    vlines!(
        ax,
        18,
        color = :black,
        linestyle = :dash
    )

end

Label(fig[7, 1], "Tasa de interés real", fontsize=20)

map(1:1) do d
    ax = Axis(
        fig[8, d],
        xgridvisible=false,
        ygridvisible=false,
        xticks = (
            (1:6:length(dates_forecast[65:end]),
            dates_forecast[65:6:end])
        ),
        xticklabelrotation=pi / 4
    )

    y = forecast_var_4_cond[65:end,4]

    lines!(
        ax,
        y
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
        color=:black
    )

end

fig
save(
    plotsdir("headline", "conditional forecast VAR4.png"),
    fig,
    px_per_unit=2.0
)

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

y = forecast_var_4_cond[65:end,2]

lines!(
    ax,
    y
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
    forecast_var_4_cond[65:end,5]
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

ax = Axis(
    fig[1,1]
)

lines!(
    ax,
    forecast_var_4[65:end, 3],
    label = "Incondicional"
)

lines!(
    ax,
    forecast_var_4_cond[65:end, 3],
    label = "Condicional"
)

axislegend()

fig
