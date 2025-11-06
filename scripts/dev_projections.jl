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

dates = GT_log_data.dates[17:end]

GT_log_data.dates = Date("2001-03"):Month(3):Date("2025-6")
x = hodrick_prescott_filter(GT_log_data.ln_y, 1600)
GT_log_data.y_gap = x[1]
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

########## Projections for VAR(4) ############
VAR_no_const = VAR(data_headline_d4_ln_mat, 2)
VAR_const = VAR(data_headline_d4_ln_mat, 2, true)


# Creating a function to create projections without constant

data_proj = VAR_no_const["Z"][:,end]

projection = Float64[]

for i in 1:40

    x_proj = VAR_no_const["A"]*data_proj

    # reorganize the projection vector

    data_proj = vcat(x_proj, data_proj[1:end-3])
    
    push!(projection, x_proj...)

end

UnicodePlots.lineplot(projection[3:3:end])

# Creating a function to create projections without constant

data_proj_const = VAR_const["Z"][:,end][2:end]

projection_const = Float64[]

for i in 1:40

    x_proj = VAR_const["A"][:,2:end]*data_proj_const

    # reorganize the projection vector

    data_proj_const = vcat(x_proj, data_proj_const[1:end-3])
    
    push!(projection_const, x_proj...)

end

UnicodePlots.lineplot(projection_const[3:3:end].+VAR_const["A"][2,1])


var_no_const_projections = Forecast(VAR_no_const, 40)
var_const_projections = Forecast(VAR_const, 40, true)

# Test manual process againts function
# both process give the same results
hcat(projection[2:3:end], var_no_const_projections[:,2])
hcat(projection_const[3:3:end].+VAR_const["A"][3,1], var_const_projections[:,3])

UnicodePlots.lineplot(StructuralForecast(VAR_no_const, 5)[:,1])

# Plot the results
# first estimate all models and estimate their projections
# non constants models
lags = 4 # Max number of lags
periods = 20
models = Dict()
forecast = Array{Float64}(undef, periods, 3, lags)

for p in 1:lags
    VAR_est = VAR(data_headline_d4_ln_mat, p)

    models["VAR_$(p)"] = VAR_est

    forecast[:,:,p] = Forecast(VAR_est, periods)

end

# Plot
fig = Figure(size=(1000, 800))

Label(fig[0, 1:lags], "Pronósticos Incondicionales de cada orden propuesto", fontsize=27)
Label(fig[1, 1:lags], "Brecha del producto", fontsize=20)

map(1:lags) do d

    ax = Axis(
        fig[2, d],
        title="Var$(d)",
        xgridvisible=false,
        ygridvisible=false,
        xticks = (
            (1:6:53,
            Dates.format.(Date("2022-3"):Month(3):Date("2035-3"), "mm-yyyy")[1:6:end])
        ),
        xticklabelrotation=pi / 4
    )

    y = forecast[:,1,d]
    y = vcat(data_headline_d4_ln_mat[70:end,1], y)

    lines!(
        ax,
        y
    )

    vlines!(
        ax,
        10,
        color = :black,
        linestyle = :dash
    )

    hlines!(
        ax,
        0,
        color=:black
    )

end

Label(fig[3, 1:lags], "Inflación", fontsize=20)

map(1:lags) do d
    ax = Axis(
        fig[4, d],
        title="Var$(d)",
        xgridvisible=false,
        ygridvisible=false,
        xticks = (
            (1:6:63,
            Dates.format.(Date("2019-12"):Month(3):Date("2035-3"), "mm-yyyy")[1:6:end])
        ),
        xticklabelrotation=pi / 4
    )

    y = forecast[:,2,d]
    y = vcat(data_headline_d4_ln_mat[60:end,2], y)

    lines!(
        ax,
        y
    )

    vlines!(
        ax,
        20,
        color = :black,
        linestyle = :dash
    )

    hlines!(
        ax,
        
    )

end

Label(fig[5, 1:lags], "Tasa de interés real", fontsize=20)

map(1:lags) do d
    ax = Axis(
        fig[6, d],
        title="Var$(d)",
        xgridvisible=false,
        ygridvisible=false,
        xticks = (
            (1:6:63,
            Dates.format.(Date("2019-12"):Month(3):Date("2035-3"), "mm-yyyy")[1:6:end])
        ),
        xticklabelrotation=pi / 4
    )

    y = forecast[:,3,d]
    y = vcat(data_headline_d4_ln_mat[60:end,3], y)

    lines!(
        ax,
        y
    )

    vlines!(
        ax,
        20,
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
    plotsdir("headline", "forecast.png"),
    fig,
    px_per_unit=2.0
)

# unconditional forecast for VAR(4)
# d4_ln_y
ln_y = trend[end] .+ forecast[:,1,4]

δ_y = trend[end] - trend[end-1]

trend_growth = [trend[end]]

for i in 1:20
    x_trend = trend_growth[end] + δ_y
    push!(trend_growth, x_trend)
end

trend_growth_gap = trend_growth[2:end] + forecast[:,1,4] 

UnicodePlots.lineplot(d4_ln_fn(vcat(trend_growth_gap)))
d4_ln_y = d4_ln_fn(vcat(GT_log_data.ln_y, trend_growth_gap))
d4_ln_y = d4_ln_y[13:end] 

# Create the data frame with the forecast
forecast_var_4 = vcat(data_headline_d4_ln_mat, forecast[:,:,4])
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
        17,
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
        17,
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
        17,
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
        17,
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
    plotsdir("headline", "forecast VAR4.png"),
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
    17,
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
    17,
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
    plotsdir("headline", "forecast headline, output and gap.png"),
    fig,
    px_per_unit = 2.0
)

