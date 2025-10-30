############# script to testing diferents number of lags ############# 
using DrWatson
@quickactivate "VAR_estimation"

using DataFrames, DataFramesMeta, CSV
using CairoMakie, UnicodePlots, Dates
using Statistics
using TerminalPager
using LinearAlgebra
using PrettyTables

include(srcdir("functions.jl"))

# Load the data
gtdata = CSV.read(
    datadir("data.csv"),
    DataFrame
)
dates = gtdata.date[13:end]
variables = names(gtdata)[2:end]

# Transform in a Matrix
gtdata = (Matrix)(gtdata[:, [:dla_gdp, :dla_cpi, :r]])

# Get the data from 2005Q1 to 2024Q4
gtdata = gtdata[13:end,:]
l = 20

############# plot the data ###############
fig = Figure(size = (1000, 400))

Label(fig[1,1:3], "Data", fontsize = 25)

map(1:size(gtdata, 2)) do d
    ax = Axis(
        fig[2,d],
        title = variables[d],
        xgridvisible = false,
        ygridvisible = false,
        xticks = (1:7:size(gtdata, 1), dates[1:7:end]),
        xticklabelrotation = pi/4
    )

    lines!(
        ax,
        gtdata[:,d]
    )

    hlines!(
        ax,
        0,
        color = :black
    )

end

fig

############# Multiple models ###############
models = Dict()
IRFs = Dict()
SR = Array{Float64}(undef, 4)
SR_log = Array{Float64}(undef, 4)
for p in 1:4
    VAR_est = VAR(gtdata, p)
    models["VAR_$(p)"] = VAR_est 
    IRFs["VAR_$(p)"] = IRF(VAR_est, l, true)

    SR[p] = compute_sacrifice_ratio(
        IRFs["VAR_$(p)"][1,2,:],
        IRFs["VAR_$(p)"][2,2,:],
        20,
        irf_is_gap=false
        )

    SR_log[p] = compute_sacrifice_ratio(
        IRFs["VAR_$(p)"][1,2,:],
        IRFs["VAR_$(p)"][2,2,:],
        20,
        irf_is_gap=true
        )

end

data = [
    "VAR(1)" SR[1] SR_log[1];
    "VAR(2)" SR[2] SR_log[2];
    "VAR(3)" SR[3] SR_log[3];
    "VAR(3)" SR[4] SR_log[4]
]


pretty_table(
    data;
    column_labels = ["Model", "SR gap", "SR log"]
)

# IRFs plots

fig = Figure(size = (900, 600))

Label(
    fig[1,1:3],
    "GDP response to a Inflation shock",
    fontsize = 25
)

map(1:3) do p
    ax = Axis(
        fig[2,p],
        title = "VAR($(p))",
        xgridvisible = false
    )

    IRF = IRFs["VAR_$(p)"]

    lines!(
        ax,
        IRF[1,2,:]
        )

    hlines!(
        ax,
        0
    )
end

Label(fig[3,:], "Inflation response to a inflation shock", fontsize = 25)

map(1:3) do p
    ax = Axis(
        fig[4,p],
        title = "VAR($(p))",
        xgridvisible = false
    )

    IRF = IRFs["VAR_$(p)"]

    lines!(
        ax,
        IRF[2,2,:]
        )

    hlines!(
        ax,
        0
    )
end

fig

data = [
    "VAR(1)" SR[1];
    "VAR(2)" SR[2];
    "VAR(3)" SR[3]        
]

header = ["Model", "sacrifice ratio"]

############ Table of sacrifice ratios ###########
pretty_table(
    data;
    column_labels = ["Model", "sacrifice ratio"]
)
