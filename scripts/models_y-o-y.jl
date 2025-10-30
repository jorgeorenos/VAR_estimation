############# script to testing models with year on year change data ############# 
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
    datadir("y-o-y-data.csv"),
    DataFrame
)
dates = gtdata.date[12:end]

# Crete a two data sets
# Normal year on year change data
gtdata_yoy = @chain gtdata begin
    @select :d4_ln_y :d4_ln_cpi :r
end
gtdata_yoy = (Matrix)(gtdata_yoy[12:end, :])

# Delta year on year change data
gtdata_delta = @chain gtdata begin
    @select :d4_ln_y :D_d4_ln_cpi :r_D
end
gtdata_delta_test = (Matrix)(gtdata_delta[12:end, :])

############# plot the data ###############
fig = Figure(size = (1000, 400))

Label(fig[1,1:3], "Data", fontsize = 25)
variables = ["GDP", "Inflation", "real interest rate"]
map(1:size(gtdata_yoy, 2)) do d
    ax = Axis(
        fig[2,d],
        title = variables[d],
        xgridvisible = false,
        ygridvisible = false,
        xticks = (1:7:size(gtdata_yoy, 1), dates[1:7:end]),
        xticklabelrotation = pi/4
    )

    lines!(
        ax,
        gtdata_yoy[:,d]
    )

    hlines!(
        ax,
        0,
        color = :black
    )

end

fig
save(
    plotsdir("Normal yoy change data.png"),
    fig,
    px_per_unit = 2.0
)

############## Models for normal data ##############
lags = 4 # Max number of lags
models_normal_data = Dict()
IRFs_normal_data = Dict()
SR_gap = Array{Float64}(undef, lags) # Sacrifice ratio by a reducction in the inflation
SR_log = Array{Float64}(undef, lags)
for p in 1:lags
    VAR_est = VAR(gtdata_yoy, p)
    
    models_normal_data["VAR_$(p)"] = VAR_est 
    
    IRFs_normal_data["VAR_$(p)"] = IRF(VAR_est, l, true)

    SR_gap[p] = compute_sacrifice_ratio(
        IRFs_normal_data["VAR_$(p)"][1,2,:],
        IRFs_normal_data["VAR_$(p)"][2,2,:],
        20,
        irf_is_gap=false
        )

    SR_log[p] = compute_sacrifice_ratio(
        IRFs_normal_data["VAR_$(p)"][1,2,:],
        IRFs_normal_data["VAR_$(p)"][2,2,:],
        20,
        irf_is_gap=true
        )

end


# Table of sacrifice ratios
pretty_table(
    ["VAR(1)" SR_gap[1] SR_log[1];
    "VAR(2)" SR_gap[2] SR_log[2];
    "VAR(3)" SR_gap[3] SR_log[3];
    "VAR(4)" SR_gap[4] SR_log[4]
    ];
    column_labels = ["Model", "SR gap", "SR log GDP"]
)

