using DrWatson
@quickactivate "VAR_estimation"

using DataFrames, DataFramesMeta, CSV
using CairoMakie, UnicodePlots, Dates
using Statistics
using TerminalPager
using LinearAlgebra
using PrettyTables

include(srcdir("functions.jl"))

# Define the length of the IRFs
l = 20

# Load the data
gtdata = CSV.read(
    datadir("gtdata.csv"),
    DataFrame
)

# data for total inflation
data_headline_inflation = @chain gtdata begin
    @transform :r = :i .- :d4_ln_cpi
    @select :date :d4_ln_y :d4_ln_cpi :r
end

# transform in a matrix
data_headline_mat = (Matrix)(data_headline_inflation[12:end, 2:end])

# data for MSE inflation
data_MSE_inflation = @chain gtdata begin
    @transform :r = :i .- :d4_ln_cpi_sub
    @select :date :d4_ln_y :d4_ln_cpi_sub :r
end

data_MSE_mat = (Matrix)(data_MSE_inflation[12:end, 2:end])

# data for the headline and nominal


################# Estimations for headline inflation ####################
lags = 4 # Max number of lags
models_headline = Dict()
IRFs_headline = Dict()
SR_cum = Array{Float64}(undef, lags) # Sacrifice ratio by a reducction in the inflation
SR_log = Array{Float64}(undef, lags)
for p in 1:lags
    VAR_est = VAR(data_headline_mat, p)
    
    models_headline["VAR_$(p)"] = VAR_est 
    
    IRFs_headline["VAR_$(p)"] = IRF(VAR_est, l, true)

    SR_cum[p] = compute_sacrifice_ratio(
        IRFs_headline["VAR_$(p)"][1,2,:],
        IRFs_headline["VAR_$(p)"][2,2,:],
        20,
        irf_is_gap=false
        )

    SR_log[p] = compute_sacrifice_ratio(
        IRFs_headline["VAR_$(p)"][1,2,:],
        IRFs_headline["VAR_$(p)"][2,2,:],
        20,
        irf_is_gap=true
        )

end

# Table of sacrifice ratios
pretty_table(
    ["VAR(1)" SR_cum[1] SR_log[1];
    "VAR(2)" SR_cum[2] SR_log[2];
    "VAR(3)" SR_cum[3] SR_log[3];
    "VAR(4)" SR_cum[4] SR_log[4]
    ];
    column_labels = ["Model", "SR cummulative", "SR log GDP"]
)

################# Estimations for headline inflation ####################
lags = 4 # Max number of lags
models_MSE = Dict()
IRFs_MSE = Dict()
SR_cum_MSE = Array{Float64}(undef, lags) # Sacrifice ratio by a reducction in the inflation
SR_log_MSE = Array{Float64}(undef, lags)
for p in 1:lags
    VAR_est = VAR(data_MSE_mat, p)
    
    models_MSE["VAR_$(p)"] = VAR_est 
    
    IRFs_MSE["VAR_$(p)"] = IRF(VAR_est, l, true)

    SR_cum_MSE[p] = compute_sacrifice_ratio(
        IRFs_MSE["VAR_$(p)"][1,2,:],
        IRFs_MSE["VAR_$(p)"][2,2,:],
        20,
        irf_is_gap=false
        )

    SR_log_MSE[p] = compute_sacrifice_ratio(
        IRFs_MSE["VAR_$(p)"][1,2,:],
        IRFs_MSE["VAR_$(p)"][2,2,:],
        20,
        irf_is_gap=true
        )

end

# Table of sacrifice ratios
pretty_table(
    ["VAR(1)" SR_cum_MSE[1] SR_log_MSE[1];
    "VAR(2)" SR_cum_MSE[2] SR_log_MSE[2];
    "VAR(3)" SR_cum_MSE[3] SR_log_MSE[3];
    "VAR(4)" SR_cum_MSE[4] SR_log_MSE[4]
    ];
    column_labels = ["Model", "SR cummulative", "SR log GDP"]
)

########### Testing headline with nominal interes rate ##############
#data_head_nominal_mat = hcat(data_headline_mat[:,1:2], gtdata[12:end, 5])

lags = 4 # Max number of lags
models_nomrate = Dict()
IRFs_nomrate = Dict()
SR_cum_nomrate = Array{Float64}(undef, lags) # Sacrifice ratio by a reducction in the inflation
SR_log_nomrate = Array{Float64}(undef, lags)
for p in 1:lags
    VAR_est = VAR(data_head_nominal_mat, p)
    
    models_nomrate["VAR_$(p)"] = VAR_est 
    
    IRFs_nomrate["VAR_$(p)"] = IRF(VAR_est, l, true)

    SR_cum_nomrate[p] = compute_sacrifice_ratio(
        IRFs_nomrate["VAR_$(p)"][1,2,:],
        IRFs_nomrate["VAR_$(p)"][2,2,:],
        20,
        irf_is_gap=false
        )

    SR_log_nomrate[p] = compute_sacrifice_ratio(
        IRFs_nomrate["VAR_$(p)"][1,2,:],
        IRFs_nomrate["VAR_$(p)"][2,2,:],
        20,
        irf_is_gap=true
        )

end

# Table of sacrifice ratios
pretty_table(
    ["VAR(1)" SR_cum_nomrate[1] SR_log_nomrate[1];
    "VAR(2)" SR_cum_nomrate[2] SR_log_nomrate[2];
    "VAR(3)" SR_cum_nomrate[3] SR_log_nomrate[3];
    "VAR(4)" SR_cum_nomrate[4] SR_log_nomrate[4]
    ];
    column_labels = ["Model", "SR cummulative", "SR log GDP"]
)

UnicodePlots.lineplot(
    data_head_nominal_mat[:,3]
)

UnicodePlots.lineplot(
    data_headline_mat[:,3]
)

