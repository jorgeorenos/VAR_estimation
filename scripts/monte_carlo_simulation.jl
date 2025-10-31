using DrWatson
@quickactivate "VAR_estimation"

using DataFrames, DataFramesMeta, CSV
using CairoMakie, UnicodePlots, Dates
using Statistics, StatsBase
using TerminalPager
using LinearAlgebra, Random
using PrettyTables

include(srcdir("functions.jl"))

# Defining some functions
d4_ln_fn = (x) -> x[5:end] - x[1:(end-4)]
function remove_outliers(data; lw = 0.10, up = 0.90)
    q1 = quantile(data, lw)
    q2 = quantile(data, up)
    lower_bound = q1
    upper_bound = q2
    return filter(x -> lower_bound < x < upper_bound, data)
end

# Define the length of the IRFs
l = 20

# Load the data
GT_log_data = CSV.read(
    datadir("data_log.csv"),
    DataFrame
)

GT_log_data.dates = Date("2001-03"):Month(3):Date("2024-12")

# data for headline inflation
data_headline = @chain GT_log_data begin
    @rsubset :dates >= Date("2003-3")
    @select :dates :ln_y :ln_cpi :i
end

# Transformations for d4_ln data
data_headline_d4_ln = copy(data_headline)
data_headline_d4_ln.d4_ln_y = [fill(NaN,4); d4_ln_fn(data_headline_d4_ln.ln_y)]
data_headline_d4_ln.d4_ln_cpi = [fill(NaN,4); d4_ln_fn(data_headline_d4_ln.ln_cpi)]
data_headline_d4_ln.r = [fill(NaN,4); data_headline_d4_ln.i[5:end] - data_headline_d4_ln.d4_ln_cpi[5:end]]

# convert to a matrix
data_headline_d4_ln_mat = (Matrix)(data_headline_d4_ln[data_headline_d4_ln.dates .> Date("2005-3"), [:d4_ln_y, :d4_ln_cpi, :r]])

############## estimation the SVAR(4) ##################
VAR_4 = VAR(data_headline_d4_ln_mat, 4)
IRFs_4 = IRF(VAR_4, l, true)
SR_point = sum(cumsum(IRFs_4[1,2,:]))/sum(IRFs_4[2,2,:])

# simulation
replications = 10000
SR_sim = Array{Float64}(undef, replications)
for j in 1:replications
    VAR_4 = VAR(data_headline_d4_ln_mat, 4)

    est_values = VAR_4["A"] * VAR_4["Z"]
    errors = VAR_4["Y"] - est_values
    errors = errors'

    errors_boots_array = Matrix{Float64}(undef, 75, 3)
    for var in 1:3
        rand_errors = Array{Float64}(undef, 75)
        for i in 1:75
            indx = rand(1:75)
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

SR_clean = remove_outliers(SR_sim, lw = 0.04, up = 0.96)
mean(SR_clean)
median(SR_clean)
std(SR_clean)
skewness(SR_clean)
kurtosis(SR_clean)

fig = Figure(size = (900, 600))

ax = Axis(
    fig[1,1]
)

hist!(ax, SR_clean)
vlines!(ax, SR_point, color = :red)
vlines!(ax, mean(SR_clean), color = :black)

fig

# bias adjusted
SR_point - (mean(SR_clean) - SR_point)

bound = mean(SR_clean) - SR_point

sum(-28.019 .< SR_clean .< 35.94)/length(SR_clean)



