using DrWatson
@quickactivate "VAR_estimation"

using DataFrames, DataFramesMeta, CSV
using CairoMakie, UnicodePlots
using Statistics
using TerminalPager
using LinearAlgebra
using Infiltrator
using Revise

include(srcdir("functions.jl"))

# Load the data
gtdata = CSV.read(
    datadir("data.csv"),
    DataFrame
)

# Transform in a Matrix
gtdata = (Matrix)(gtdata[:, [:dla_gdp, :dla_cpi, :r]])
# The data comes from 2002Q1 to 2024Q4

############### Estimate the VAR model by LS estimation ###############
# Following illian and Lütkepohl the LS estimator is:
# A = [v, A1, A2,..., Ap] = YZ'(ZZ')^-1
# Where Y = [y1, y2, ..., yT] is KxT and Z = [Z0,...,ZT-1] is (Kp+1)xT
# Zt-1 = (1, y't-1, ... , y't-p)

# out propose is estimate a VAR(2) to the data of Guatemala

# First create the Y KxT matrix
# yt in the Y matrix are the variables in the model for its aproppiate t periods
gtdata[3, :] # is the first element of Y. Follows this process:

Y = mapreduce(hcat, 3:92) do x
    gtdata[x, :]
end

# Second create the Z KpxT matrix
hcat(gtdata[2, :]', gtdata[1, :]') # Z0

Z_array = Array{Float64}(undef, 90, 6)
for i in 1:90
    Z_array[i, :] = hcat(gtdata[i+1, :]', gtdata[i, :]')
end
Z = Z_array'

A = (Z' \ Y')'

Sigma_u = (Y - A * Z) * (Y - A * Z)' / (92 - 3 * 2 - 1)

############### Estimate the VAR model by LS estimation function ###############
VAR_estimation = VAR(gtdata, 4)

# get the estimated parameters
A = get_params(VAR_estimation)
A1 = A[1]
A2 = A[2]
A3 = A[3]
A4 = A[4]

############### estimations ###############
estimations = VAR_estimation["A"] * VAR_estimation["Z"]
errors = VAR_estimation["Y"] - estimations

fig = Figure(size=(900, 600))

ax = Axis(fig[1, 1])

lines!(ax, estimations[1, :], label="estimation")
lines!(ax, estimations[1, :] + errors[1, :], label="estimation + errors", color=(:green), linestyle=:dash, linewidth=2)
lines!(ax, VAR_estimation["Y"][1, :], label="observed")

axislegend(ax)

fig

############### Longterm restrictions ###############
# Applyin the Blanchad and Qua method
# we need the Theta matrix: chol(A(1)^1ΣuA(1)^1')
# where A(1) = (I_k - A1 - A2 - ... - Ap)
A_1 = Matrix(I, 3, 3) - A1 - A2 - A3 - A4
A_1_inv = inv(A_1)
Θ_1 = cholesky(Hermitian(A_1_inv * VAR_estimation["Σu"] * A_1_inv')).L

B_0_inv = A_1 * Θ_1
B_0 = inv(B_0_inv)

# structural parameters
B1 = B_0 * A1
B2 = B_0 * A2
B3 = B_0 * A3
B4 = B_0 * A4

B = hcat(B1, B2, B3, B4)

# comparing with the BQ_VAR function
BQ_VAR(VAR_estimation)["B_0_inv"] .- B_0_inv # produce the same results

############### Impulse Response Functions ###############
# J matrix J = [I_k, 0kxk(p-1)]
J = hcat(Matrix(I, 3, 3), zeros(3, 9))
# Matrix
A = vcat(
    hcat(A1, A2, A3, A4),
    hcat(Matrix(I, 9, 9), zeros(9, 3))
)

periods = 20
IRF_struct = Array{Float64}(undef, 3, 3, periods)

for i in 1:periods
    IRF_struct[:, :, i] = (J * A^i * J') * B_0_inv
end

s = 20

fig = Figure(size=(900, 600))

ax = Axis(
    fig[1, 1],
    title="GDP Structural Impulse Function\nTo inflation shock",
    xgridvisible=false,
    ygridvisible=false
)

lines!(
    IRF_struct[1, 2, 1:s]
)

hlines!(
    ax,
    0
)

ax = Axis(
    fig[2, 1],
    title="Inflatión Structural Impulse Function\nTo inflation shock",
    xgridvisible=false,
    ygridvisible=false
)

lines!(
    IRF_struct[2, 2, 1:s]
)

hlines!(
    ax,
    0
)

fig

#Sacrifice ratio
sum(IRF_struct[1, 2, 1:s]) / sum(IRF_struct[2, 2, 1:s])

# comparing with the IRF function
IRF_function = IRF(VAR_estimation, periods, true)
# Sacrifice ratio
sum(IRF_function[1, 2, :]) / sum(IRF_function[2, 2, :])


