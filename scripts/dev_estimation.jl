using DrWatson
@quickactivate "VAR_estimation"

using DataFrames, DataFramesMeta, CSV
using CairoMakie
using Statistics
using TerminalPager

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
    Z_array[i,:] = hcat(gtdata[i+1, :]', gtdata[i, :]')
end
Z = Z_array'

A = (Z' \ Y')' 

Sigma_u = (Y - A*Z) * (Y - A*Z)' / (92 - 3*2 - 1)

############### Estimate the VAR model by LS estimation function ###############
VAR_estimation = VAR(gtdata, 4)

# get the estimated parameters
A = get_params(VAR_estimation[1], 3, 4)
A1 = A[1]
A2 = A[2]
A3 = A[3]
A4 = A[4]

############### Long term restrictitons ###############


