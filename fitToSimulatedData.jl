include("initializeDifferentialEquation.jl")
include("generateSyntheticData.jl")
include("bayes.jl")

# Specify model
model = fit_metab(synthDataDF, prob, nDays, nrow(synthDataDF), p["dayStart"], p["Q"], 
getPAR, getOSat, tspan[1], tspan[2], 
6.0, 1.0, -0.4, 1.0, 0.05, 1.0, 1.0, 0.1, 0.1)

# Sample posterior
# using Zygote, DiffEqSensitivity
# Turing.setadbackend(:zygote)
chain = mapreduce(c -> sample(model, NUTS(.65),10), chainscat, 1:3)