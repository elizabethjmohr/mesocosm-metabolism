include("metabolismDifferentialEquation.jl")

# TODO: make this into a function in order to initialize DE based on priors, nDays, etc.


# Specify model parameters, output time range, and initial condition
nDays = 3
p = Dict("dailyGPP"=> [6.33, 6.33, 6.33],
         "ER" => [-0.375, -0.375, -0.375], 
         "kGas" => 1.05, 
         "Q" => 1.794, 
         "getPAR" => getPAR, 
         "getOSat" => getOSat, 
         "dayStart" => 4, 
         "nDays"=> nDays) 
tspan = (4.0, nDays * 24.0 + p["dayStart"] - 0.01) # Model time units = hours
u0 = [6.66411]

# Define and solve ODE problem
prob = ODEProblem(metab_SAM!, u0, tspan, p)
