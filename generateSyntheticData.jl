include("metabolismDifferentialEquation.jl")
using Plots, DataFrames

# Define function that returns PAR for any time t
function getPAR(t)
    if 6 <=rem(t, 24)< 18
        return 1/12
    else
        return 0.0
    end
end

# Define function that returns saturated O2 concentration for any time t
function getOSat(t)
    return(7.0)
end

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

# Generate synthetic data by solving the model and adding noise
sol1 = solve(prob, AutoVern7(Rodas4()), reltol = 1e-5, saveat = 0.15)
# Note - wierd artifact appears if relative tolerance is set to default
synthData =  sol1 + 0.01*randn(size(Array(sol1)))
plot(sol1, alpha = 0.3, legend = false); scatter!(sol1.t, synthData')
synthDataDF = DataFrame(t = sol1.t, DO = vec(synthData))

# Plot solution
plot(sol1)