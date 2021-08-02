using DifferentialEquations

# Define differential equation describing change in DO concentration with respect to time
function metab_SAM!(du,u,p,t)
    PAR = p["getPAR"](t)
    OSat = p["getOSat"](t)
    day = convert(Int64, floor(t/24.0 - (p["dayStart"])/24.0)) +1
    du[1] = PAR*(p["dailyGPP"])[day] + (p["ER"])[day] + p["kGas"]*(OSat - u[1]) + p["Q"]*(OSat - u[1])/27.0
end


