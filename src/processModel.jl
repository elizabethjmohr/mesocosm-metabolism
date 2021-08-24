# Define differential equation describing change in DO concentration with respect to time
function mesocosm_metabolism!(du,u,p,t)
    PAR = p["PAR"](t)
    OSat = p["OSat"](t)
    day = p["day"](t)
    # day = convert(Int64, floor(t/24.0 - (p["dayStart"])/24.0)) + 1
    du[1] = PAR*(p["GPP_daily"])[day] + (p["ER_daily"])[day] + p["kGas"]*(OSat - u[1]) + p["Q"]*(OSat - u[1])/27.0
end


