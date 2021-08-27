# Define differential equation describing change in DO concentration with respect to time
function mesocosm_metabolism!(du,u,p,t)
    du[1] = p["PAR"](t)*(p["GPP_daily"])[p["day"](t)] + 
            (p["ER_daily"])[p["day"](t)] + 
            (p["kGas"] + (p["Q"]/p["V"]))*(p["OSat"](t) - u[1]) 
end


