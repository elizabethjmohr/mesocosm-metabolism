# initialize DE based on priors, nDays, etc.
function initialize_process_model(;times, PAR, OSat, day, GPP_daily, ER_daily, kGas = 1.05, Q = 1.794, V = 27.0, dayStart = 4, u0)

    if convert(Int64, floor(maximum(times)/24.0 - dayStart/24.0) + 1) > length(GPP_daily)
        error("The simulation timespan contains more days than are provided in GPP_daily")
    end
    
    p = Dict(
        "PAR" => PAR,
        "OSat" => OSat,
        "day" => day,
        "GPP_daily"=> GPP_daily,
        "ER_daily" => ER_daily, 
        "kGas" => kGas, 
        "Q" => Q,
        "V" => V)

    tspan = (minimum(times), maximum(times))

    # Define and solve ODE problem
    prob = ODEProblem(mesocosm_metabolism!, u0, tspan, p)
    return prob
end



