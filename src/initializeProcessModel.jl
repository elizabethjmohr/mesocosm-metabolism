# initialize DE based on priors, nDays, etc.
function initialize_process_model(;times, PARData, OSatData, GPP_daily, ER_daily, kGas = 1.05, Q = 1.794, dayStart = 4, u0)
    
    PAR = interpolate((times,), PARData, Gridded(Linear()))
    OSat = interpolate((times,), OSatData, Gridded(Linear()))

    if convert(Int64, floor(maximum(times)/24.0 - dayStart/24.0) + 1) > length(GPP_daily)
        error("The simulation timespan contains more days than are provided in GPP_daily")
    end
    
    p = Dict(
        "PAR" => PAR,
        "OSat" => OSat,
        "GPP_daily"=> GPP_daily,
        "ER_daily" => ER_daily, 
        "kGas" => kGas, 
        "Q" => Q, 
        "dayStart" => dayStart
        )

    tspan = (minimum(times), maximum(times))

    # Define and solve ODE problem
    prob = ODEProblem(mesocosm_metabolism!, u0, tspan, p)
end



