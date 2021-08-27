@model function simulate_data(::Type{T} = Float64; 
    DOdata = missing, 
    ODE, 
    times, 
    nDays, 
    GPP_daily_mu, 
    GPP_daily_sigma, 
    ER_daily_mu, 
    ER_daily_sigma, 
    k600_meanlog, 
    k600_sdlog, 
    K02_conv, 
    err_obs_iid_sigma_scale, 
    err_proc_iid_sigma_scale
    ) where {T}

    # Priors on daily GPP and ER, k600, observation error, and process error
    GPP_daily ~ MvNormal(FillArrays.Fill(GPP_daily_mu, nDays), GPP_daily_sigma) 
    ER_daily ~ MvNormal(FillArrays.Fill(ER_daily_mu, nDays), ER_daily_sigma)  
    k600 ~ LogNormal(k600_meanlog, k600_sdlog)
    kGas = K02_conv*k600
    err_obs_iid_sigma ~ truncated(Cauchy(0, err_obs_iid_sigma_scale), 0, Inf)
    err_proc_iid_sigma ~ truncated(Cauchy(0, err_proc_iid_sigma_scale), 0, Inf)

    p = Dict(
        "PAR" => ODE.p["PAR"],
        "OSat" => ODE.p["OSat"],
        "day" => ODE.p["day"],
        "GPP_daily"=> GPP_daily,
        "ER_daily" => ER_daily, 
        "kGas" => kGas, 
        "Q" => ODE.p["Q"],
        "V" => ODE.p["V"])

    if DOdata === missing
        DOdata = Vector{T}(undef, length(times))
        DOdata[1] = ODE.u0[1]
    end

    DO_true = Vector{T}(undef, (length(times)))
    DO_true[1] ~ Normal(DOdata[1], err_obs_iid_sigma)

    for i = 2:length(times)
        problem = remake(ODE, p = p, u0 = [DO_true[i-1]], tspan = (times[i-1], times[i]))
        predicted = solve(problem, alg = Trapezoid())
        DO_true[i] ~ Normal(predicted(times[i])[1], err_proc_iid_sigma)
        DOdata[i] ~ Normal(DO_true[i], err_obs_iid_sigma)
    end

    return DOdata
end