@model function simulate_metabolism(::Type{T} = Float64; 
    DOdata = missing, 
    ODE, 
    times, 
    nDays, 
    PAR, 
    OSat, 
    day,
    GPP_daily_mu, 
    GPP_daily_sigma, 
    ER_daily_mu, 
    ER_daily_sigma, 
    k600_meanlog, 
    k600_sdlog, 
    Q, 
    dayStart, 
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

    # Specify model initial condition prior
    z0 ~ truncated(Normal(ODE.u0[1], err_obs_iid_sigma), 0, Inf)

    p = Dict(
        "PAR" => PAR,
        "OSat" => OSat,
        "day" => day,
        "GPP_daily"=> GPP_daily,
        "ER_daily" => ER_daily, 
        "kGas" => kGas, 
        "Q" => Q, 
        "dayStart" => dayStart
        )
   
    problem = remake(
        ODE, 
        p = p, 
        u0 = [z0],
        tspan = (minimum(times), maximum(times))
    )

    predicted = solve(
        problem, 
        OrdinaryDiffEq.AutoVern7(OrdinaryDiffEq.Rodas4P()), 
        reltol = 1e-5, 
        tstops = repeat([6,18], outer = nDays).+24*repeat([i for i in 0:(nDays-1)], inner = 2),
        #maxiters = 1e6,
        saveat = times
    )

    DO_true = Vector{T}(undef, (length(times) - 1))

    if DOdata === missing
        DOdata = Vector{T}(undef, length(times))
        DOdata[1] = ODE.u0[1]
    end

    for i = 2:length(predicted)
        DO_true[i-1] ~ Normal(predicted[i][1], err_proc_iid_sigma)
        DOdata[i] ~ Normal(DO_true[i-1], err_obs_iid_sigma)
    end

    return DOdata
end