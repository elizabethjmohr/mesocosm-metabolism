@model function fit_metabolism(::Type{T} = Float64; DOdata = missing, ODE, times, nDays, PARData, OSatData, GPP_daily_mu = 6.33, GPP_daily_sigma = 1.0, ER_daily_mu = -0.375, ER_daily_sigma = 1.0, k600_meanlog = 0.05, k600_sdlog = 0.2, Q = 1.794, dayStart = 4, K02_conv = 1.0, err_obs_iid_sigma_scale = 0.1, err_proc_iid_sigma_scale=0.1) where {T}
    PAR = Interpolations.interpolate((times,), PARData, Interpolations.Gridded(Interpolations.Linear()))
    OSat = Interpolations.interpolate((times,), OSatData, Interpolations.Gridded(Interpolations.Linear()))

    if DOdata === missing
        DOdata = Vector{T}(undef, length(times))
        DOdata[1] = ODE.u0[1]
    end

    DO_true = Vector{T}(undef, length(times))
    
    # Priors on daily GPP and ER, k600, observation error, and process error
    GPP_daily = Vector{T}(undef, nDays)
    ER_daily = Vector{T}(undef, nDays)

    for i = 1:nDays
        GPP_daily[i] ~ truncated(Normal(GPP_daily_mu, GPP_daily_sigma), 0, Inf)
        ER_daily[i] ~ truncated(Normal(ER_daily_mu, ER_daily_sigma), -Inf, 0)
    end

    k600 ~ LogNormal(k600_meanlog, k600_sdlog)
    kGas = K02_conv*k600

    err_obs_iid_sigma_scaled ~ truncated(Cauchy(0, 1), 0, Inf)
    err_proc_iid_sigma_scaled ~ truncated(Cauchy(0, 1), 0, Inf)

    err_obs_iid_sigma = err_obs_iid_sigma_scale * err_obs_iid_sigma_scaled
    err_proc_iid_sigma = err_proc_iid_sigma_scale * err_proc_iid_sigma_scaled

    # Specify model initial condition and parameters
    DO_true[1] ~ Normal(DOdata[1], err_obs_iid_sigma)

    p = Dict(
        "PAR" => PAR,
        "OSat" => OSat,
        "GPP_daily"=> GPP_daily,
        "ER_daily" => ER_daily, 
        "kGas" => kGas, 
        "Q" => Q, 
        "dayStart" => dayStart
        )
   
    problem = remake(
        ODE, 
        p = p, 
        u0 = [DO_true[1]],
        tspan = (minimum(times), maximum(times))
    )

    predicted = solve(
        problem, 
        OrdinaryDiffEq.AutoVern7(OrdinaryDiffEq.Rodas4()), 
        reltol = 1e-5, 
        saveat = times
    )

    for i = 2:length(predicted)
        DO_true[i] ~ Normal(predicted[i][1], err_proc_iid_sigma)
        DOdata[i] ~ Normal(DO_true[i], err_obs_iid_sigma)
    end
end