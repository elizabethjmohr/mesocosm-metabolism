@model function fit_metabolism(;ODE, times, nDays, DOdata = missing, PARData, OSatData, GPP_daily_mu = 6.33, GPP_daily_sigma = 1.0, ER_daily_mu = -0.375, ER_daily_sigma = 1.0, k600_meanlog = 1.05, k600_sdlog = 0.5, Q = 1.794, dayStart = 4, K02_conv = 1.0, err_obs_iid_sigma_scale = 0.1, err_proc_iid_sigma_scale=0.1) 
    PAR = interpolate((times,), PARData, Gridded(Linear()))
    OSat = interpolate((times,), OSatData, Gridded(Linear()))
    
    if DOdata === missing
        DOdata = Vector{Float64}(undef, length(times))
    end

    DO_true = Vector{Float64}(undef, length(times))

    # Priors on daily GPP and ER, k600, observation error, and process error
    GPP_daily = Vector{Float64}(undef, nDays)
    ER_daily = Vector{Float64}(undef, nDays)
    for i = 1:nDays
        GPP_daily[i] ~ truncated(Normal(GPP_daily_mu, GPP_daily_sigma), 0, Inf)
        ER_daily[i] ~ truncated(Normal(ER_daily_mu, ER_daily_sigma), -Inf, 0)
    end

    # Need to figure out how to truncate these multivariate distributions
    # GPP_daily ~ MvNormal([GPP_daily_mu for i in 1:nDays], GPP_daily_sigma.*Diagonal(ones(nDays)))
    # ER_daily ~ MvNormal([ER_daily_mu for i in 1:nDays], ER_daily_sigma.*Diagonal(ones(nDays)))

    #GPP_daily1 ~truncated(Normal(GPP_daily_mu, GPP_daily_sigma), 0, Inf)
    #GPP_daily2 ~truncated(Normal(GPP_daily_mu, GPP_daily_sigma), 0, Inf)
    #GPP_daily3 ~truncated(Normal(GPP_daily_mu, GPP_daily_sigma), 0, Inf)
    #ER_daily1 ~truncated(Normal(ER_daily_mu, ER_daily_sigma), -Inf, 0)
    #ER_daily2 ~truncated(Normal(ER_daily_mu, ER_daily_sigma), -Inf, 0)
    #ER_daily3 ~truncated(Normal(ER_daily_mu, ER_daily_sigma), -Inf, 0)

    k600 ~ LogNormal(k600_meanlog, k600_sdlog)
    kGas = K02_conv*k600

    err_obs_iid_sigma_scaled ~ truncated(Cauchy(0, 1), 0, Inf)
    err_proc_iid_sigma_scaled ~ truncated(Cauchy(0, 1), 0, Inf)

    err_obs_iid_sigma = err_obs_iid_sigma_scale * err_obs_iid_sigma_scaled
    err_proc_iid_sigma = err_proc_iid_sigma_scale * err_proc_iid_sigma_scaled

    # Specify model initial condition and parameters
    DO_true[1] ~ Normal(DOdata[1], err_obs_iid_sigma)

    # eltype(GPP_daily1).(u0)

    p = Dict(
        "PAR" => PAR,
        "OSat" => OSat,
        "GPP_daily"=> GPP_daily,
        "ER_daily" => ER_daily, 
        "kGas" => kGas, 
        "Q" => Q, 
        "dayStart" => dayStart
        )
   
    problem = remake(ODE, p = p, u0 = DO_true[1], tspan = (minimum(times), maximum(times)))
    predicted = solve(problem, AutoVern7(Rodas4()), reltol = 1e-5, saveat = DOdata[:,1])

    for i = 2:length(predicted)
        DO_true[i] ~ Normal(predicted[i][1], err_proc_iid_sigma)
        DOdata[i] ~ Normal(DO_true[i], err_obs_iid_sigma)
    end
end