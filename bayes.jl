using Turing, Distributions
include("initializeDifferentialEquation.jl")

@model function fit_metab(data, prob, nDays, n, dayStart, Q, getPAR, getOSat, t0, tf, GPP_daily_mu, GPP_daily_sigma, ER_daily_mu, ER_daily_sigma, k600_meanlog, k600_sdlog, K02_conv, err_obs_iid_sigma_scale, err_proc_iid_sigma_scale)
    # Priors on daily GPP and ER, k600, observation error, and process error
    GPP_daily = zeros(nDays)
    for i = 1:nDays
        GPP_daily[i] ~ truncated(Normal(GPP_daily_mu, GPP_daily_sigma), 0, Inf)
    end
    
    ER_daily = zeros(nDays)
    for i = 1:nDays
        ER_daily[i] ~ truncated(Normal(ER_daily_mu, ER_daily_sigma), -Inf, 0)
    end

    k600 ~ LogNormal(k600_meanlog, k600_sdlog)
    kGas = K02_conv*k600

    err_obs_iid_sigma_scaled ~ truncated(Cauchy(0, 1), 0, Inf)
    #err_proc_iid_sigma_scaled ~ truncated(Cauchy(0, 1), 0, Inf)

    err_obs_iid_sigma = err_obs_iid_sigma_scale * err_obs_iid_sigma_scaled
    #err_proc_iid_sigma = err_proc_iid_sigma_scale * err_proc_iid_sigma_scaled

    # Specify model initial condition and parameters
    u0 = data[1,2]
    p = Dict("dailyGPP"=> GPP_daily,
            "ER" => ER_daily, 
            "kGas" => kGas, 
            "Q" => Q, 
            "getPAR" => getPAR, 
            "getOSat" => getOSat, 
            "dayStart" => dayStart, 
            "nDays"=> nDays) 
   
    problem = remake(prob, p = p, u0 = [u0], tspan = (t0, tf))
    predicted = solve(problem, AutoVern7(Rodas4()), reltol = 1e-5, saveat = data[:,1])
    
    #DO_mod = zeros(n)

    for i = 1:n
        #DO_mod[i] ~ Normal(predicted[i], err_proc_iid_sigma)
        data[i,2] ~ Normal(predicted[i][1], err_obs_iid_sigma[1])
    end

end