function generate_synthetic_data(
    ;times = collect(4.0:1.0:75.0), 
    nDays, 
    PAR, 
    OSat, 
    day = day,
    GPP_daily_mu, 
    GPP_daily_sigma, 
    ER_daily_mu, 
    ER_daily_sigma, 
    k600_meanlog, 
    k600_sdlog, 
    Q , 
    dayStart, 
    K02_conv , 
    err_obs_iid_sigma_scale, 
    err_proc_iid_sigma_scale,
    u0)

    # TODO: Put in a check to make sure that times and nDays are compatible. Throw error if not.

    problem = initialize_process_model(
        times = times, 
        PAR = PAR, 
        OSat = OSat, 
        day = day,
        GPP_daily = fill(GPP_daily_mu, nDays), 
        ER_daily = fill(ER_daily_mu, nDays), 
        kGas = exp(k600_meanlog), 
        Q = Q, 
        dayStart = dayStart, 
        u0 = u0
    )
    
    model = simulate_metabolism(
        DOdata = missing, 
        ODE = problem, 
        times = times, 
        nDays = nDays, 
        PAR = PAR, 
        OSat = OSat, 
        day = day,
        GPP_daily_mu = GPP_daily_mu, 
        GPP_daily_sigma = GPP_daily_sigma, 
        ER_daily_mu = ER_daily_mu, 
        ER_daily_sigma = ER_daily_sigma, 
        k600_meanlog = k600_meanlog, 
        k600_sdlog = k600_sdlog, 
        Q = Q, 
        dayStart = dayStart, 
        K02_conv = K02_conv, 
        err_obs_iid_sigma_scale = err_obs_iid_sigma_scale, 
        err_proc_iid_sigma_scale= err_proc_iid_sigma_scale)
    
    samp = model()
    return(samp)
end