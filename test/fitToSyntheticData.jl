using DifferentialEquations, Turing, Distributions

# Specify model
model = fit_metabolism(
    DOdata = missing, 
    ODE = problem, 
    times = times, 
    nDays = 3, 
    PARData = PARData, 
    OSatData = OSatData, 
    GPP_daily_mu = 6.33, 
    GPP_daily_sigma = 1.0, 
    ER_daily_mu = -0.375, 
    ER_daily_sigma = 1.0, 
    k600_meanlog = 0.05, 
    k600_sdlog = 0.2, 
    Q = 1.794, 
    dayStart = 4, 
    K02_conv = 1.0, 
    err_obs_iid_sigma_scale = 0.1, 
    err_proc_iid_sigma_scale= 0.1
)

# Sample posterior
chains = sample(model, NUTS(0.65), 1)
synthData = DataFrames.stack(select(DataFrame(chains), r"DOdata"))[!,"value"]
