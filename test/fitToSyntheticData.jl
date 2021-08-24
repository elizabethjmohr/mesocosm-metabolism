using MesocosmMetabolism, CSV, Turing, DataFrames, MCMCChains, StatsPlots, Interpolations
times = collect(4.0:0.1:75.0)
nDays = 3
dayStart = 4
days = [convert(Int64, floor(times[i]/24.0 - (dayStart/24.0)) + 1) for i in 1:length(times)]
PARData= [if (7<= (times[i] - 24*(days[i]-1)) < 19) 1 else 0 end for i in 1:length(times)]
OSatData = repeat([8.2], length(times))
day = Interpolations.ConstantInterpolation((times,), days)
PAR = Interpolations.ConstantInterpolation((times,), PARData)
OSat = Interpolations.interpolate((times,), OSatData, Interpolations.Gridded(Interpolations.Linear()))
GPP_daily_mu = 1.72 # g/(L*hr) - only applies to 12 hours when lights are on
GPP_daily_sigma = 0.01
ER_daily_mu = -2.0 # g/(L*hr)
ER_daily_sigma = 0.01
k600_meanlog = -0.2 # hr^-1
k600_sdlog = 0.01
Q = 1.794
dayStart = 4
K02_conv = 1.0
err_obs_iid_sigma_scale = 0.01
err_proc_iid_sigma_scale = 0.01
u0 = [6.0]

using Random
Random.seed!(1234)
synthData = generate_synthetic_data(
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
    err_proc_iid_sigma_scale = err_proc_iid_sigma_scale,
    u0 = u0)

plot(times, synthData, seriestype = :scatter, legend = false)

# Change priors
GPP_daily_mu = 1.72 # g/(L*hr) - only applies to 12 hours when lights are on
GPP_daily_sigma = 3.33
ER_daily_mu = -2.0 # g/(L*hr)
ER_daily_sigma = 2.0
k600_meanlog = -0.69 
k600_sdlog = 0.53
err_obs_iid_sigma_scale = 0.01
err_proc_iid_sigma_scale = 0.01

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
    u0 = u0)

# Specify model
model = fit_metabolism(
    DOdata = synthData, 
    ODE = problem, 
    times = times, 
    nDays = nDays, 
    PAR = PAR, 
    OSat = OSat, 
    day = day,
    GPP_daily_mu = GPP_daily_mu, 
    GPP_daily_sigma = GPP_daily_sigma, 
    ER_daily_mu = ER_daily_mu, 
    ER_daily_sigma = ER_daily_sigma , 
    k600_meanlog = k600_meanlog, 
    k600_sdlog = k600_sdlog, 
    Q =  Q, 
    dayStart = dayStart, 
    K02_conv = K02_conv, 
    err_obs_iid_sigma_scale = err_obs_iid_sigma_scale, 
    err_proc_iid_sigma_scale = err_proc_iid_sigma_scale)

# Sample posterior
chains = sample(model, NUTS(0.80), MCMCThreads(), 10, 3)

# Evaluate chain mixing
gelmandiag(chains)

# Trace plot for GPP on day 2
df = DataFrame(chains)
@df df plot(
    :iteration,
    cols(4), 
    group = (:chain),
    lw = 3,
    xguidefontsize=15, 
    xlabel = "Iteration", 
    ylabel = "GPP: Day 2", 
    yguidefontsize = 15,
    legendtitle = "Chain", 
    background_color = :black,
    color_palette = palette([colorant"#3D8D59", colorant"#77BDC6",colorant"#18537E"]),
    legend = :false)

# Plot posterior and prior for GPP on day 2
@df df density(
    cols(4), 
    group = (:chain),
    lw = 3,
    xguidefontsize=15, 
    xlabel = "GPP: Day 2", 
    ylabel = "Probability Density",
    yguidefontsize = 15,
    legendtitle = "Chain", 
    background_color = :black,
    color_palette = palette([colorant"#3D8D59", colorant"#77BDC6",colorant"#18537E"]),
    legend = :false)


