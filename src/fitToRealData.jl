using MesocosmMetabolism, CSV, Turing, DataFrames, MCMCChains, StatsPlots, Interpolations

data = CSV.read("/Users/elizabethmohr/Documents/MSU/RProjects/stream-analogue-mesocosms/data/test1.csv", DataFrame)

Q = 1.794
dayStart = 4
times = data.modelTime # hours
PAR = Interpolations.ConstantInterpolation((times,), data.PAR)
OSat = Interpolations.interpolate((times,), data."DO.sat", Interpolations.Gridded(Interpolations.Linear()))
nDays = 1
days = [convert(Int64, floor(times[i]/24.0 - (dayStart/24.0)) + 1) for i in 1:length(times)]
day = Interpolations.ConstantInterpolation((times,), days)


# Note: these priors are taken from the default values in streamMetabolizer,
# converted to units of hours and converted to a per volume basis
# by dividing by mean depth in SAMS of 0.15 m.
GPP_daily_mu = 1.72 # g/(L*hr) - only applies to 12 hours when lights are on
GPP_daily_sigma = 3.33
ER_daily_mu = -2.0 # g/(L*hr)
ER_daily_sigma = 2.0
k600_meanlog = -0.69 # hr^-1
k600_sdlog = 0.53
K02_conv = 1.0
err_obs_iid_sigma_scale = 0.03
err_proc_iid_sigma_scale = 5.0
u0 = [data.DO_mgL[1]]

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

model = fit_metabolism(
    DOdata = data.DO_mgL, 
    ODE = problem, 
    times = times, 
    nDays = 1, 
    PAR = PAR, 
    OSat = OSat, 
    day = day,
    GPP_daily_mu = GPP_daily_mu, 
    GPP_daily_sigma = GPP_daily_sigma, 
    ER_daily_mu = ER_daily_mu, 
    ER_daily_sigma = ER_daily_sigma, 
    k600_meanlog = k600_meanlog, 
    k600_sdlog = k600_sdlog, 
    Q =  Q, 
    dayStart = dayStart, 
    K02_conv = K02_conv, 
    err_obs_iid_sigma_scale = err_obs_iid_sigma_scale, 
    err_proc_iid_sigma_scale= err_proc_iid_sigma_scale
)

setadbackend(:forwarddiff)
chains = sample(model, NUTS(0.80), MCMCThreads(), 10, 3)

gelmandiag(chains)

df = DataFrame(chains)

# GPP trace plot
@df df plot(
    :iteration,
    cols(columnindex(df,"GPP_daily[1]")), 
    group = (:chain),
    lw = 3,
    xlabel = "Iteration", 
    ylabel = "GPP", 
    legendtitle = "Chain",
    xguidefontsize=15, 
    yguidefontsize = 15)

# GPP posterior and prior  
@df df density(
    cols(columnindex(df,"GPP_daily[1]")), 
    group = (:chain),
    lw = 3,
    xlabel = "GPP", 
    ylabel = "Probability Density",
    legendtitle = "Chain",
    xguidefontsize=15, 
    yguidefontsize = 15)

plot!(
    truncated(Normal(GPP_daily_mu, GPP_daily_sigma), 0, Inf),
    lw = 3)

# ER trace plot
@df df plot(
    :iteration,
    cols(columnindex(df,"ER_daily[1]")), 
    group = (:chain),
    lw = 3,
    xlabel = "Iteration", 
    ylabel = "ER", 
    legendtitle = "Chain",
    xguidefontsize=15, 
    yguidefontsize = 15)

# ER posterior and prior  
@df df density(
    cols(columnindex(df,"ER_daily[1]")), 
    group = (:chain),
    lw = 3,
    xlabel = "ER", 
    ylabel = "Probability Density",
    legendtitle = "Chain",
    xguidefontsize= 15, 
    yguidefontsize = 15)

plot!(
    truncated(Normal(ER_daily_mu, ER_daily_sigma), 0, Inf),
    lw = 3)

# k600 trace plot
@df df plot(
    :iteration,
    cols(columnindex(df,"k600")), 
    group = (:chain),
    lw = 3,
    xlabel = "Iteration", 
    ylabel = "k600", 
    legendtitle = "Chain",
    xguidefontsize=15, 
    yguidefontsize = 15)

# k600 posterior and prior  
@df df density(
    cols(columnindex(df,"k600")), 
    group = (:chain),
    lw = 3,
    xlabel = "k600", 
    ylabel = "Probability Density",
    legendtitle = "Chain",
    xguidefontsize= 15, 
    yguidefontsize = 15)

plot!(
    LogNormal(k600_meanlog, k600_sdlog),
    lw = 3)

# Plot data with modeled DO from multiple samples
using DifferentialEquations
DOPlot = StatsPlots.scatter(
    times, 
    data.DO_mgL, 
    ylabel = "DO (mg/L)",
    xguidefontsize= 15, 
    yguidefontsize = 15);
n = 100
rows = rand(100:nrow(df), n)
for k in 1:n
    GPP = Array(select(df, "GPP_daily[1]")[rows[k],:])
    ER = Array(select(df, "ER_daily[1]")[rows[k], :])
    k600 = Array(select(df, "k600")[rows[k],:])[1]
    problem = initialize_process_model(
        times = times, 
        PARData = PARData, 
        OSatData = OSatData, 
        GPP_daily = GPP, 
        ER_daily = ER, 
        kGas = k600, 
        Q = Q, 
        dayStart = dayStart, 
        u0 = u0)
    sol = solve(
        problem, 
        OrdinaryDiffEq.AutoVern7(OrdinaryDiffEq.Rodas4P2()), 
        reltol = 1e-5, 
        maxiters = 1e6,
        tstops = repeat([7,19], outer = nDays).+24*repeat([i for i in 0:(nDays-1)], inner = 2),
        saveat = times
    )
    plot!(sol, alpha=0.2, color = "#BBBBBB", legend = :false)
end
xlabel!("Hour of day")

display(DOPlot)

# Plot prior distributions
using Plots, Distributions, StatsPlots
plot(Normal(GPP_daily_mu, GPP_daily_sigma))
plot(Normal(ER_daily_mu, ER_daily_sigma))
plot(LogNormal(k600_meanlog, k600_sdlog))
plot(truncated(Cauchy(0, err_obs_iid_sigma_scale), 0, Inf))
plot(truncated(Cauchy(0, err_proc_iid_sigma_scale), 0, Inf))