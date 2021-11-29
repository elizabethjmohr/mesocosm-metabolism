using MesocosmMetabolism, CSV, Turing, DataFrames, MCMCChains, StatsPlots, Interpolations, ReverseDiff

data = CSV.read("/Users/elizabethmohr/Documents/MSU/RProjects/stream-analogue-mesocosms/data/DO.csv", DataFrame)
Q = 1.794 #L/hr
V = 27.0 # L
dayStart = 4

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

setadbackend(:reversediff)

results = Dict()
for i in unique(data[:,:loggerID])[1:4]
    subset = filter(r -> r[:loggerID] == i, data)
    times = subset.modelTime
    PAR = Interpolations.ConstantInterpolation((times,), subset.PAR)
    OSat = Interpolations.interpolate((times,), subset."DO.sat", Interpolations.Gridded(Interpolations.Linear()))
    days = [convert(Int64, floor(times[i]/24.0 - (dayStart/24.0)) + 1) for i in 1:length(times)]
    day = Interpolations.ConstantInterpolation((times,), days)
    u0 = [subset.DO_mgL[1]]
    nDays = convert(Int64,ceil((maximum(times) - minimum(times))/24))
    problem = initialize_process_model(
        times = times, 
        PAR = PAR, 
        OSat = OSat, 
        day = day,
        GPP_daily = fill(GPP_daily_mu, nDays), 
        ER_daily = fill(ER_daily_mu, nDays), 
        kGas = exp(k600_meanlog), 
        Q = Q, 
        V = V,
        dayStart = dayStart, 
        u0 = u0)
    model = statisticalModel(
        nDays, 
        GPP_daily_mu, 
        GPP_daily_sigma, 
        ER_daily_mu, 
        ER_daily_sigma, 
        k600_meanlog, 
        k600_sdlog, 
        K02_conv, 
        err_obs_iid_sigma_scale, 
        err_proc_iid_sigma_scale)

    sampleMe = model(DOdata = subset.DO_mgL, ODE = problem,times = times)
    chains = sample(sampleMe, NUTS(), MCMCThreads(), 1000, 3) 
    results[i] = chains
end

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

plot!(Normal(GPP_daily_mu, GPP_daily_sigma),lw = 3)

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

plot!(Normal(ER_daily_mu, ER_daily_sigma),lw = 3)

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

plot!(LogNormal(k600_meanlog, k600_sdlog),lw = 3)

# Plot data with modeled DO from multiple samples
using DifferentialEquations
logger = 415148
MesocosmData = filter(r -> r[:loggerID] == logger, data)
df = DataFrame(results[logger])
times = MesocosmData.modelTime
PAR = Interpolations.ConstantInterpolation((times,), MesocosmData.PAR)
OSat = Interpolations.interpolate((times,), MesocosmData."DO.sat", Interpolations.Gridded(Interpolations.Linear()))
days = [convert(Int64, floor(times[i]/24.0 - (dayStart/24.0)) + 1) for i in 1:length(times)]
day = Interpolations.ConstantInterpolation((times,), days)
u0 = [MesocosmData.DO_mgL[1]]
nDays = convert(Int64,ceil((maximum(times) - minimum(times))/24))

DOPlot = StatsPlots.scatter(
    MesocosmData.modelTime, 
    MesocosmData.DO_mgL, 
    ylabel = "DO (mg/L)",
    xguidefontsize= 15, 
    yguidefontsize = 15);

n = 100
rows = rand(1:nrow(df), n)

for k in 1:n
    GPP = Array(df[k, r"^GPP_daily"])
    ER = Array(df[k, r"^ER_daily"])
    k600 = df[k, "k600"]
    problem = initialize_process_model(
        times = MesocosmData.modelTime, 
        PAR = PAR, 
        OSat = OSat,
        day = day,  
        GPP_daily = GPP, 
        ER_daily = ER, 
        kGas = k600, 
        Q = Q, 
        V = V,
        dayStart = dayStart, 
        u0 = u0)
    sol = solve(problem, alg = Euler(), dt = 0.01)
    plot!(sol, alpha=0.3, color = "#BBBBBB", legend = :false)
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