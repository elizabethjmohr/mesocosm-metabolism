using Interpolations, DifferentialEquations, Turing, FillArrays, ReverseDiff

# Define ODE
function mesocosm_metabolism!(du,u,p,t)
    du[1] = p["PAR"](t)*(p["GPP_daily"])[p["day"](t)] + 
            (p["ER_daily"])[p["day"](t)] + 
            (p["kGas"] + (p["Q"]/27.0))*(p["OSat"](t) - u[1]) 
end

# Generate some fake data
times = collect(4.0:0.1:75.0)
nDays = 3
dayStart = 4
days = [convert(Int64, floor(times[i]/24.0 - (dayStart/24.0)) + 1) for i in 1:length(times)]
PARData= [if (7<= (times[i] - 24*(days[i]-1)) < 19) 1 else 0 end for i in 1:length(times)]
OSatData = repeat([8.2], length(times))

p = Dict(
    "PAR" => Interpolations.ConstantInterpolation((times,), PARData),
    "OSat" => Interpolations.interpolate((times,), OSatData, Interpolations.Gridded(Interpolations.Linear())),
    "day" => Interpolations.ConstantInterpolation((times,), days),
    "GPP_daily"=> fill(1.8, nDays),
    "ER_daily" => fill(-1.8, nDays), 
    "kGas" => 0.6, 
    "Q" => 1.794)

tspan = (minimum(times), maximum(times))
u0 = [5.5]

problem = ODEProblem(mesocosm_metabolism!, u0, tspan, p)
sol = solve(problem, alg = AutoVern7(Rodas4()),reltol = 1e-5, maxiters = 1e6, tstops = [7,19], saveat = times)
fakeData = vec(sol) + 0.2*rand(length(times))

# Specify statistical model
@model function f(::Type{T} = Float64; DOdata, ODE) where {T}
    # Priors on daily GPP and ER, k600, observation error, and process error
    GPP_daily ~ MvNormal(FillArrays.Fill(1.80, $nDays), 2.0) 
    ER_daily ~ MvNormal(FillArrays.Fill(-1.80, $nDays), 2.0)  
    k600 ~ LogNormal(-0.69, 0.53)
    err_obs_iid_sigma ~ truncated(Cauchy(0, 0.03), 0, Inf)
    err_proc_iid_sigma ~ truncated(Cauchy(0, 5.0), 0, Inf)
    
    # ODE parameters
    p = Dict(
        "PAR" => ODE.p["PAR"],
        "OSat" => ODE.p["OSat"],
        "day" => ODE.p["day"],
        "GPP_daily"=> GPP_daily,
        "ER_daily" => ER_daily, 
        "kGas" => k600, 
        "Q" => ODE.p["Q"])
        
    DO_true = Vector{T}(undef, (length(times)))
    DO_true[1] ~ Normal(DOdata[1], err_obs_iid_sigma)

    for i = 2:length(times)
        problem = remake(ODE, p = p, u0 = [DO_true[i-1]], tspan = (times[i-1], times[i]))
        predicted = solve(problem, alg = Euler(), dt = 0.1)
        DO_true[i] ~ Normal(predicted(times[i])[1], err_proc_iid_sigma)
        DOdata[i] ~ Normal(DO_true[i], err_obs_iid_sigma)
    end
end

model = f(DOdata = fakeData, ODE = problem)
setadbackend(:reversediff)
samples = sample(model, NUTS(0.65), 10)

@model function g(::Type{T} = Float64; DOdata, ODE, times, PAR, OSat, day, Q, V) where {T}
    
   # Priors on daily GPP and ER, k600, observation error, and process error
   GPP_daily ~ MvNormal(FillArrays.Fill(1.80, $nDays), 2.0) 
   ER_daily ~ MvNormal(FillArrays.Fill(-1.80, $nDays), 2.0)  
   kGas ~ LogNormal(-0.69, 0.53)
   err_obs_iid_sigma ~ truncated(Cauchy(0, 0.03), 0, Inf)
   err_proc_iid_sigma ~ truncated(Cauchy(0, 5.0), 0, Inf)
   
   DO_mod = Vector{T}(undef, (length(times)))
   DO_true = Vector{T}(undef, (length(times)))
   DO_true[1] ~ Normal(DOdata[1], err_obs_iid_sigma)

   for i = 2:length(times)
       DO_mod[i] = DO_true[i-1] + (times[i] - times[i-1])*
        (PAR(times[i-1])*GPP_daily[day(times[i-1])] + 
        ER_daily[day(times[i-1])] + 
        (kGas + Q/V)*(OSat(times[i-1]) - DO_true[i-1]))
       DO_true[i] ~ Normal(DO_mod[i], err_proc_iid_sigma)
       DOdata[i] ~ Normal(DO_true[i], err_obs_iid_sigma)
   end

end

V = 27.0
model = g(DOdata = fakeData, ODE = problem, times = times, PAR = PAR, OSat = OSat, day = day, Q = Q, V = V)
samples = sample(model, NUTS(0.65), 10)