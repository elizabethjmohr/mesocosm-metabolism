times = collect(4.0:1.0:75.0)
PARData = [0; 0; repeat([1/12], 12); repeat([repeat([0.0], 12); repeat([1/12], 12)], 2); repeat([0.0], 10)]
OSatData = repeat([7.0], length(times))
GPP_daily = [6.33, 6.33, 6.33]
ER_daily = [-0.375, -0.375, -0.375]
kGas = 1.05
Q = 1.794
dayStart = 4
u0 = [6.66411]

problem = initialize_process_model(
    times = times, 
    PARData = PARData, 
    OSatData = OSatData, 
    GPP_daily = GPP_daily, 
    ER_daily = ER_daily, 
    kGas = kGas, 
    Q = Q, 
    dayStart = dayStart, 
    u0 = u0
)

sol = solve(
    problem, 
    OrdinaryDiffEq.AutoVern7(OrdinaryDiffEq.Rodas4()), 
    reltol = 1e-5, 
    saveat = 0.15
)

synthData =  vec(sol + 0.01*randn(size(Array(sol))))

# Plot synthetic data with model solution
Plots.plot(sol, alpha = 0.3, legend = false); 
Plots.scatter!(sol.t, synthData)
