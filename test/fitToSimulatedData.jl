# Specify model
model = fit_metab(synthDataDF, prob, nDays, nrow(synthDataDF), p["dayStart"], p["Q"], 
PAR, OSat, tspan[1], tspan[2], 
6.33, 0.01, -0.375, 0.01, 0.05, 0.01, 1.0, 0.1, 0.1)

# Sample posterior
chains = sample(model, NUTS(), 10)