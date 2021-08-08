# MesocosmMetabolism.jl

This Julia package provides tools for: 
1. Simulating dissolved oxygen concentration time-series in stream-analogue mesocosms (SAMS) based on daily gross primary production, daily ecosystem respiration, gas exchange, photosynthetically active radiation, saturated oxygen concentrations, and the volume of water flowing into and out from the SAM. 
2. Estimating daily gross primary production, daily ecosystem respiration, and gas exchange given times series of dissolved oxygen concentrations, photosynthetically active radiation, saturated oxygen concentrations, and the volume of water flowing into and out from the SAM.  

![Screen Shot 2021-08-08 at 3 31 41 PM](https://user-images.githubusercontent.com/34284337/128646339-4786c4df-8605-41cd-8983-3f590e45d42b.png)

Cutaway of a stream-analogue mesocosm, or SAM

## Functions

**mesocosm_metabolism!**: defines the ODE describing changes in oxygen concentration with respect to time

**initialize_process_model**: instantiates an "ODEProblem" given a set of parameters, an initial condition, and a simulation timespan

**fit_metabolism**: defines a statistical model that can sampled using any of the MCMC algorithms available in Turing. If DO data is supplied as an argument, the model resulting from this function can be sampled for the purposed of obtaining estimates of GPP, ER, and gas exchange. If DO data is not supplied, the model resulting form this function will generate "synthetic data". 




