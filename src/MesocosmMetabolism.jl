module MesocosmMetabolism

using DifferentialEquations, Interpolations, Turing, DataFrames, Plots, Distributions, LinearAlgebra

export mesocosm_metabolism!
export initialize_process_model
export fit_metabolism

include("processModel.jl")
include("initializeProcessModel.jl")
include("statisticalModel.jl")

end
