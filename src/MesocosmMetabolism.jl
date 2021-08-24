module MesocosmMetabolism

using DifferentialEquations
using Interpolations
using Turing
using DataFrames
using Plots
using Distributions
using LinearAlgebra
using OrdinaryDiffEq
using FillArrays

import Distributions.truncated
import Distributions.Normal
import Distributions.Cauchy
import Distributions.LogNormal
import OrdinaryDiffEq.AutoVern7
import OrdinaryDiffEq.Rodas4

include("processModel.jl")
include("initializeProcessModel.jl")
include("statisticalModel.jl")
include("generateSyntheticData.jl")
include("simulateMetabolism.jl")

export mesocosm_metabolism!
export initialize_process_model
export fit_metabolism
export generate_synthetic_data
export simulate_metabolism

end
