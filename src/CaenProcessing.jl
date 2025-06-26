
module CaenProcessing

using AstroParticleUnits
using TypedTables
using StructArrays
using StatsBase
# using Printf: @sprintf 

using FHist
using Distributions
using UnitfulDistributions
using LandauDistribution

using Minuit2

const ns = u"ns"
const mV = u"mV"
const C = u"C"

include("read_and_process.jl")
include("coincidence.jl")

include("langau.jl")
include("distributions.jl")
include("fitting.jl")

export ns, mV, C
export read_waveforms, process_data
export get_coincidence_mask

export fit_muon_distribution, find_muon_Vmin
export Gaussian, Langau

end