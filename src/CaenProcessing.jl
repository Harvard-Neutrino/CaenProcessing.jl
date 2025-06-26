
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

export ns, mV, C
export read_waveforms, process_data

end