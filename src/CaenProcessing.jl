
module CaenProcessing

using AstroParticleUnits
using TypedTables
using StructArrays
using StatsBase

const ns = u"ns"
const mV = u"mV"
const C = u"C"

include("read_and_process.jl")

export ns, mV, C
export read_waveforms, process_data

export volts_per_ADC, time_per_sample

end