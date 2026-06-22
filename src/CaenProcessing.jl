
module CaenProcessing

using AstroParticleUnits
using StructArrayTables
using StatsBase

const ns = u"ns"
const mV = u"mV"
const C = u"C"
const pC = u"pC"

include("read_and_process.jl")
include("runs.jl")

export ns, mV, C, pC
export read_waveforms, process_data
export read_run, read_run_config, Run

export VOLTS_PER_ADC, TIME_PER_SAMPLE

end