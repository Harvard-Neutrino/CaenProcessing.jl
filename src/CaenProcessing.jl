
module CaenProcessing

using AstroParticleUnits
using StructArrayTables
using StatsBase

const ns = u"ns"
const mV = u"mV"
const C = u"C"
const pC = u"pC"

include("waveforms.jl")
include("read_and_process.jl")
include("runs.jl")

export ns, mV, C, pC
export read_waveforms, process_data
export Waveforms, Waveform, sample_times, nsamples, volts
export read_run, read_run_config, Run, reprocess

export VOLTS_PER_ADC, TIME_PER_SAMPLE

end