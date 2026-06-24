
# Containers for the raw waveform readout.
#
# `read_waveforms` returns a `Waveforms` (one channel's ADC matrix plus per-event
# metadata); indexing it with `w[i]` yields a single-event `Waveform` that carries
# its own 0-based, unit-aware time axis, so plotting and overlays line up with the
# times reported in the processed event table (`waveform_max_time`, `event_time_CFD`).

"""
    Waveform

One event's raw waveform: a 0-based, unit-carrying `sample_times` axis
(`0, dt, 2dt, …`) paired with its raw ADC `samples`, plus provenance (`channel`,
`event_id`, `trigger_time`). `samples` is a view into the parent `Waveforms`, so no
data is copied; use [`volts`](@ref) for the samples in physical units.
"""
struct Waveform{T<:AbstractVector, V<:AbstractVector, Q}
    sample_times :: T      # sample times (with units)
    samples      :: V      # raw ADC samples
    channel      :: Int
    event_id     :: Int32
    trigger_time :: Q
end

"""
    Waveforms

Raw waveforms from one channel: the ADC matrix `samples` (`n_samples × n_evts`)
together with per-event `event_id` and `trigger_time`.

Index with `w[i]` to get the i-th event as a [`Waveform`](@ref). `length(w)` is the
number of events, `nsamples(w)` the samples per event, and `sample_times(w)` the
shared 0-based time axis.
"""
struct Waveforms{Q<:AbstractVector}
    samples      :: Matrix{Int16}   # n_samples × n_evts, raw ADC
    channel      :: Int
    event_id     :: Vector{Int32}
    trigger_time :: Q               # per-event, with units
end

Base.length(w::Waveforms)  = size(w.samples, 2)
nsamples(w::Waveforms)     = size(w.samples, 1)

"""
    sample_times(w::Waveforms)

The shared sample-time axis, 0-based and unit-carrying: `0, dt, 2dt, …`.
"""
sample_times(w::Waveforms) = (0:nsamples(w) - 1) .* TIME_PER_SAMPLE

function Base.getindex(w::Waveforms, i::Integer)
    Waveform(sample_times(w), (@view w.samples[:, i]), w.channel, w.event_id[i], w.trigger_time[i])
end

Base.length(wf::Waveform) = length(wf.samples)

"""
    volts(wf::Waveform)

The waveform's ADC samples converted to physical units (volts).
"""
volts(wf::Waveform) = wf.samples .* VOLTS_PER_ADC

function Base.show(io::IO, w::Waveforms)
    print(io, "Waveforms(ch", w.channel, ", ", length(w), " events × ",
          nsamples(w), " samples)")
end

function Base.show(io::IO, wf::Waveform)
    print(io, "Waveform(ch", wf.channel, ", event_id ", wf.event_id, ", ",
          length(wf.samples), " samples)")
end