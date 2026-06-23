
# Utilities for reading `caen-run`-style run directories.
#
# A run directory (produced by the caen-daq `caen-run` wrapper) contains:
#   run.toml        -- machine-parseable metadata (parse THIS, not filenames)
#   used.cfg        -- exact config snapshot passed to wavedump
#   run.dat-chN     -- one binary data file per active channel
#
# `read_run_config` parses run.toml into a NamedTuple; `read_run` parses the
# config and processes every per-channel data file into a table.

import TOML

# regex matching the per-channel data files, e.g. "run.dat-ch0", "run.dat-ch10"
const _CH_FILE_REGEX = r"ch(\d+)$"

"""
    Run

The result of reading a `caen-run` run directory (see [`read_run`](@ref)).

Fields:
- `config` : the parsed `run.toml` (a NamedTuple; see [`read_run_config`](@ref)).
             Per-channel acquisition settings (e.g. `trigger_threshold`) live
             here, under `config.channels[ch]`.
- `data`   : `Dict{Int, table}` mapping each channel number to its processed
             waveform table.

Channel data tables are accessible via `run.data[ch]` or the shorthand `run[ch]`;
`keys(run)` / `haskey(run, ch)` work as well.
"""
struct Run
    config
    data
end

Base.getindex( r::Run, ch::Integer ) = r.data[ch]
Base.keys( r::Run ) = keys( r.data )
Base.haskey( r::Run, ch::Integer ) = haskey( r.data, ch )

# tolerant field access for the config NamedTuple (fields may be absent)
_get( nt, k, default=nothing ) = hasproperty(nt, k) ? getproperty(nt, k) : default

function Base.show( io::IO, r::Run )
    n = length( r.data )
    print( io, "Run(", repr(_get(r.config, :label, "?")), ", ",
        n, n == 1 ? " channel)" : " channels)" )
end

function Base.show( io::IO, ::MIME"text/plain", r::Run )
    c = r.config
    ts = _get( c, :timestamp, nothing )
    println( io, "Run ", repr(_get(c, :label, "?")), ts === nothing ? "" : "  ($ts)" )

    # ── .config branch: a compact subset of the parsed run.toml ──
    println( io, "├─ .config" )
    cfg_lines = String[]
    bias = _get( c, :bias, nothing )
    if bias !== nothing
        v = _get( bias, :voltage_V, nothing )
        v !== nothing && push!( cfg_lines, "bias  = $v V" )
    end
    acq = _get( c, :acquisition, nothing )
    if acq !== nothing
        rt = _get( acq, :run_time_s, nothing )
        rt !== nothing && push!( cfg_lines, "time  = $rt s" )
    end
    notes = _get( c, :notes, nothing )
    (notes !== nothing && !isempty(notes)) && push!( cfg_lines, "notes = $notes" )
    # per-channel trigger thresholds are config, not data
    cfg_chans = _get( c, :channels, nothing )
    if cfg_chans !== nothing && !isempty( cfg_chans )
        thrs = [ "ch$ch=$(_get(cfg_chans[ch], :trigger_threshold, "?"))"
                 for ch in sort(collect(keys(cfg_chans))) ]
        push!( cfg_lines, "thresholds: " * join(thrs, ", ") )
    end
    for l in cfg_lines
        println( io, "│    ", l )
    end

    # ── .data branch: one line per channel (event counts only) ──
    chans = sort( collect(keys(r.data)) )
    print( io, "└─ .data" )
    for ch in chans
        print( io, "\n     ch", ch, " → ", length(r.data[ch]), " events" )
    end
end

"""
    _to_namedtuple(x)

Recursively convert a `Dict{String}` (as returned by `TOML.parsefile`) into a
NamedTuple, so fields are discoverable via `propertynames` and accessible with
dot syntax. Nested tables become nested NamedTuples; arrays are preserved (with
any table elements likewise converted). Non-Dict values pass through unchanged.
"""
function _to_namedtuple(d::AbstractDict)
    pairs = ( Symbol(k) => _to_namedtuple(v) for (k, v) in d )
    return (; pairs...)
end
_to_namedtuple(v::AbstractVector) = map(_to_namedtuple, v)
_to_namedtuple(x) = x

"""
    read_run_config(path) -> NamedTuple

Parse a `run.toml` into a NamedTuple. `path` may be either the run.toml file
itself or a directory containing it.

Fields mirror the TOML structure (e.g. `cfg.bias.voltage_V`, `cfg.acquisition`,
`cfg.label`). The `channels` table is special-cased to a `Dict{Int, NamedTuple}`
keyed by channel number, since its numeric keys ("0", "1", …) are not valid
NamedTuple field names.
"""
function read_run_config( path )
    toml_path = isdir(path) ? joinpath(path, "run.toml") : path
    raw = TOML.parsefile( toml_path )

    # special-case `channels`: numeric string keys -> Dict{Int, NamedTuple}
    channels = nothing
    if haskey(raw, "channels")
        raw = copy(raw)
        chans = pop!(raw, "channels")
        channels = Dict{Int, Any}(
            parse(Int, k) => _to_namedtuple(v) for (k, v) in chans )
    end

    nt = _to_namedtuple( raw )
    return channels === nothing ? nt : merge(nt, (; channels))
end

"""
    read_run(dir; n_evts=Inf) -> Run

Read and process every per-channel data file in a `caen-run` run directory
`dir`. Returns a [`Run`](@ref): its `config` field is the parsed `run.toml` (see
[`read_run_config`](@ref)) and its `data` field is a `Dict{Int, table}` mapping
each channel number to the table returned by [`process_data`](@ref).

The data files are taken from `config.files.data` when present, and otherwise
discovered by globbing `run.dat-ch*` in `dir`.
"""
function read_run( dir; n_evts=Inf, kwargs... )
    config = read_run_config( dir )

    # prefer the authoritative file list in run.toml; fall back to globbing.
    data_files = _run_data_files( config, dir )

    data = Dict{Int, Any}()
    for fname in data_files
        m = match( _CH_FILE_REGEX, fname )
        m === nothing && continue
        ch = parse( Int, m.captures[1] )
        data[ch] = process_data( joinpath(dir, fname); n_evts=n_evts, kwargs... )
    end

    return Run( config, data )
end

"""
    _run_data_files(config, dir) -> Vector{String}

Resolve the list of per-channel data filenames for a run, preferring
`config.files.data` and falling back to globbing `run.dat-ch*` in `dir`
(sorted by channel number).
"""
function _run_data_files( config, dir )
    if hasproperty(config, :files) && hasproperty(config.files, :data)
        return collect( String.(config.files.data) )
    end
    files = filter( f -> occursin(_CH_FILE_REGEX, f), readdir(dir) )
    return sort( files; by = f -> parse(Int, match(_CH_FILE_REGEX, f).captures[1]) )
end
