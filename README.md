# CaenProcessing.jl

Read and process binary `.dat` waveform files from our CAEN digitizer (written by
`wavedump` / the `caen-run` acquisition wrapper) into tidy, unit-aware tables.

Install using Julia's package manager: 
```julia
pkg> add https://github.com/kcarloni/CaenProcessing.jl.git
```

## Quick start

```julia
using CaenProcessing

# process one channel's data file -> a table, one row per event
t = process_data("run.dat-ch0")

# or read a whole caen-run directory (all channels + metadata) at once
run = read_run("path/to/run_dir")
run.config            # parsed run.toml (NamedTuple)
run[0]                # channel 0's table  (shorthand for run.data[0])
```

## Reading data

| Function | Returns | Notes |
|---|---|---|
| `read_waveforms(f_dat, n_evts=Inf)` | `Matrix{Int16}` (`n_samples × n_evts`) | raw ADC samples, no processing |
| `process_data(f_dat; n_evts=Inf)` | table (one row per event) | computes the quantities below |
| `read_run_config(path)` | `NamedTuple` | parse a `run.toml` (file *or* directory) |
| `read_run(dir; n_evts=Inf)` | [`Run`](#the-run-object) | parse config + process every channel |

`n_evts` caps how many events are read (default: all).

## The event table

`process_data` returns a `StructArray` (from
[StructArrayTables.jl](https://github.com/kcarloni/StructArrayTables.jl)); access
a column with `t.colname`. Physical-unit columns carry
[AstroParticleUnits](https://github.com/kcarloni/AstroParticleUnits.jl) units.

| Column | Unit | Description |
|---|---|---|
| `event_id` | — | digitizer event counter |
| `trigger_time` | s | trigger timestamp (overflow-corrected) |
| `is_saturated` | — | `true` if any sample hit the ADC rail |
| `baseline` / `baseline_ADC` | mV / ADC | pre-pulse pedestal level |
| `baseline_rms` / `baseline_rms_ADC` | mV / ADC | noise RMS over the baseline window |
| `waveform_max` | mV | **baseline-subtracted** pulse height |
| `peak_ADC` | ADC | **absolute** ADC reading at the peak (pedestal *not* subtracted) |
| `waveform_max_time` | ns | time of the peak sample |
| `event_time_CFD` | ns | pulse-start time (constant-fraction discriminator) |
| `charge_integral` | C | integrated charge over the pulse window |

**Raw vs. physical.** For convenience some amplitude quantities are provided in both
physical units and raw ADC counts. The `_ADC` columns are exactly the
pre-scaling values, so
```
physical == x_ADC * VOLTS_PER_ADC
```

holds for every `_ADC` pair (e.g. `baseline == baseline_ADC * VOLTS_PER_ADC`).

`peak_ADC` is deliberately *separate* from this rule: it is the literal ADC value
at the peak with the pedestal left in, so unlike `waveform_max` it is **not** just
the pulse height in ADC units. It's handy for saturation / DC-offset checks.

## The `Run` object

`read_run` returns a `Run` wrapping the run's metadata and per-channel data:

| Access | Meaning |
|---|---|
| `run.config` | parsed `run.toml` (`run.config.bias.voltage_V`, `run.config.channels[ch]`, …) |
| `run.data` | `Dict{Int, table}` — channel number → event table |
| `run[ch]` | shorthand for `run.data[ch]` |
| `keys(run)`, `haskey(run, ch)` | the available channel numbers |

It prints a compact summary:

```
Run "bias31V_5min"  (20260622T175105)
├─ .config
│    bias  = 31.0 V
│    time  = 300 s
│    thresholds: ch0=1200, ch1=1170
└─ .data
     ch0 → 467552 events
     ch1 → 467552 events
```

Per-channel trigger thresholds come from the config (`run.config.channels[ch]`),
not the data — `run.data` holds only the recorded waveforms.

### `caen-run` directory layout

`read_run` expects a directory produced by the `caen-run` script (https://github.com/kcarloni/caen-daq):

```
run_dir/
    run.toml        # metadata: label, bias, acquisition, per-channel config
    used.cfg        # exact wavedump config snapshot
    run.dat-ch0     # one binary data file per active channel
    run.dat-ch1
```

The channel files are taken from `run.toml`'s `[files].data` list when present,
otherwise discovered by globbing `run.dat-ch*`.

## Digitizer constants

Setup constants live in `read_and_process.jl`; `VOLTS_PER_ADC` and
`TIME_PER_SAMPLE` are exported:

| Constant | Value | Meaning |
|---|---|---|
| `VOLTS_PER_ADC` | 2 V / 2¹⁴ | ADC count → volts |
| `TIME_PER_SAMPLE` | 2 ns | sampling period |
| `TIME_PER_CLOCK_CYCLE` | 8 ns | trigger-timestamp clock (125 MHz) |
| `R_LOAD` | 50 Ω | load resistor (for charge integral) |
