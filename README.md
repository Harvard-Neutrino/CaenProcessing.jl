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

## The processed event table

`process_data` returns a table (one row per event) of high-level summary quantities, some with physical units, as described below:

| Column | Unit | Description | Source |
|---|---|---|---|
| `event_id` | — | digitizer event counter | from header |
| `trigger_time` | s | trigger timestamp (overflow-corrected) | from header |
| `is_saturated` | — | `true` if any sample hit the ADC rail | from raw |
| `baseline` / `baseline_ADC` | mV / ADC | pre-pulse pedestal level | from processing |
| `baseline_rms` / `baseline_rms_ADC` | mV / ADC | noise RMS over the baseline window | from processing |
| `waveform_max` | mV | **baseline-subtracted** pulse height | from processing |
| `peak_ADC` | ADC | **absolute** ADC reading at the peak (pedestal *not* subtracted) | from raw |
| `waveform_max_time` | ns | time of the peak sample | from raw |
| `event_time_CFD` | ns | pulse-start time (constant-fraction discriminator) | from processing |
| `charge_integral` | pC | integrated charge over the pulse window | from processing |

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

### Processing details:

Most of the quantities above come straight from the data's header or the raw samples, but some quantities (in particular `waveform_max` and `charge_integral`) depend on the waveform's baseline, which itself requires first determining the onset of the pulse. 

The pulse onset and the baseline/charge windows are found event-by-event *from the waveform itself*, as follows:

1. **Onset** — reported as `event_time_CFD`. A provisional pedestal (the median
   of the whole trace) is subtracted to estimate the peak height; the onset is
   the last pre-peak sample below `cfd_fraction` of that height (a
   constant-fraction discriminator).
2. **Baseline** — sets `baseline` / `baseline_ADC`, `baseline_rms` /
   `baseline_rms_ADC`, and hence `waveform_max` (the baseline-subtracted peak
   height). The pedestal is the mean over a `base_width`-sample window ending
   `base_guard` samples before the onset, kept clear of the first `skip_initial`
   samples (which can carry a digitizer settling transient); `baseline_rms` is
   the noise RMS over that same window.
3. **Charge** — sets `charge_integral`: the baseline-subtracted samples summed
   from `charge_pre` samples before the onset to `charge_len` samples after it.

The window sizes (in samples) are keyword arguments of the data loading functions, with the following defaults:

| Keyword | Default | Meaning |
|---|---|---|
| `cfd_fraction` | 0.2 | onset threshold, as a fraction of the peak height |
| `skip_initial` | 16 | samples excluded at the start of the record |
| `base_guard` | 8 | gap between the baseline window and the onset |
| `base_width` | 40 | width of the baseline window |
| `charge_pre` | 8 | charge integral starts this many samples before the onset |
| `charge_len` | 180 | length of the charge-integral window |

Pass overrides if needed, e.g:

```julia
process_data("run.dat-ch0"; charge_len=220)     # one channel
read_run("path/to/run_dir"; charge_len=220)     # every channel in the run
```

## `read_run` and the `Run` object:

For data taken using the `caen-run` script from `caen-daq` (https://github.com/kcarloni/caen-daq), the `read_run` convenience function can be used to package together the processed data with its metadata:

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

The `Run` object has the following properties and accessors:

| Access | Meaning |
|---|---|
| `run.config` | parsed `run.toml` (`run.config.bias.voltage_V`, `run.config.channels[ch]`, …) |
| `run.data` | `Dict{Int, table}` — channel number → event table |
| `run[ch]` | shorthand for `run.data[ch]` |
| `keys(run)`, `haskey(run, ch)` | the available channel numbers |


## Digitizer constants

CaenProcessing.jl relies on some fixed properties of the Caen digitizer: 

| Constant | Value | Meaning |
|---|---|---|
| `VOLTS_PER_ADC` | 2 V / 2¹⁴ | ADC count → volts |
| `TIME_PER_SAMPLE` | 2 ns | sampling period |
| `TIME_PER_CLOCK_CYCLE` | 8 ns | trigger-timestamp clock (125 MHz) |
| `R_LOAD` | 50 Ω | load resistor (for charge integral) |

`VOLTS_PER_ADC` and `TIME_PER_SAMPLE` are occasionally useful for converting raw quantities into physical units, and are thus exported.