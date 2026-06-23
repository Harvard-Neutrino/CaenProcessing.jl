
# properties of digitizer setup
const VOLTS_PER_ADC = 2u"V" / 2^14
const TIME_PER_SAMPLE = 2u"ns"
const TIME_PER_CLOCK_CYCLE = 8ns #125 mHz clock cycle
const R_LOAD = 50u"Ω" # load resistor

#const n_samples = 1000
const header_size = 6 

const T_HEADER = UInt32
const T_DATA = Int16

"""
    calc_n_evts( fname )
since we know the number of samples + the dtypes, 
we can pre_calc the number of events

NB: could add a check for divisibility? = check on whether file data is somehow non-standard
"""
calc_n_evts( fname, n_samples ) = filesize( fname ) ÷ ( 
    header_size * sizeof(T_HEADER) + n_samples * sizeof(T_DATA) )

# ====================================================

"""
    get_n_event_samples( fname )
Get the number of samples for every event in the file.
Since all events have the same number of samples,
just need to read the header of the first event
"""
function get_n_event_samples( fname )
    open( fname, "r" ) do io
        header = [ read(io, UInt32) for i in 1:6 ]
        return ( header[1] - 24 ) ÷ 2
    end
end

"""
    read_next_trigger!( io, waveform, event_id  )

Read the next trigger data from the binary input stream `io`,updating the values of `waveform` and `event_id`. 
"""
function read_next_trigger!( io, waveform::AbstractArray{Int16},
    event_id, trigger_time, num_overflows )

    try
        header = [ read(io, UInt32) for i in 1:6 ]

        event_size = ( header[1] - 24 ) ÷ 2 
        # board_id = header[2]
        # pattern = header[3]
        # channel = header[4]
        event_counter = header[5]

        trigger_time_tag = header[6]
        # problem: the digitizer saves "trigger_time_tag" as UInt32, which will overflow w/o warning.
        # to solve this, check if the time since prev. trigger is of the same order of mag. as typemax( Int32 ) ≃ 2e9
        t = convert(Float64, trigger_time_tag) + num_overflows[1] * typemax(Int32)
        # and then increment
        if ( trigger_time[1] - t ) > 2e9
            num_overflows[1] += 1
            t += typemax(Int32)
        end

        # if length(waveform) != event_size
        #     error("waveform vector has length $(length(waveform)), but the event length is $(event_size)!")
        # end

        read!(
            io,  
            waveform
        )
        event_id[1] = event_counter
        trigger_time[1] = t
        return true 

    catch error
        if ! isa( error, EOFError ) 
            display(error)
        end
        return false
    end
end

"""
    instantiate_data_vecs( n_samples )

Instantiate vectors for waveforms, event_id, trigger_time, etc. 
to be modified in `read_next_trigger!`. 
"""
function instantiate_data_vecs( n_samples )

    w = Vector{T_DATA}(undef, n_samples)
    event_id = Vector{T_HEADER}( undef, 1 )
    trigger_time = Vector{Float64}( undef, 1 )
    num_overflows = [0]

    return w, event_id, trigger_time, num_overflows
end

"""
    read_waveforms( f_dat; n_evts=Inf )
"""
function read_waveforms( f_dat, n_evts=Inf )

    n_samples = get_n_event_samples( f_dat )

    n_evts = ( isinf(n_evts) ) ? calc_n_evts(f_dat, n_samples) : n_evts 
    waveform_array = Matrix{T_DATA}(undef, n_samples, n_evts )

    # instantiate for while loop
    event_id, trigger_time, num_overflows = instantiate_data_vecs( n_samples )[2:end]

    open( f_dat, "r" ) do io     
    # read + process data evt-by-evt from file   
        for i in axes( waveform_array, 2 )
            read_next_trigger!( io, 
                (@view waveform_array[:, i]), 
                event_id, trigger_time, num_overflows )
        end
    end
    return waveform_array
end


begin
# calculation functions 

    function check_saturation( waveform )
        return any( waveform .== 0 )
    end

    function calc_baseline( waveform, xmin, xmax )
        mean( waveform[xmin:xmax] )
    end

    function calc_baseline_rms( waveform, xmin, xmax )
        rms = sqrt( mean( (waveform[xmin:xmax]).^2 ) )
    end


    """
        calc_CFD_sample( waveform, pct )

    Return the index of the last sample pre-peak for which the waveform is at pct of the peak amplitude. 
    Algorithm = constant fraction discriminator.
    """
    function calc_CFD_sample_index( waveform, pct )

        k = argmax( waveform )
        max_w = waveform[k]

        # go backwards from peak until amplitude = pct * max_w 

        for i in k:-1:1
            ( waveform[i] < max_w * pct ) && return i 
        end

        # no such sample found (?)
        return -1 
    end

end

"""
    process_waveform( waveform; kwargs... )

Compute per-event summary quantities from a single `waveform`.

The baseline and charge-integral windows are anchored to the detected pulse
*onset*, so they adapt to pulses arriving anywhere in the record. All parameters
are keyword arguments with defaults, so no run configuration is required.

Keyword arguments (lengths in samples):
- `cfd_fraction` : onset is the last pre-peak sample below this fraction of the peak height
- `skip_initial` : leave the first this-many samples out of the baseline window
- `base_guard`   : gap between the baseline window and the onset
- `base_width`   : width of the pre-pulse baseline window
- `charge_pre`   : begin the charge integral this many samples before the onset
- `charge_len`   : length of the charge-integral window
- `polarity`     : `"POSITIVE"` or `"NEGATIVE"` pulse direction
"""
function process_waveform( waveform;
        cfd_fraction=0.2, skip_initial=16, base_guard=8, base_width=40,
        charge_pre=8, charge_len=180, polarity="POSITIVE" )

    n = length( waveform )
    sgn = startswith( uppercase(string(polarity)), "POS" ) ? 1.0 : -1.0
    is_saturated = check_saturation( waveform )

    s = sgn .* waveform                       # positive-going pulse

    # locate the pulse using a provisional whole-trace median pedestal
    prov = median( s )
    k = argmax( s )                           # peak sample
    peak = s[k] - prov

    # pulse onset: last pre-peak sample below `cfd_fraction` of the peak height.
    # `onset = firstindex(s) - 1` is the sentinel for "no crossing found".
    onset = firstindex(s) - 1
    for i in k:-1:firstindex(s)
        if ( s[i] - prov ) < cfd_fraction * peak
            onset = i
            break
        end
    end

    # baseline: mean over a `base_width`-sample window placed `base_guard`
    # samples before the onset, kept clear of the first `skip_initial` samples
    bstop = clamp( onset - base_guard - 1, skip_initial + base_width, n )
    bwin  = (bstop - base_width + 1):bstop
    baseline = mean( @view s[bwin] )
    baseline_rms = std( @view s[bwin]; corrected=false, mean=baseline )

    # baseline-subtracted pulse height
    max_w = s[k] - baseline

    # charge integral over a window bracketing the pulse
    cwin = clamp( onset - charge_pre, firstindex(s), n ):clamp( onset + charge_len - 1, firstindex(s), n )
    charge_integral = sum( @view s[cwin] ) - baseline * length(cwin)

    return (
        is_saturated = is_saturated,

        # physical-unit values, with their ADC-unit counterparts (the `_ADC`
        # columns) saved alongside for convenience; the invariant
        # `physical == x_ADC * VOLTS_PER_ADC` holds for each such pair.
        # `peak_ADC` is deliberately separate: it is the literal ADC reading at
        # the peak with the pedestal NOT subtracted, so (unlike waveform_max) it
        # is not just waveform_max expressed in ADC units.
        baseline = sgn * baseline * VOLTS_PER_ADC,
        baseline_ADC = sgn * baseline,
        baseline_rms = baseline_rms * VOLTS_PER_ADC,
        baseline_rms_ADC = baseline_rms,

        waveform_max = max_w * VOLTS_PER_ADC,
        peak_ADC = sgn * ( max_w + baseline ),
        waveform_max_time = ( k - firstindex(s) ) * TIME_PER_SAMPLE,

        event_time_CFD = ( onset - firstindex(s) ) * TIME_PER_SAMPLE,

        charge_integral = charge_integral * VOLTS_PER_ADC / R_LOAD * TIME_PER_SAMPLE,
    )
end


"""
    process_data( f_dat; n_evts=Inf )

Read a .dat file named `f_dat` from the Caen digitizer to a Table object.
"""
function process_data( f_dat; n_evts=Inf, kwargs... )

    n_samples = get_n_event_samples( f_dat )
    n_evts = convert( Int, min( n_evts, calc_n_evts(f_dat, n_samples) ) )

    # instantiate for while loop
    waveform =  Vector{float(T_DATA)}(undef, n_samples)
    w, event_id, trigger_time, num_overflows = instantiate_data_vecs( n_samples )

    # properties to save 
    blank_col(T::Type) = Vector{T}(undef, n_evts)
    blank_col(T::Quantity) = StructArray{T}( undef, n_evts )
    blank_col(T::Type, u::FreeUnits ) = blank_col( typeof(one(T) * u) )
    t = StructArray(
        event_id = blank_col( Int32 ),
        trigger_time = blank_col( Float64, s ),
        is_saturated = blank_col( Bool ),
        
        baseline = blank_col( Float64, mV ),
        baseline_ADC = blank_col( Float64 ),
        baseline_rms = blank_col( Float64, mV ),
        baseline_rms_ADC = blank_col( Float64 ),

        waveform_max_time = blank_col( Int32, ns ),
        # waveform_max_time_raw = blank_col( Int32 ),
        waveform_max = blank_col( Float64, mV ),
        peak_ADC = blank_col( Float64 ),
        
        event_time_CFD = blank_col( Int32, ns ),
        # event_time_CFD_raw = blank_col( Int32 ),

        charge_integral = blank_col( Float64, pC),
    )

    open( f_dat, "r" ) do io     
    # read + process data evt-by-evt from file   
        for i in 1:n_evts

            not_done = read_next_trigger!( io, w, 
                event_id, trigger_time, num_overflows )

            # convert type for calculations
            waveform .= w 

            # calculations
            calcs = process_waveform( waveform; kwargs... )

            # copy to table
            begin
                t.event_id[i] = event_id[1]
                t.trigger_time[i] = trigger_time[1] * TIME_PER_CLOCK_CYCLE

                t.is_saturated[i] = calcs.is_saturated
                t.baseline[i] = calcs.baseline
                t.baseline_ADC[i] = calcs.baseline_ADC
                t.baseline_rms[i] = calcs.baseline_rms
                t.baseline_rms_ADC[i] = calcs.baseline_rms_ADC

                t.waveform_max[i] = calcs.waveform_max
                t.peak_ADC[i] = calcs.peak_ADC
                t.waveform_max_time[i] = calcs.waveform_max_time
                # t.waveform_max_time_raw[i] = calcs.waveform_max_time_raw

                t.event_time_CFD[i]    = calcs.event_time_CFD
                # t.event_time_CFD_raw[i] = calcs.event_time_CFD_raw
                t.charge_integral[i] = calcs.charge_integral
            end
        end
    end

    return t
end