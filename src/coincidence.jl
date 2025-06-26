
"""
    get_coincidence_mask( tables::AbstractArray{Table}; V_threshold=160mV )
"""
function get_coincidence_mask( tables::AbstractArray{<:Table}; V_threshold=160mV )
    mask = .&( 
        ( t.waveform_max .> V_threshold for t in tables)...
    )
end

"""
    get_coincidence_mask( files::AbstractArray{AbstractString}; V_threshold=160mV )

Given a set of data files, returns a mask indicating only events with coincident triggers.

For now, the trigger condition is fixed to be
    `evt.waveform_max .> V_threshold`

The digitizer natively triggers all  even (0,2,4,6) or odd (1,3,5,7) channels whenever one channel in the set triggers. Because of this, coincidence checks can be immediately implemented for a subset of files from any given multi-channel data run. 
""" 
function get_coincidence_mask( files; kwargs... )
    get_coincidence_mask( [process_data(f) for f in files]; kwargs... )
end
