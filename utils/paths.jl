
using .CaenProcessing.AstroParticleUnits
using .CaenProcessing.TypedTables

# filename parsing ------

function get_scintillator_num( fname )
    "Z" * match(
        r"Z(?<ZNUM>\d+)",
        fname
    ).captures[1]
end

function get_bias_voltage( fname )
    parse( Float64, match( 
        r"bias-(?<bias>\d+.\d+)",
        fname
    ).captures[1] ) * u"V"
end

function get_temperature( fname )
    parse( Int, match(
        r"(?<temp>-?\d+)C",
        fname
    ).captures[1] ) * u"°C"
end

function get_run_duration( fname )
    parse( Int, match(
        r"(?<duration>\d+)min",
        fname
    ).captures[1] ) * u"minute"
end

function get_digitizer_threshold( fname )
    parse( Int, match(
        r"thres-(?<threshold>\d+)",
        fname
    ).captures[1] )
end

function parse_datetime( fname )
    match(
        r"(?<year>\d+)-(?<month>\d+)-(?<day>\d+)T(?<hour>\d+):(?<minute>\d+):(?<second>\d+)",
        fname
    )
end

function get_timestamp( fname )
    parse_datetime( fname ).match
end

function get_channel( fname )
    parse( Int, match(
        r"_ch(?<channel>\d+)",
        fname
    )["channel"] )
end

# datafile finding ------------

"""
    find_data_files( basedir )

Iterates through `basedir` and finds all the data files, i.e. those ending in `dat`. 
"""
function find_data_files( basedir )
    runs = String[]
    for ( path, dirs, files ) in walkdir( basedir )
        for file in files
            ( endswith(file, ".dat") ) && push!( 
                runs, joinpath(path, file)[(length(basedir)+1):end] )
        end
    end
    return runs
end


"""
    parse_files( files )

Parses the filenames provided and returns a table of the filenames and their properties (timestamp, channel, submodule name, etc.)
"""
function parse_files( files )
    t = Table(
        ZNUM=get_scintillator_num.(files),
        ch=get_channel.(files),
        timestamp=get_timestamp.(files),
        duration=get_run_duration.(files),
        threshold=get_digitizer_threshold.(files),
        bias=get_bias_voltage.(files),
        temp=get_temperature.(files),
        fname=files,
    )
end

"""
    find_and_parse_runs( basedir )

Finds each data run, corresponding to a unique timestamp, and returns a table of the run properties (temperature, bias, etc.)
"""
function find_and_parse_runs( basedir )
    file_table = parse_files( find_data_files(basedir) )

    run_timestamps = unique( file_table.timestamp )
    sort!( run_timestamps )

    ixs = [ findfirst(x -> x == ts, file_table.timestamp) for ts in run_timestamps ]

    @Select(timestamp, temp, bias, duration, threshold)(file_table[ixs])
end




# ===============================

function find_runs_at_timestamp( time_str, basedir )
    runs = String[]
    for ( path, dirs, files ) in walkdir( basedir )
        for file in files
            (occursin( time_str, file ) & endswith(file, ".dat")) && push!( runs, joinpath(path, file) )
        end
    end
    return runs
end

"""
    function find_files_from_run( time_str, basedir, use_even_set=true )

The digitizer natively triggers all  even (0,2,4,6) or odd (1,3,5,7) channels whenever one channel in the set triggers. Because of this, coincidence checks can be immediately implemented for a subset of files from any given multi-channel data run. 
"""
function find_files_from_run( time_str, basedir, use_even_set=true )

    all_runs = find_runs_at_timestamp( time_str, basedir )

    if use_even_set
        return all_runs[ iseven.(get_channel.(all_runs)) ]
        
    else
        return all_runs[ isodd.(get_channel.(all_runs)) ]
        
    end
end


