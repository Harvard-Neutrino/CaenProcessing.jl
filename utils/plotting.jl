
# module CaenPlotting
using Printf: @sprintf 

using FHist
using .CaenProcessing.UnitfulDistributions
using .CaenProcessing.Distributions
using .CaenProcessing.AstroParticleUnits

function plot_dist_Vmax!( ax, t; bine=nothing,
    color=:black, markersize=4, label=nothing, )

    ax.xlabel = "max voltage of pulse [mV]"
    ax.ylabel = "counts"

    data = t.waveform_max / 1mV;
    bine = isnothing(bine) ? (minimum(data):20:maximum(data)) : bine 
    h = Hist1D( data; binedges=bine )

    x = bincenters(h)
    y = bincounts(h)
    y_err = sqrt.(y)
    # norm = sum(y) * step(bine)
    ixs = (y .> 0)
    x, y, y_err = x[ixs], y[ixs], y_err[ixs]

    y_err_low = copy(y_err)
    y_err_low[ y .== 1] .= (1 - 1e-5)

    scatter!( ax, x, y;
        color,
        markersize,
        label
        )
    errorbars!( ax, x, y, y_err_low, y_err;
        color
        )
    #
end

function plot_dist_Vmax( t; fname, kwargs... )

    fig = Figure( size=(4inch, 3inch), fontsize=10pt)
    ax = Axis( fig[1,1],
        xgridvisible=false,
        ygridvisible=false,
        yscale=log10
    )

    plot_dist_Vmax!( ax, t; kwargs... )
    !isnothing(fname) && Label( fig[1,1,Top()], text=fname, fontsize=6pt )

    return fig, ax 
end

function plot_muon_fit!( ax, m, errs; 
    norm, Vcut, 
    xmax=1600,
    color, label="" )

    # x = minimum(data):1:maximum(data)
    Vcut /= 1mV
    x = Vcut:xmax

    lines!( ax, x, pdf(m, x) * norm;    
        color, linewidth=2, alpha=0.5,
        label=label * pretty_repr(m * 1.0mV, errs),

    )
    vlines!( ax, Vcut; 
        color, linestyle=(:dash, 2),
        linewidth=1, alpha=0.5
    )

end

# -----------------

function pretty_repr( m::Distribution )
    s = "$(typeof(m).name.wrapper):\n"
    for sym in propertynames(m)
        s *= "$(string(sym)) = " * (@sprintf "%.1f" getproperty(m, sym)) * "\n"
    end
    return s
end
function pretty_repr( m::UnitfulDistribution{Q} ) where Q
    s = "$(typeof(m.D).name.wrapper):\n"
    for sym in propertynames(m.D)
        x = getproperty(m, sym)
        s *= "$(string(sym)) = " * (@sprintf "%.1f" x/oneunit(x)) * " $(unit(x))\n"
    end
    return s
end
function pretty_repr( m::Distribution, errs )
    s = "$(typeof(m).name.wrapper):\n"
    for (i, sym) in enumerate(propertynames(m))
        s *= "$(string(sym)) = " * (@sprintf "%.1f ± %.1f" getproperty(m, sym) errs[i]) * "\n"
    end
    return s
end
function pretty_repr( m::UnitfulDistribution{Q}, errs ) where Q
    s = "$(typeof(m.D).name.wrapper):\n"
    for (i, sym) in enumerate( propertynames(m.D) )
        x = getproperty(m, sym)
        s *= "$(string(sym)) = " * (@sprintf "%.1f ± %.1f" x/oneunit(x) errs[i]) * " $(unit(x))\n"
    end
    return s
end

# end