
"""
    fit_muon_distribution( waveform_max )

Fits a Langau (Landau-Gaussian convolution) to the distribution of peak waveform voltages.

Raises an error if Minuit reports that the fit failed. 
"""
function fit_muon_distribution( waveform_max; verbose=false, Vcut=200mV )

    # chop off the low-voltage background
    data = waveform_max[waveform_max .> Vcut] / 1mV;
    model( μ, σ1, σ2 ) = Langau(μ, σ1, σ2)

    eval_model( x, args... ) = pdf( model(args...), x )
    cost = UnbinnedNLL( data,
        eval_model,
        names = Minuit2.get_argument_names(model)
    )

    m = Minuit( cost, 
        # ~ whatever for initial conditions
        μ=230., σ1=30., σ2=15.,
        strategy=1
    )
    m.limits["σ1", "σ2"] = (1e-6, Inf)
    migrad!( m )

    # print results
    if verbose;
        display( m )
    end

    # raise error if fit fails? 
    if !m.is_valid
        error("minimization failed!")
    end

    return model( m.values... ), m 
end

function find_muon_Vmin( 
    waveform_max; Vmax=500, ΔV=20 )

    data = waveform_max ./ 1mV
    Vmax = Vmax ./ 1mV

    data = data[ data .< Vmax ]
    h = Hist1D( data; binedges=0:ΔV:Vmax )
    bine = binedges(h)
    y = bincounts(h)
    δy = diff(y)

    # find the first positive step, from the end
    j1 = findlast( δy[1:end] .> 0 )

    # find the first negative step, from j1
    j2 = findlast( δy[1:j1] .< 0 )

    return bine[2:(end-1)][j2] * 1mV
end

# function fit_muons_and_bkg( waveform_max )

#     data = t.waveform_max[t.waveform_max .> 10mV] / 1mV

#     model( x, muon_μ, muon_σ1, muon_σ2, bkg_μ, bkg_σ, frac_muon ) = pdf(
#         MixtureModel(
#             [
#                 Langau(muon_μ, muon_σ1, muon_σ2),
#                 Normal(bkg_μ, bkg_σ)
#             ],
#             [ frac_muon, 1 - frac_muon ]
#         ),
#         data
#     ) 

#     cost = UnbinnedNLL( data, model )
#     length( cost.parameters )
#     m = Minuit( cost, 
#         muon_μ=230., muon_σ1=30., muon_σ2=15., 
#         bkg_μ=180, bkg_σ=10., frac_muon=0.95;
#     )
#     m.limits["muon_σ1", "muon_σ2", "bkg_σ"] = (1e-6, Inf)
#     m.limits["frac_muon"] = (0, 1)
#     migrad!( m, strategy=1; ncall=1 )

#     return 
# end


# m_bkg = Gaussian(0, 30)
# m_bkg = MixtureModel([
#     (Gaussian(0, 6), 0.6),
#     (Gaussian(20, 30), 0.4 * 0.95),
#     # (Gaussian(0, 20), 0.9),
#     # (Gaussian(55, 30), 0.1 * 0.85),
#     (Gaussian(120, 40), 0.4 *  0.05)
# ])
# m_bkg = Exponential(25)
