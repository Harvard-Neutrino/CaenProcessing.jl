
using Distributions: @check_args

"""
    Langau( μ, σ_lan, σ_gau )

Distribution given by the convolution of a Gaussian kernel over a Landau.
Inspired by the ROOT implementation, see: https://root.cern/doc/v636/langaus_8C.html. 
"""
struct Langau{T<:Real} <: Distributions.ContinuousUnivariateDistribution
    μ::T
    σ_lan::T
    σ_gau::T
    Langau{T}( μ::T, σ_lan::T, σ_gau::T ) where {T <: Real} = new{T}(μ, σ_lan, σ_gau)
end

function Langau( μ::T, σ_lan::T, σ_gau::T; check_args::Bool=true ) where {T<:Real}
    @check_args Langau (σ_lan, σ_lan > zero(σ_lan)) (σ_gau, σ_gau > zero(σ_gau)) 
    return Langau{T}(μ, σ_lan, σ_gau)
end
Langau( μ::Real, σ_lan::Real, σ_gau::Real; check_args::Bool=true) = Langau(
    promote(μ, σ_lan, σ_gau)...; check_args=check_args
)

Base.broadcastable( D::Langau ) = Ref(D)
Distributions.StatsAPI.params( D::Langau ) = ( D.μ, D.σ_lan, D.σ_gau )

# using DSP 

# # fast calculation for ranges, using DSP.conv
# # NVM: requires some finicking for the case where step( x ) ≫ σ_gau. 
# function Distributions.pdf( 
#     D::Langau{Q}, x::AbstractRange{Q} ) where {Q}
#     lan = Landau( D.μ, D.σ_lan )
#     Δx = step( x )

#     # compute padding factor 
#     n_pad = 5 * convert( Int, div( D.σ_gau, Δx, RoundUp ) )

#     x_min = minimum(x) - n_pad * Δx 
#     x_max = maximum(x) + n_pad * Δx
    
#     u_x = oneunit(Q)
#     f1 = pdf( lan, x_min:Δx:x_max ) 
#     f2 = pdf( 
#         Normal(0*u_x, D.σ_gau), 
#         (-n_pad*Δx):Δx:(n_pad*Δx) 
#     ) * Δx

#     out = ( conv( f1 * u_x, f2 ) / u_x )[(1+2n_pad):(end-2n_pad)]
# end

# simple dot-product for single points
function Distributions.pdf( 
    D::Langau, x::Real )

    lan = Landau( D.μ, D.σ_lan )

    n = 50
    Δx = D.σ_gau / n
    n_pad = 5 * n 

    u_x = 1.0 #oneunit(Q)
    v = ( -n_pad * Δx ):Δx:( n_pad * Δx )

    f1 = pdf( lan, ( x .+ v ) )
    f2 = pdf( Normal( 0 * u_x, D.σ_gau ), v ) * Δx 

    return f1' * f2 
end