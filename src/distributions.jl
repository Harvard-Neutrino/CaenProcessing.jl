
const Gaussian = Normal

# convenience constructor
Distributions.MixtureModel( v::AbstractVector{Tuple{D, T}} ) where {D <: Distribution, T} = MixtureModel( getindex.(v, 1), getindex.(v, 2) )

# ======================================

# Distributions.MixtureModel( models::Vararg{Tuple{D, T}, N}) where {D <: Distribution, T, N} = 

# Much more stripped down than Distributions.MixtureModel,
# struct MixedModel{T,R}
#     components::T
#     weights::AbstractVector{R}
# end
# Distributions.pdf( D::MixedModel, x ) = ( D.weights )' * pdf.( D.components, Ref(x) ) 