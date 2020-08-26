module BesselFunctions

using AndExport
# using GSL
# using SpecialFunctions
    
# @export sphj(n, x) = GSL.sf_bessel_jl(n, x)
# @export sphn(n, x) = GSL.sf_bessel_yl(n, x)
# @export sphh₊(n, x) = sphj(n, x) + im*sphn(n, x)
# @export sphh₋(n, x) = sphj(n, x) - im*sphn(n, x)
# # @export sphj(n, x) = sqrt(π/2x) * besselj(n+0.5, x)
# # @export sphn(n, x) = sqrt(π/2x) * bessely(n+0.5, x)
# # @export sphh₊(n, x) = sqrt(π/2x) * besselh1(n+0.5, x)
# # @export sphh₋(n, x) = sqrt(π/2x) * besselh2(n+0.5, x)

# * Changed to use the complex versions in SpecialFunctions
using SpecialFunctions

@export sphj(n,x) = sqrt(π/2x) * besselj(n+1/2, x)
@export sphn(n,x) = sqrt(π/2x) * bessely(n+1/2, x)
@export sphh₊(n, x) = sqrt(π/2x) * hankelh1(n+1/2, x)
@export sphh₋(n, x) = sqrt(π/2x) * hankelh2(n+1/2, x)

#####

@export dsphjdx(n, x) = n/x * sphj(n, x) - sphj(n+1,x)
@export dsphndx(n, x) = n/x * sphn(n, x) - sphn(n+1,x)
@export dsphh₊dx(n, x) = n/x * sphh₊(n, x) - sphh₊(n+1,x)
@export dsphh₋dx(n, x) = n/x * sphh₋(n, x) - sphh₋(n+1,x)

@export dsphjdr(n,k,r) = dsphjdx(n,k*r) * k
@export dsphndr(n,k,r) = dsphndx(n,k*r) * k

@export dsphh₊dr(n,k,r) = dsphh₊dx(n,k*r) * k
@export dsphh₋dr(n,k,r) = dsphh₋dx(n,k*r) * k

@export ricj(n,x) = x*sphj(n,x)
# Note: this -ve makes ricn(0,x) = cos(x)
@export ricn(n,x) = -x*sphn(n,x)
# Note: these negatives are to recover the proper hankel behaviour
# (which is to obtain i^(-n-1) exp(iz) and i^(n+1) exp(-iz)
@export rich₊(n,x) = ricj(n,x) - im*ricn(n,x)
@export rich₋(n,x) = ricj(n,x) + im*ricn(n,x)

@export dricjdx(n,x) = sphj(n,x) + x*dsphjdx(n,x)
@export dricndx(n,x) = -sphn(n,x) + -x*dsphndx(n,x)
@export drich₊dx(n,x) = dricjdx(n,x) - im*dricndx(n,x)
@export drich₋dx(n,x) = dricjdx(n,x) + im*dricndx(n,x)

@export dricjdr(n,k,r) = k * dricjdx(n,k*r)
@export dricndr(n,k,r) = k * dricndx(n,k*r)
@export drich₊dr(n,k,r) = k * drich₊dx(n,k*r)
@export drich₋dr(n,k,r) = k * drich₋dx(n,k*r)

end
