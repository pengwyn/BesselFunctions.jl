module BesselFunctions

using AndExport
# using GSL
# using SpecialFunctions
    
# @AndExport sphj(n, x) = GSL.sf_bessel_jl(n, x)
# @AndExport sphn(n, x) = GSL.sf_bessel_yl(n, x)
# @AndExport sphh₊(n, x) = sphj(n, x) + im*sphn(n, x)
# @AndExport sphh₋(n, x) = sphj(n, x) - im*sphn(n, x)
# # @AndExport sphj(n, x) = sqrt(π/2x) * besselj(n+0.5, x)
# # @AndExport sphn(n, x) = sqrt(π/2x) * bessely(n+0.5, x)
# # @AndExport sphh₊(n, x) = sqrt(π/2x) * besselh1(n+0.5, x)
# # @AndExport sphh₋(n, x) = sqrt(π/2x) * besselh2(n+0.5, x)

# * Changed to use the complex versions in SpecialFunctions
using SpecialFunctions

@AndExport sphj(n,x) = sqrt(π/2x) * besselj(n+1/2, x)
@AndExport sphn(n,x) = sqrt(π/2x) * bessely(n+1/2, x)
@AndExport sphh₊(n, x) = sqrt(π/2x) * hankelh1(n+1/2, x)
@AndExport sphh₋(n, x) = sqrt(π/2x) * hankelh2(n+1/2, x)

#####

@AndExport dsphjdx(n, x) = n/x * sphj(n, x) - sphj(n+1,x)
@AndExport dsphndx(n, x) = n/x * sphn(n, x) - sphn(n+1,x)
@AndExport dsphh₊dx(n, x) = n/x * sphh₊(n, x) - sphh₊(n+1,x)
@AndExport dsphh₋dx(n, x) = n/x * sphh₋(n, x) - sphh₋(n+1,x)

@AndExport dsphjdr(n,k,r) = dsphjdx(n,k*r) * k
@AndExport dsphndr(n,k,r) = dsphndx(n,k*r) * k

@AndExport dsphh₊dr(n,k,r) = dsphh₊dx(n,k*r) * k
@AndExport dsphh₋dr(n,k,r) = dsphh₋dx(n,k*r) * k

@AndExport ricj(n,x) = x*sphj(n,x)
# Note: this -ve makes ricn(0,x) = cos(x)
@AndExport ricn(n,x) = -x*sphn(n,x)
# Note: these negatives are to recover the proper hankel behaviour
# (which is to obtain i^(-n-1) exp(iz) and i^(n+1) exp(-iz)
@AndExport rich₊(n,x) = ricj(n,x) - im*ricn(n,x)
@AndExport rich₋(n,x) = ricj(n,x) + im*ricn(n,x)

@AndExport dricjdx(n,x) = sphj(n,x) + x*dsphjdx(n,x)
@AndExport dricndx(n,x) = -sphn(n,x) + -x*dsphndx(n,x)
@AndExport drich₊dx(n,x) = dricjdx(n,x) - im*dricndx(n,x)
@AndExport drich₋dx(n,x) = dricjdx(n,x) + im*dricndx(n,x)

@AndExport dricjdr(n,k,r) = k * dricjdx(n,k*r)
@AndExport dricndr(n,k,r) = k * dricndx(n,k*r)
@AndExport drich₊dr(n,k,r) = k * drich₊dx(n,k*r)
@AndExport drich₋dr(n,k,r) = k * drich₋dx(n,k*r)

end
