module BesselFunctions

using AndExport
# using GSL
# using SpecialFunctions
    
# @xport sphj(n, x) = GSL.sf_bessel_jl(n, x)
# @xport sphn(n, x) = GSL.sf_bessel_yl(n, x)
# @xport sphh₊(n, x) = sphj(n, x) + im*sphn(n, x)
# @xport sphh₋(n, x) = sphj(n, x) - im*sphn(n, x)
# # @xport sphj(n, x) = sqrt(π/2x) * besselj(n+0.5, x)
# # @xport sphn(n, x) = sqrt(π/2x) * bessely(n+0.5, x)
# # @xport sphh₊(n, x) = sqrt(π/2x) * besselh1(n+0.5, x)
# # @xport sphh₋(n, x) = sqrt(π/2x) * besselh2(n+0.5, x)

# * Changed to use the complex versions in SpecialFunctions
using SpecialFunctions

@xport sphj(n,x) = sqrt(π/2x) * besselj(n+1/2, x)
@xport sphn(n,x) = sqrt(π/2x) * bessely(n+1/2, x)
@xport sphh₊(n, x) = sqrt(π/2x) * hankelh1(n+1/2, x)
@xport sphh₋(n, x) = sqrt(π/2x) * hankelh2(n+1/2, x)

#####

@xport dsphjdx(n, x) = n/x * sphj(n, x) - sphj(n+1,x)
@xport dsphndx(n, x) = n/x * sphn(n, x) - sphn(n+1,x)
@xport dsphh₊dx(n, x) = n/x * sphh₊(n, x) - sphh₊(n+1,x)
@xport dsphh₋dx(n, x) = n/x * sphh₋(n, x) - sphh₋(n+1,x)

@xport dsphjdr(n,k,r) = dsphjdx(n,k*r) * k
@xport dsphndr(n,k,r) = dsphndx(n,k*r) * k

@xport dsphh₊dr(n,k,r) = dsphh₊dx(n,k*r) * k
@xport dsphh₋dr(n,k,r) = dsphh₋dx(n,k*r) * k

@xport ricj(n,x) = x*sphj(n,x)
# Note: this -ve makes ricn(0,x) = cos(x)
@xport ricn(n,x) = -x*sphn(n,x)
# Note: these negatives are to recover the proper hankel behaviour
# (which is to obtain i^(-n-1) exp(iz) and i^(n+1) exp(-iz)
@xport rich₊(n,x) = ricj(n,x) - im*ricn(n,x)
@xport rich₋(n,x) = ricj(n,x) + im*ricn(n,x)

@xport dricjdx(n,x) = sphj(n,x) + x*dsphjdx(n,x)
@xport dricndx(n,x) = -sphn(n,x) + -x*dsphndx(n,x)
@xport drich₊dx(n,x) = dricjdx(n,x) - im*dricndx(n,x)
@xport drich₋dx(n,x) = dricjdx(n,x) + im*dricndx(n,x)

@xport dricjdr(n,k,r) = k * dricjdx(n,k*r)
@xport dricndr(n,k,r) = k * dricndx(n,k*r)
@xport drich₊dr(n,k,r) = k * drich₊dx(n,k*r)
@xport drich₋dr(n,k,r) = k * drich₋dx(n,k*r)

end
