using LegendrePolynomials, Printf
import OffsetArrays: no_offset_view

# adaptively compute analytical covariances from power spectral density
function sphere_cov(g, t; rtol=eps(), atol=0, blksz=100, maxterms=10000, verb=false)
    cov  = 0.0
    l0   = 0
    diff = 0.0
    converged = false
    while !converged && l0 < maxterms 
        l = l0:l0+blksz-1
        diff = sum(
            (2l.+1)/(4pi) .* g.(l.*(l.+1)).^2 .* no_offset_view(
                collectPl.(cos(t), lmin=l[1], lmax=l[end])
                )
            )
        cov += diff

        converged = (abs((diff / blksz) / cov) < rtol) || (abs(diff / blksz) < atol)
        l0 += blksz
    end

    if verb > 0 && l0 >= maxterms
        @printf "t=%.2f did not converge to rtol=%.2e atol=%.2e after max %i terms, remaining terms of size %.2e (relative %.2e)\n" t rtol atol maxterms abs(diff / blksz) abs((diff / blksz) / cov)
    elseif verb > 1
        @printf "used %i terms with block size %i, rtol=%.2e atol=%.2e acheived\n" l0 blksz rtol atol 
    end

    return cov
end