# Draw nested extreme value random variables:
#   P[X < x] = exp(∑_{c∈C} (∑_{f∈c} exp(-x_f/(1-σ))^(1-σ))

library(fourierin)
library(pracma)
library(extraDistr)

# σ should be a scalar.
# nests should be a length-num_products vector listing each product's nest.
rnestedlogit = function(N, σ, nests,
                        tol=1e-3, plot_pdf=FALSE, plot_imaginary=FALSE,
                        lower_int = -50, upper_int = 50,
                        lower_eval = -5, upper_eval = 10,
                        resolution=2^15
                        ) {
    # Check nests is a vector of integers.
    stopifnot(is.vector(nests))
    stopifnot(all(as.integer(nests) == nests))
    # Assign each nest to a sequential integer.
    num_nests = length(unique(nests))
    xwalk = cbind(unique(nests), 1:num_nests)
    new_nests = sapply(nests, function(nest) xwalk[which(xwalk[,1] == nest), 2])
    # Check that...
    stopifnot(all(new_nests >= 1))
    stopifnot(num_nests == max(new_nests))
    # Count products
    num_products = length(new_nests)
    # Define the characteristic function.
    characteristic = function(z) gammaz(1 - 1i*z)/gammaz(1 - 1i*(1-σ)*z)
    # Integrate to find the PDF.
    out = fourierin(f = characteristic, lower_int = lower_int, upper_int = upper_int,
                    lower_eval = lower_eval, upper_eval = upper_eval,
                    const_adj = -1, freq_adj = -1, resolution = resolution)
    # Plot the pdf, if requested.
    if (plot_pdf) plot(out$w, Re(out$values))
    if (plot_imaginary) plot(out$w, Im(out$values))
    # Check that...
    stopifnot(max(abs(Im(out$values)))                           < tol) # imaginary values are small
    stopifnot(min(Re(out$values))                                > -tol) # pdf is almost non-negative
    stopifnot(max(Re(out$values[1]), Re(out$values[resolution])) < tol) # pdf end points are zero
    # Censor negative and imaginary parts out of the pdf.
    out$values_censored = sapply(out$values, function(v) max(0, Re(v)))
    # Integrate over the pdf to find the CDF.
    out$cdf = cumsum(out$values_censored)
    out$cdf = out$cdf/out$cdf[resolution]
    # Draw using inverse transform sampling, interpolating the inverse CDF.
    u = runif(N*num_nests)
    index_before = sapply(u, function(ui) which(out$cdf < ui)[sum(out$cdf < ui)] )
    index_after  = index_before + 1
    cdf_before = out$cdf[index_before]
    cdf_after  = out$cdf[index_after]
    weight_on_after = (u - cdf_before)/(cdf_after - cdf_before)
    ζ = weight_on_after*out$w[index_after] + (1-weight_on_after)*out$w[index_before]
    # Allocate ζ over draws.
    ζ_reshape = matrix(ζ, nrow=N, ncol=num_nests)
    # for some reason this doesn't work??
    #       outer(1:N, 1:num_products, function(draw, p) ζ_reshape[draw, new_nests[p]])
    ζ_allocated = matrix(nrow=N, ncol=num_products)
    for (p in 1:num_products) {
        ζ_allocated[,p] = ζ_reshape[,new_nests[p]]
    }
    # Now just add the standard GEV.
    gumbels = matrix(rgumbel(num_products*N), nrow=N)
    return((1-σ)*gumbels + ζ_allocated)
}

