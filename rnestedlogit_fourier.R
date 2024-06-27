# Draw GEV random variables of the nested logit form:
#	P[X < x] = exp(∑_{c∈C} (∑_{f∈C} exp(-x_f/(1-σ))^(1-σ))
#
# To run:
# 	source("/Users/wtownsend/Documents/GitHub/r-nested-logit/rnestedlogit_fourier.R")

library(fourierin)
library(pracma)

# σ should be a scalar.
# nests should be a length-product vector of each product's nest.
rnestedlogit = function(N, σ, nests, tol=1e-2, plot_pdf=FALSE) {
	# Check nests is a vector of positive integers, with none missing.
	stopifnot(is.vector(nests))
	stopifnot(all(as.integer(nests) == nests))
	stopifnot(all(nests >= 1))
	num_nests    = max(nests)
	num_products = length(nests)
	stopifnot(num_nests == length(unique(nests)))
	# Define the characteristic function.
	characteristic = function(z) gammaz(1 - 1i*z)/gammaz(1 - 1i*(1-σ)*z)
	# Integrate to find the PDF -- should check these parameters.
	resolution = 2^12
	out = fourierin(f = characteristic, lower_int = -20, upper_int = 20,
	                lower_eval = -5, upper_eval = 5,
	                const_adj = -1, freq_adj = -1, resolution = resolution)
	# Plot the pdf, if requested.
	if (plot_pdf) plot(out$w, Re(out$values))
	# Check that...
	stopifnot(max(abs(Im(out$values)))                           < tol) # imaginary values are small
	stopifnot(min(Re(out$values))                                < tol) # pdf is almost non-negative
	stopifnot(max(Re(out$values[1]), Re(out$values[resolution])) < tol) # pdf end points are zero
	# Censor negative and imaginary parts out of the pdf.
	out$values_censored = sapply(out$values, function(v) max(0, Re(v)))
	# Integrate over the pdf to find the CDF.
	out$cdf = cumsum(out$values_censored)
	out$cdf = out$cdf/out$cdf[resolution]
	# Draw.
	u = runif(N*num_nests)
	index_before = sapply(u, function(ui) which(out$cdf < ui)[sum(out$cdf < ui)] )
	index_after  = index_before + 1
	cdf_before = out$cdf[index_before]
	cdf_after  = out$cdf[index_after]
	weight_on_after = (u - cdf_before)/(cdf_after - cdf_before)
	ζ = weight_on_after*out$w[index_after] + (1-weight_on_after)*out$w[index_before]
	# Allocate over draws.
	ζ_reshape = matrix(ζ, nrow=N, ncol=num_nests)
	# for some reason this doesn't work??
	#		outer(1:N, 1:num_products, function(draw, p) ζ_reshape[draw, nests[p]])
	ζ_allocated = matrix(nrow=N, ncol=num_products)
	for (p in 1:num_products) {
		ζ_allocated[,p] = ζ_reshape[,nests[p]]
	}
	# Now just add the standard GEV.
	gumbels = matrix(rgumbel(num_products*N), nrow=N)
	return(gumbels + ζ_allocated)
}


set.seed(1)

# Test.
N = 1000000
num_nests = 4
F_per_nest = 1 + rpois(num_nests, 5)
nests = do.call(c, lapply(1:num_nests, function(n) matrix(n, ncol=F_per_nest[n])))
F = length(nests)
delta = rnorm(F)
σ = 0.25

f_by_nest = lapply(1:num_nests, function(n) which(nests==n))

u = rnestedlogit(N, σ, nests, plot_pdf=TRUE)
v = u + t(delta)%x%matrix(1, nrow=N)
f = sapply(1:N, function(i) which.max(v[i,]))
props = sapply(1:F, function(i) mean(f == i))


group_props = sapply(f_by_nest, function(nest) sum(props[nest]))
group_props_of_f = group_props[nests]
conditional_props = props/group_props_of_f

# Find analytic values...
nest_value = sapply(f_by_nest, function(nest) sum(exp(delta[nest]/(1-σ))))
conditional_props_analytic = exp(delta/(1-σ))/nest_value[nests]
group_probs_analytic = nest_value^(1-σ)/sum(nest_value^(1-σ))

# Compare.
print(cbind(conditional_props, conditional_props_analytic))
print(cbind(group_props, group_probs_analytic))


