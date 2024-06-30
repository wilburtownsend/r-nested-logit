# Test my nested logit code by comparing simulated shares to actual shares.
source("rnestedlogit.R")
set.seed(1)

# Set parameters.
N = 100000
num_nests = 4
F_per_nest = 1 + rpois(num_nests, 5)
nests = do.call(c, lapply(1:num_nests, function(n) matrix(n, ncol=F_per_nest[n])))
F = length(nests)
delta = rnorm(F)
σ = 0.25
f_by_nest = lapply(1:num_nests, function(n) which(nests==n))

# Draw.
u = rnestedlogit(N, σ, nests, plot_pdf=TRUE)

# Allocate across f.
v = u + t(delta)%x%matrix(1, nrow=N)
f = sapply(1:N, function(i) which.max(v[i,]))
props = sapply(1:F, function(i) mean(f == i))
group_props = sapply(f_by_nest, function(nest) sum(props[nest]))
group_props_of_f = group_props[nests]
conditional_props = props/group_props_of_f

# Find analytic values.
nest_value = sapply(f_by_nest, function(nest) sum(exp(delta[nest]/(1-σ))))
conditional_props_analytic = exp(delta/(1-σ))/nest_value[nests]
group_probs_analytic = nest_value^(1-σ)/sum(nest_value^(1-σ))

# Compare.
print(cbind(conditional_props, conditional_props_analytic))
print(cbind(group_props, group_probs_analytic))