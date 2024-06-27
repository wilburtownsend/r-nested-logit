# Draw GEV random variables of the nested logit form:
#	P[X < x] = exp(∑_{c∈C} (∑_{f∈C} exp(-x_f/λ))^λ)
# Convenience wrapper for rmvevd() from 
# https://search.r-project.org/CRAN/refmans/evd/html/mvevd.html
# 	xx citation("evd")
#
# nests should be a disjoint, exhaustive list of nests.
# λ should be a scalar.
# 	source("/Users/wtownsend/Documents/GitHub/r-nested-logit/rnestedlogit EVD.R")
# library(evd)
library(fourierin)

rnestedlogit = function(N, λ, nests, as_in_my_model=FALSE) {
	# Check nests 
	all_nests = do.call(c, nests)
	F = max(all_nests)
	stopifnot(length(all_nests) == F)
	stopifnot(min(all_nests) == 1)
	stopifnot(sort(unique(all_nests)) == sort(all_nests))
	# Ensure all nests are sorted.
	nests = lapply(nests, function(nest) sort(nest))
	# Create asy.
	all_subsets = do.call(c, lapply(1:F, function(l) combn(F, l, simplify = FALSE)))
	num_subsets = length(all_subsets)
	is_nest = sapply(all_subsets, function(subset) any(sapply(nests, function(nest) length(subset) == length(nest) && all(subset == nest))))
	asy = lapply(1:num_subsets, function(i) if (is_nest[i]) all_subsets[[i]]*0+1 else all_subsets[[i]]*0)
	# Create dep.
	dep = rep(λ, 2^F-F-1)
	deviates = rmvevd(N, dep = dep, asy = asy, model = "alog", d = F)
	# If as_in_my_model, we do the convention that we divide by lambda before returning
	if (as_in_my_model) deviates = deviates / λ
	return(deviates)
}



# Test.
N = 1
num_nests = 10
F_per_nest = 1 + rpois(num_nests, 9)
nests = list()
f = 1
for (n in 1:num_nests) {
	nests = c(nests, list(seq(from=f, to = f + F_per_nest[n]-1)))
	f = f + F_per_nest[n]
}
#nests = list(c(1), c(2,3), c(4,6,7,8), c(5, 9), c(10, 11))
F = max(do.call(c, nests))
delta = rnorm(F)
λ = 0.25
nest_of_f = sapply(1:F, function(f) which(sapply(nests, function(nest) f %in% nest)))
u = rnestedlogit(N, λ, nests, TRUE)
stop
v = u + t(delta)%x%matrix(1, nrow=N)
f = sapply(1:N, function(i) which.max(v[i,]))
props = sapply(1:F, function(i) mean(f == i))
group_props = sapply(nests, function(nest) sum(props[nest]))
group_props_of_f = group_props[nest_of_f]
conditional_props = props/group_props_of_f
print(conditional_props)
print(log(conditional_props[4]) - log(conditional_props[2]))

# Now find group props analytically...
delta_group = sapply(nests, function(nest) λ*log(sum(exp(delta[nest]))))
group_probs_analytic = exp(delta_group)/sum(exp(delta_group))
print(cbind(group_props, group_probs_analytic))



