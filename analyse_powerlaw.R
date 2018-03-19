# install.packages("devtools")
# devtools::install_github("csgillespie/poweRlaw", subdir="pkg")

library("poweRlaw")

test <- relvar$MTE
test <- test[which(test>0)]

## continuous power-law
d_cpl = conpl$new(test)
plot(d_cpl)

## continuous log-normal
d_cln = conlnorm$new(test)
plot(d_cln)
lines(d_cln)

# ## test if data is drawn from a power law distribution
# bs_cpl = bootstrap_p(d_cpl, no_of_sims=1000, threads=2)
