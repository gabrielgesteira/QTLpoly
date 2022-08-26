## Installing packages
## install.packages("qtlpoly")
##install.packages("mappoly")
## ##install.packages("devtools")
## devtools::install_github("gabrielgesteira/qtlpoly")

## Loading packages
library(qtlpoly) 
library(mappoly)
plot_map_list(maps4x)
plot(maps4x[[1]])
summary_maps(maps4x)
head(pheno4x)

## Calculating genotype probabilities
genoprob4x = lapply(maps4x, calc_genoprob)

## Putting geno + pheno together
data = read_data(ploidy = 4, geno.prob = genoprob4x, pheno = pheno4x, step = 1) 

## Viewing dataset
print(data, detailed = TRUE) 

## Adjusting null model 
null.mod = null_model(data = data, n.clusters = 4)

## Viewing null model
print(null.mod)

## Performing first search round
search.mod = search_qtl(data = data, model = null.mod, w.size = 15, sig.fwd = 0.01, n.clusters = 4)

## Visualizing results
print(search.mod)

## Performing model optimization
optimize.mod = optimize_qtl(data = data, model = search.mod, sig.bwd = 0.0001, n.clusters = 4)

## Visualizing results
print(optimize.mod)

## Performing another search round
search.mod2 = search_qtl(data=data, model = optimize.mod, sig.fwd = 0.0001, n.clusters = 4)

## Visualizing results
print(search.mod2)

## Genomewide profiling with final set of QTL
profile.mod = profile_qtl(data = data, model = optimize.mod, d.sint = 1.5, polygenes = FALSE, n.clusters = 4)

## Visualizing results
print(profile.mod)

## Checking lower SI
print(profile.mod, sint = "lower") 

## Checking upper SI
print(profile.mod, sint = "upper")

## Performing automatic search
remim.mod = remim(data = data, w.size = 15, sig.fwd = 0.01, sig.bwd = 0.0001, d.sint = 1.5, n.clusters = 4)

## Visualizing results
print(remim.mod)

## Checking lower SI
print(remim.mod, sint = "lower")

## Checking upper SI
print(remim.mod, sint = "upper")

## Plotting profile
for(p in remim.mod$pheno.col) plot_profile(data = data, model = remim.mod, pheno.col = p, ylim = c(0,10))

## Plotting profile
plot_profile(data = data, model = remim.mod, grid = TRUE) 
plot_profile(data = data, model = remim.mod, grid = FALSE) 

## Plotting QTL SI
plot_sint(data = data, model = remim.mod)

## Calculate metrics
fitted.mod = fit_model(data=data, model=remim.mod, probs="joint", polygenes="none")

## Viewing results
summary(fitted.mod)

## Plotting QTL info
plot_qtl(data = data, model = remim.mod, fitted = fitted.mod, drop.pheno = FALSE)

## Account for allele effects
est.effects = qtl_effects(ploidy = 4, fitted = fitted.mod)

## Plot allele effects
plot(est.effects) 

## Calculate breeding values
y.hat = breeding_values(data = data, fitted = fitted.mod)

## Plot BV
plot(y.hat) 

## FEIM
## Perform permutations
perm = permutations(data = data, n.sim = 1000, n.clusters = 4)

## Check permutations
print(perm)

## Define significance threshold
(sig.lod = perm$sig.lod$`0.95`)

## Adjust FEIM model
feim.mod = feim(data = data, w.size = 15, sig.lod = sig.lod)

## Viewing FEIM results
print(feim.mod)

## Plotting FEIM results
plot_profile(data = data, model = feim.mod, grid = TRUE)

## Exporting results to VIEWpoly
## save(maps4x, file="mappoly.maps.RData")
## save(data, file="qtlpoly.data.RData")
## save(remim.mod, file="qtlpoly.remim.mod.RData")
## save(fitted.mod, file="qtlpoly.fitted.mod.RData")
## save(est.effects, file="qtlpoly.est.effects.RData")
