setwd("~/Dropbox/pmdvlab/xoi/")

library(tidyverse)
library(ggplot2)

source("helpers/xo_helpers.R")
source("helpers/gamma_mcmc_2.R")

## set global ggplot2 theme
old_theme <- theme_set(theme_minimal())

## function to process MCMC result into object that returns both estimates in terms of original covariates
## (like `predict()` for a linear model) and the actual coefficient estimates (to allow post-hoc contrasts)
process_mcmc_result <- function(ff, chains = LETTERS[1:4], covars = c("age","cross"), ...) {
	
	coefs <- tibble()
	df <- tibble()
	
	for (cc in chains) {
		ff_f <- paste0(ff, "_", cc, ".rds")
		message("reading ", ff_f, " ...")
		xx <- readRDS(ff_f)
		colnames(xx$trace) <- xx$data$coefs
		
		for (ii in 1:length(xx$data$data)) {
			d <- xx$data$data[[ii]][1,]
			g <- as.vector(xx$trace %*% xx$data$design[[ii]][1,])
			this_df <- tibble(nu = g, .iter = seq_along(g), .chain = cc)
			for (covar in covars) {
				this_df[[covar]] <- d[[covar]][1]
			}
			df <- bind_rows(df, this_df)
		}
		
		coefs <- bind_rows(coefs,
						   as_tibble(xx$trace) %>% mutate(.iter = 1:nrow(xx$trace), .chain = cc))
	}
	
	return(list(predicted = df,
				coefs = coefs))
	
}

## function to calculate pointwise log-likelihood from MCMC chain and the input data
## (this is needed to calculate WAIC and similar quantities, for formal model comparison)
calc_pointwise_loglik <- function(x, idx = seq(5001,25000,40), ...) {

	message("Preparing data ...")
	XX <- Reduce(rbind, x$data$design)
	data <- Reduce(rbind, x$data$data)
	nobs <- nrow(XX)
	nsamples <- length(idx)
	L <- matrix(0, nrow = nsamples, ncol = nobs)
	
	message("Calculating pointwise log-likelihoods for ", nobs, " observations x ", nsamples, " posterior samples ...")
	
	pb <- txtProgressBar(min = 0, max = nsamples, style = 3)
	for (ii in seq_len(nsamples)) {
		nu <- XX %*% x$trace[ idx[ii], ] %>% as.vector() %>% exp()
		#print(unique(nu))
		for (jj in seq_len(nobs)) {
			lambda <- loglik( nu[jj], data[jj,] )
			L[ii,jj] <- lambda
		}
		setTxtProgressBar(pb, ii)
	}
	
	message("Done.\n")
	return(L)
	
}

## calculate pointwise log-likelihoods
ff <- paste0("gamma_fits_agecross/out_sigma005_", LETTERS[1:4], ".rds")
raw_rez <- lapply(ff, readRDS)
pw <- lapply(raw_rez, calc_pointwise_loglik)

## calculate pointwise log-likelihoods UNDER NULL MODEL
ff <- paste0("gamma_null_agecross/out_sigma005_", LETTERS[1:4], ".rds")
raw_rez <- lapply(ff, readRDS)
pw_null <- lapply(raw_rez, calc_pointwise_loglik)

## save pointwise likelihoods for later
pointwise <- list(full = pw, null = pw_null)
saveRDS(pointwise, "gamma_fits_agecross/pointwise_ll_sigma005.rds")

waic_full <- loo::loo(Reduce(rbind, pointwise$full))
waic_null <- loo::loo(Reduce(rbind, pointwise$null))
loo::loo_compare(waic_null, waic_full)

## make table of ELPD values for full and null models
list(full = waic_full, null = waic_null) %>%
	lapply(., "[[", "estimates") %>%
	lapply(., as_tibble, rownames = "term", .name_repair = tolower) %>%
	bind_rows(.id = "model") %>% 
	subset(term == "elpd_loo")

## how many highly influential obs under each model
table( waic_null$diagnostics$pareto_k <= 0.5 ) %>% prop.table()
table( waic_full$diagnostics$pareto_k <= 0.5 ) %>% prop.table()

list(full = waic_full, null = waic_null) %>%
	lapply(., "[[", "estimates") %>%
	lapply(., as_tibble, rownames = "term", .name_repair = tolower) %>%
	bind_rows(.id = "model") %>% 
	subset(term == "elpd_loo")

## load MCMC results
rez <- process_mcmc_result("gamma_fits_agecross/out_sigma005", covars = c("age","cross"))

## save summary stats for later
x <- subset(rez$predicted, .iter > 5e3) %>%
	group_by(age, cross) %>%
	do(calc_hpdi_from_vector(exp(.$nu)))
saveRDS(x, "gamma_fits_agecross/group_hpds_sigma005.rds")

## some prespecified marginal means
mmeans <- list(
	# by genotype, averaging over age effect
	"cross" = matrix( c(1,0,0,0,0,0,1/2,
						  1,1,0,0,0,0,1/2,
						  1,0,1,0,0,0,1/2,
						  1,0,0,1,0,0,1/2,
						  1,0,0,0,1,0,1/2,
						  1,0,0,0,0,1,1/2),
						byrow = TRUE, nrow = 6,
						dimnames = list(c("FG","FH","GF","GH","HF","HG"), NULL)),
	# by age, averaging over genotypes
	"ages"  =  matrix( c(1,1/6,1/6,1/6,1/6,1/6,0,
						 1,1/6,1/6,1/6,1/6,1/6,1),
						 byrow = TRUE, nrow = 2,
						 dimnames = list(c("E","Y"), NULL)),
	# by genotype x age (ie marginal means for all groups)
	"crosses_age" = matrix( c(1,0,0,0,0,0,0,
						  1,1,0,0,0,0,0,
						  1,0,1,0,0,0,0,
						  1,0,0,1,0,0,0,
						  1,0,0,0,1,0,0,
						  1,0,0,0,0,1,0,
						  1,0,0,0,0,0,1,
						  1,1,0,0,0,0,1,
						  1,0,1,0,0,0,1,
						  1,0,0,1,0,0,1,
						  1,0,0,0,1,0,1,
						  1,0,0,0,0,1,1),
						byrow = TRUE, nrow = 12,
						dimnames = list(c("FG:E","FH:E","GF:E","GH:E","HF:E","HG:E",
										  "FG:Y","FH:Y","GF:Y","GH:Y","HF:Y","HG:Y"), NULL)),
	# contrasts by X chromosome genotype -- note that these need cross-wise marginal means first
	"X_chrom" =  matrix( c(-1/4,-1/4,-1/4,-1/4,+1/2,+1/2,
						   -1/4,-1/4,+1/2,+1/2,-1/4,-1/4,
						   +1/2,+1/2,-1/4,-1/4,-1/4,-1/4),
						 byrow = TRUE, nrow = 3,
						 dimnames = list(c("dom","mus","cas"), NULL)),
	# contasts between reciprocal genotypes -- note that these need cross-wise marginal means first
	"reciprocals" = matrix( c(+1/2,0,-1/2,0,0,0,
							  0,+1/2,0,0,-1/2,0,
							  0,0,0,1/2,0,-1/2),
							byrow = TRUE, nrow = 3,
							dimnames = list(c("FG-GF","FH-HF","GH-HG"), NULL)),
	# contasts between reciprocal genotypes -- note that these need cross-wise marginal means first
	"all_pairwise" = all_pairwise_comparisons_matrix(c("FG","FH","GF","GH","HF","HG")),
	# contrasts old vs young -- note that these need age-wise marginal means first
	"age_contrast" = matrix( c(1,-1),
							 byrow = TRUE, nrow = 1,
							 dimnames = list(c("E-Y"), NULL))
)

subset(rez$predicted, .iter > 5e3) %>%
	ggplot(., aes(x = cross, y = nu, colour = age)) +
	tidybayes::stat_dist_halfeye(position = position_dodge(width = 1)) +
	coord_flip()

## extract coefficients, for estimating marginal means and contrasts
betas <- subset(rez$coefs, .iter > 5e3) %>% select(-c(.iter,.chain)) %>% as.matrix()
coda::as.mcmc(betas) %>% coda::HPDinterval() %>% exp()

## marginal means by cross 
by_cross <- betas %*% t(mmeans$cross) %>%
	exp() %>%
	as_tibble() %>%
	pivot_longer(-c(), names_to = "group", values_to = "nu")

## HPDIs for nu in each cross
by_cross_hpds <- betas %*% t(mmeans$cross) %>%
	exp() %>%
	calc_hpdi(burnin = 0) %>%
	rename(group = term)

## marginal means by age 
by_age <- betas %*% t(mmeans$ages) %>%
	exp() %>%
	as_tibble() %>%
	pivot_longer(-c(), names_to = "group", values_to = "nu")

## HPDIs for nu in each age
by_age_hpds <- betas %*% t(mmeans$ages) %>%
	exp() %>%
	calc_hpdi(burnin = 0) %>%
	rename(group = term)

## combine marginal means into one df
by_both <- bind_rows("by genotype" = by_cross, "by age" = by_age, .id = "term")
by_both_hpds <- bind_rows("by genotype" = by_cross_hpds, "by age" = by_age_hpds, .id = "term")

## for purposes of ordering things properly on the axis ...
levs <- c(names(CROSSES), names(AGES))
lev_labs <- c(CROSSES, AGES)
by_both$group <- factor(by_both$group, levs)
by_both$term <- relevel( factor(by_both$term), "by genotype" )
by_both_hpds$term <- relevel( factor(by_both_hpds$term), "by genotype" )

## show group means by cross and age
subset(by_both) %>%
	ggplot(., aes(x = group, y = nu)) +
	tidybayes::stat_halfeye() +
	geom_hline(yintercept = 11.3, lty = "dashed", colour = "grey70") +
	geom_text(data = by_both_hpds,
			  aes(x = group, y = estimate, label = round(estimate, 1)),
			  position = position_nudge(x = 0.4, y = 0.3), size = 3) +
	scale_y_continuous(expression(nu), trans = "log2", limits = c(5,40)) +
	scale_x_discrete(labels = lev_labs, limits = rev) +
	lemon::facet_rep_grid(term ~ ., scale = "free_y", space = "free_y", repeat.tick.labels = TRUE) +
	#facet_grid(term ~ ., scale = "free_y", space = "free_y") +
	coord_flip() +
	theme_minimal() +
	theme(axis.title.y = element_blank())
ggsave("figures/gamma_mcmc/posteriors_by_cross_or_age.pdf", width = 4, height = 6)

## show group means by age only
subset(by_both) %>%
	subset(term == "by age") %>%
	ggplot(., aes(x = group, y = nu)) +
	tidybayes::stat_halfeye() +
	scale_y_continuous(expression(nu), trans = "log2",  limits = c(5,40)) +
	scale_x_discrete(labels = c("E" = "old male", "Y" = "young male"), limits = rev) +
	coord_flip() +
	theme_minimal() +
	theme(axis.title.y = element_blank())
ggsave("figures/gamma_mcmc/posteriors_by_age.pdf", width = 3, height = 2)

## what about relationship to recombination rate?
maplen <- read_csv("xos_inferred/maplen_by_cross.csv")
colnames(maplen) <- c("group", "cM", "cM_lo", "cM_hi", "cM_se")
left_join(maplen, by_cross_hpds) %>%
	ggplot() +
	geom_linerange(aes(x = cM, ymin = conf_lo, ymax = conf_hi)) +
	geom_segment(aes(x = cM_lo, xend = cM_hi, y = estimate, yend = estimate)) +
	geom_label(aes(x = cM, y = estimate, label = CROSSES[group]), size = 2.5) +
	scale_y_continuous(expression(nu), trans = "log2") +
	scale_x_continuous("\nmap length (cM)")
ggsave("figures/map_length_vs_nu.pdf", width = 4, height = 4)

## now do contrasts by X chromosome genotype 
# by_Xchr <- (betas %*% t(mmeans$crosses)) %*% t(mmeans$X_chrom) %>%
# 	#exp() %>%
# 	as_tibble() %>%
# 	pivot_longer(-c(), names_to = "group", values_to = "nu")
# 
# ## axis labels for X chromosome contrasts
# x_contrasts <- c("dom" = expression("X"^italic("dom")-"others"),
# 				 "mus" = expression("X"^italic("mus")-"others"),
# 				 "cas" = expression("X"^italic("cas")-"others"))
# 
# ## show contrasts by X chromosome genotype
# p1 <- subset(by_Xchr) %>%
# 	mutate(term = "by X chromosome") %>%
# 	mutate(group = factor_mus(group)) %>%
# 	ggplot(., aes(x = group, y = nu, colour = group)) +
# 	tidybayes::stat_halfeye() +
# 	geom_hline(yintercept = 0, lty = "dashed", colour = "grey50") +
# 	scale_x_discrete(labels = x_contrasts, limits = rev) +
# 	scale_y_continuous(expression(log(nu[i] - bar(nu))), limits = c(-1,1)) +
# 	scale_color_mus(guide = "none") +
# 	#facet_grid(term ~ ., scale = "free_y", space = "free_y") +
# 	coord_flip() +
# 	theme_minimal() +
# 	theme(axis.title.y = element_blank())
# print(p1)
# ggsave("figures/gamma_mcmc/posterior_contrasts_by_X_chrom.pdf", width = 3, height = 3)

## now do all pairwise contrasts between crosses
make_contrast_labels <- function(x) {
	s <- c("F" = "CAST", "G" = "PWK", "H" = "WSB")
	m <- stringr::str_match(x, "(\\w)(\\w)-(\\w)(\\w)")
	m <- m[ ,-1, drop=FALSE ]
	m <- apply(m, 2, function(f) s[f])
	xx <- paste0(m[,1],"x",m[,2], " / ", m[,3],"x",m[,4])
	return(xx)
}
by_pairwise <- (betas %*% t(mmeans$cross)) %*% t(mmeans$all_pairwise) %>%
	as_tibble() %>%
	pivot_longer(-c(), names_to = "group", values_to = "nu")
by_pairwise_hpds <- (betas %*% t(mmeans$cross)) %*% t(mmeans$all_pairwise) %>%
	calc_hpdi(burnin = 0) %>%
	rename(group = term)

by_pairwise_hpds %>%
	mutate(signif = (conf_lo > 0 | conf_hi < 0)) %>%
	ggplot(., aes(x = reorder(group, estimate), shape = signif)) +
	geom_hline(yintercept = 1, lty = "dashed", colour = "red") +
	geom_pointrange(aes(y = exp(estimate), ymin = exp(conf_lo), ymax = exp(conf_hi)), fill = "white") +
	scale_shape_manual(values = c(19,21), guide = "none") +
	scale_y_continuous(expression(nu[i]/nu[j]), trans = "log2") +
	scale_x_discrete("\ncomparison", labels = make_contrast_labels) +
	theme(axis.text.x = element_text(angle = 90, hjust = 1))

## fancier version for paper
by_pairwise %>%
	left_join( mutate(by_pairwise_hpds,	signif = (conf_lo > 0 | conf_hi < 0)) ) %>%
	ggplot(., aes(x = reorder(group, nu))) +
	geom_hline(yintercept = 1, lty = "dashed", colour = "grey50") +
	tidybayes::stat_pointinterval(aes(y = exp(nu), shape = signif), fill = "white") +
	scale_shape_manual(values = c(19,21), guide = "none") +
	scale_y_continuous(expression(nu[i]/nu[j]), trans = "log2") +
	scale_x_discrete("", labels = make_contrast_labels) +
	coord_flip()
ggsave("figures/gamma_mcmc/posterior_contrasts_all_pairwise.pdf", width = 4, height = 4)

## now do contrasts between reciprocal crosses
# axis labels for reciprocal crosses
recip_contrasts <- c("FG-GF" = "CASTxPWK / PWKxCAST",
					 "FH-HF" = "CASTxWSB / WSBxCAST",
					 "GH-HG" = "PWKxWSB / WSBxPWK")
by_reciprocals <- (betas %*% t(mmeans$cross)) %*% t(mmeans$reciprocals) %>%
	#exp() %>%
	as_tibble() %>%
	bind_cols(rez$coefs[ ,c(".iter",".chain") ]) %>%
	pivot_longer(-c(.iter, .chain), names_to = "group", values_to = "nu")
by_reciprocals_hpds <- (betas %*% t(mmeans$cross)) %*% t(mmeans$reciprocals) %>%
	#exp() %>%
	calc_hpdi(burnin = 5e3) %>%
	rename(group = term)

by_reciprocals <- by_reciprocals %>%
	mutate(group = factor(group, names(recip_contrasts)))
by_reciprocals_hpds <- by_reciprocals_hpds %>%
	mutate(group = factor(group, names(recip_contrasts)))

by_reciprocals_hpds %>% mutate(across(-group, exp))

## show contrasts between reciprocal crosses
p2 <- subset(by_reciprocals, .iter > 5e3) %>%
	mutate(term = "by reciprocal crosses") %>%
	ggplot(., aes(x = group, y = exp(nu))) +
	tidybayes::stat_halfeye() +
	geom_hline(yintercept = 1, lty = "dashed", colour = "grey50") +
	geom_text(data = by_reciprocals_hpds,
			  aes(x = group, y = exp(estimate), label = round(exp(estimate), 2)),
			  position = position_nudge(x = 0.4, y = 0.3), size = 3) +
	scale_x_discrete(labels = recip_contrasts, limits = rev) +
	scale_y_continuous(expression(nu[i]/nu[j]), limits = c(0.5,2), trans = "log2") +
	coord_flip() +
	theme_minimal() +
	theme(axis.title.y = element_blank())
print(p2)
ggsave("figures/gamma_mcmc/posterior_contrasts_by_reciprocal_crosses.pdf", width = 4, height = 3)

cowplot::plot_grid(p1, p2, nrow = 2, align = "v")
ggsave("figures/gamma_mcmc/posterior_contrasts_by_X_chrom_and_reciprocals.pdf", width = 4, height = 4)

## now do contrasts old vs young
by_age_contrast <- (betas %*% t(mmeans$ages)) %*% t(mmeans$age_contrast) %>%
	#exp() %>%
	as_tibble() %>%
	bind_cols(rez$coefs[ ,c(".iter",".chain") ]) %>%
	pivot_longer(-c(.iter, .chain), names_to = "group", values_to = "nu")

## axis labels for age contrast
age_contrasts <- c("E-Y" = "old male - young male")

## show contrasts between age groups
p3 <- subset(by_age_contrast, .iter > 5e3) %>%
	mutate(term = "by age contrast") %>%
	mutate(group = factor(group, names(age_contrasts))) %>%
	ggplot(., aes(x = group, y = nu)) +
	tidybayes::stat_halfeye() +
	geom_hline(yintercept = 0, lty = "dashed", colour = "grey50") +
	scale_x_discrete(labels = age_contrasts, limits = rev) +
	scale_y_continuous(expression(log(nu[i] - bar(nu))), limits = c(-1,1)) +
	coord_flip() +
	theme_minimal() +
	theme(axis.title.y = element_blank())
print(p3)
ggsave("figures/gamma_mcmc/posterior_contrasts_by_age.pdf", width = 3, height = 1.5)

cowplot::plot_grid(p1 + theme(axis.title.x = element_blank()),
				   p2 + theme(axis.title.x = element_blank()),
				   p3,
				   nrow = 3, rel_heights = c(3,3,2), align = "v")
ggsave("figures/gamma_mcmc/posterior_contrasts_omnibus.pdf", width = 4, height = 5)





### === old code, no longer used ###
load_mcmc_result <- function(ff, params = c("nu"), burnin = 5000, ...) {
	
	x <- readRDS(ff)
	theta <- as.matrix(x$trace)
	if (!is.null(params))
		colnames(theta) <- params
	rez <- coda::as.mcmc(theta)
	
	return(rez)
	
}

crosses <- tibble(cross = c("FG","GF","FH","HF","GH","HG")) %>%
	mutate(path = paste0("gamma_fits/chains_", cross, ".rds"))

chains <- lapply(crosses$path, load_mcmc_result) %>%
	lapply(as_tibble) %>%
	setNames(., crosses$cross) %>%
	bind_rows(., .id = "cross") %>%
	group_by(cross) %>%
	mutate(.iter = 1:n())

subset(chains, .iter > 1000) %>%
	ggplot(., aes(x = cross, y = nu)) +
	tidybayes::stat_dist_halfeye() +
	scale_y_continuous(trans = "log2") +
	coord_flip()

load_multivar_mcmc <- function(ff, ...) {
	
	x <- readRDS(ff)
	chains <- as.matrix(x$trace)
	colnames(chains) <- x$data$coefs
	chains <- coda::as.mcmc(chains)
	
}

## convert intercept+diffs to group means
make_groupwise_estimates <- function(x, ...) {
	
	thecols <- colnames(x)[-1]
	thecols <- gsub("^cross", "",thecols, perl = TRUE)
	base_group <- setdiff(names(CROSSES), thecols)
	new_cols <-  c(base_group, thecols)
	
	k <- ncol(x)
	beta <- diag(k)
	beta[,1] <- 1
	rez <- x %*% t(beta)
	colnames(rez) <- new_cols
	
	coda::as.mcmc(rez)
	
}

make_pairwise_contrasts <- function(x, exponentiate = FALSE, ...) {
	
	if (exponentiate)
		x <- exp(x)
	
	k <- ncol(x)
	pp <- combn(k, 2)
	beta <- matrix(0, nrow = ncol(pp), ncol = k)
	nn <- paste0(colnames(x)[ pp[1,] ], " - ", colnames(x)[ pp[2,] ])
	rownames(beta) <- nn
	for (ii in seq_len(ncol(pp))) {
		beta[ ii,pp[1,ii] ] <- 1
		beta[ ii,pp[2,ii] ] <- -1
	}
	
	rez <- x %*% t(beta)
	coda::as.mcmc(rez)
	
}

apply_contrast <- function(x, contr_matrix, contr_names = NULL, ...) {
	
	if (!is.matrix(contr_matrix))
		contr_matrix <- matrix(contr_matrix, nrow = 1)
		
	rez <- x %*% t(contr_matrix)
	if (!is.null(contr_names)) {
		if (length(contr_names) == ncol(rez))
			colnames(rez) <- contr_names
	}
	
	coda::as.mcmc(rez)
	
}
	
	
x <- lapply(LETTERS[1:4], function(f) load_multivar_mcmc(paste0("gamma_fits_cross/out_sigma005_", f, ".rds")))
x <- lapply(x, make_groupwise_estimates)
bayesplot::mcmc_combo(x)

chains <- lapply(x, as_tibble) %>%
	bind_rows(.id = "chain") %>% 
	pivot_longer(-chain, names_to = "var", values_to = "value") %>% 
	group_by(chain) %>% 
	mutate(.iter = 1:n(),
		   var = factor(var, names(CROSSES)))

subset(chains, var != "(Intercept)" & .iter > 5e3) %>%
	ggplot(aes(x = var, y = exp(value))) + 
	tidybayes::stat_dist_halfeye() + 
	scale_y_continuous(trans = "log2") + 
	scale_x_crosses() +
	coord_flip()

pw_chains <- lapply(x, make_pairwise_contrasts) %>%
	lapply(., as_tibble) %>%
	bind_rows(.id = "chain") %>% 
	pivot_longer(-chain, names_to = "contrast", values_to = "value") %>% 
	group_by(chain) %>% 
	mutate(.iter = 1:n())

subset(pw_chains, .iter > 5e3) %>%
	ggplot(aes(x = contrast, y = (value))) + 
	tidybayes::stat_dist_halfeye() + 
	geom_hline(yintercept = 0, lty = "dashed", colour = "darkred") +
	coord_flip()

contr_matrix <- matrix(
	c(-1/4,-1/4,-1/4,-1/4,1/2,1/2,
	  -1/4,-1/4,1/2,1/2,-1/4,-1/4,
	  1/2,1/2,-1/4,-1/4,-1/4,-1/4),
	byrow = TRUE, nrow = 3 )

mus_X_chains <- lapply(x, apply_contrast,
					   contr_matrix = contr_matrix,
					   contr_names = c("dom_X","mus_X","cas_X")) %>%
	lapply(., as_tibble) %>%
	bind_rows(.id = "chain") %>% 
	pivot_longer(-chain, names_to = "contrast", values_to = "value") %>% 
	group_by(chain) %>% 
	mutate(.iter = 1:n())

subset(mus_X_chains, .iter > 5e3) %>%
	ggplot(aes(x = contrast, y = (value))) + 
	tidybayes::stat_dist_halfeye() + 
	geom_hline(yintercept = 0, lty = "dashed", colour = "darkred") +
	coord_flip()



### ======== ###
## posterior predictive sims
sims <- group_by(rez$predicted, cross) %>%
	slice_sample(n = 1000) %>%
	mutate(draw = 1:n()) %>%
	group_by(cross, draw) %>%
	do(sim_gamma_model(nu = exp(.$nu), nchrom = 20))

sims2 <- group_by(rez$predicted, cross) %>%
	subset(.iter > 5000) %>%
	summarise(nu = exp(mean(nu))) %>%
	group_by(cross) %>%
	do(sim_gamma_model(nu = .$nu, nchrom = 1000, obligate_chiasma = TRUE, verbose = TRUE))

mirror_sims <- bind_rows(
	select(xoloc_all, cross, d, censor) %>% mutate(src = "observed"),
	select(sims2, cross, d, censor) %>% mutate(src = "simulated") )

subset(mirror_sims, censor == 0) %>%
	ggplot(., aes(x = cross, y = d, color = src, fill = src)) +
	tidybayes::stat_dist_halfeye(position = position_dodge(width = 0.3), slab_alpha = 0.2) +
	theme_bw()
	
