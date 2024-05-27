setwd("~/Dropbox/pmdvlab/xoi/")

library(tidyverse)
library(ggplot2)

devtools::load_all("~/Dropbox/pmdvlab/mouser/")
devtools::load_all("~/Dropbox/pmdvlab/popcorn/")

source("helpers/xo_helpers.R")

## set global ggplot2 theme
old_theme <- theme_set(theme_minimal())

scale_color_interference_models <- function(...) {
	MetBrewer::scale_color_met_d("Isfahan1", ...)
}

calc_model_aic <- function(x, ...) {
	
	k_gamma <- nrow(x)
	k_stahl <- 2*nrow(x)
	ll_gamma <- sum(x$loglik_0)
	ll_stahl <- sum(x$loglik)
	
	aic_gamma <- 2*k_gamma - 2*ll_gamma
	aic_stahl <- 2*k_stahl - 2*ll_stahl
	return( tibble(k_gamma, k_stahl, ll_gamma, ll_stahl, aic_gamma, aic_stahl) )
	
}

## read parameter estimates from bootstrap replicates (computed on the cluster)
tmp <- readRDS("stahl_fits_cross/bootstraps_all.rds")
tmp$cross <- factor(tmp$cross, levels = names(CROSSES))
point_est <- tmp %>%
	subset(rep_local == 0) %>%
	group_by(cross) %>%
	#group_by(age) %>%
	slice(1)

## read parameter estimates from null model (ie one set of interference params fit to whole dataset)
x <- readRDS("stahl_fits_null.rds")
list(full = point_est, null = x) %>%
	lapply(calc_model_aic) %>%
	bind_rows(.id = "model")

## calculate all pairwise differences in specified columns "which_columns"
## pairs are defined as distinct rows, presumed indexed by a variable in "group"
pairwise_rowwise_diffs <- function(reps, group, which_cols = c("p","nu","nu_0"), ...) {
	
	.calc_diffs <- function(rows) {
		diff(d[rows,])
	}
	
	nr <- nrow(reps)
	d <- as.matrix(reps[ ,which_cols ])
	pp <- t(combn(nr, 2))
	pw <- t( apply(pp, 1, .calc_diffs) )
	colnames(pw) <- which_cols
	nn <- tibble(x1 = group[ pp[,1] ], x2 = group[ pp[,2] ])
	rez <- as_tibble(pw)
	bind_cols(nn, rez)
	
}

## not used
bootstrap_ci <- function(reps, alpha = 0.05, ...) {
	theta <- reps$nu[ which(reps$rep == 0)[1] ]
	theta_star <- reps$nu[ reps$rep > 0 ]
	deltas <- theta - theta_star
	tibble( estimate = theta,
			lo = theta + quantile(deltas, alpha/2),
			hi = theta + quantile(deltas, 1-alpha/2))
}

## calculate non-parametric bootstrap CIs for parameters of Stahl (gamma-escape) model
make_quantile_intervals_for_params <- function(reps, alpha = 0.05, ...) {
	
	rez <- summarise(reps,
					 p_lo = quantile(p, alpha/2),
					 p_lomid = quantile(p, 1/4),
					 p_himid = quantile(p, 3/4),
					 p_hi = quantile(p, 1-alpha/2),
					 nu_lo = quantile(nu, alpha/2),
					 nu_lomid = quantile(nu, 1/4),
					 nu_himid = quantile(nu, 3/4),
					 nu_hi = quantile(nu, 1-alpha/2),
	)
	rez$model <- "gamma-escape"
	rez_null <- summarise(reps,
						  p_lo = NA_real_,
						  p_lomid = NA_real_,
						  p_himid = NA_real_,
						  p_hi = NA_real_,
						  nu_lo = quantile(nu_0, alpha/2),
						  nu_lomid = quantile(nu_0, 1/4),
						  nu_himid = quantile(nu_0, 3/4),
						  nu_hi = quantile(nu_0, 1-alpha/2))
	rez_null$model <- "gamma"
	bind_rows(rez, rez_null)
	
}

## calculate non-parametric bootstrap CIs for pairwise differences in parameters of Stahl (gamma-escape) model
make_quantile_intervals_for_diffs <- function(reps, alpha = 0.05, ...) {
	
	rez <- summarise(reps,
					 p_med = quantile(p, 0.5),
					 p_lo = quantile(p, alpha/2),
					 p_lomid = quantile(p, 1/4),
					 p_himid = quantile(p, 3/4),
					 p_hi = quantile(p, 1-alpha/2),
					 nu_med = quantile(nu, 0.5),
					 nu_lo = quantile(nu, alpha/2),
					 nu_lomid = quantile(nu, 1/4),
					 nu_himid = quantile(nu, 3/4),
					 nu_hi = quantile(nu, 1-alpha/2),
					 nu0_med = quantile(nu_0, 0.5),
					 nu0_lo = quantile(nu_0, alpha/2),
					 nu0_lomid = quantile(nu_0, 1/4),
					 nu0_himid = quantile(nu_0, 3/4),
					 nu0_hi = quantile(nu_0, 1-alpha/2)
	)
	
	rez <- mutate(rez,
				  p_signif = (p_lo > 0 | p_hi < 0),
				  nu_signif = (nu_lo > 0 | nu_hi < 0),
				  nu0_signif = (nu0_lo > 0 | nu0_hi < 0))
	return(rez)
}

## calculate bootstrap non-parametric CIs, join with point estimates
tmp2 <- tmp %>%
	subset(rep_local > 0) %>%
	group_by(cross) %>%
	#group_by(age) %>%
	do(make_quantile_intervals_for_params(.)) %>%
	left_join(point_est) %>%
	mutate(nu_hat = case_when(model == "gamma" ~ nu_0,
							  model == "gamma-escape" ~ nu),
		   p_hat = case_when(model == "gamma" ~ NA_real_,
		   				  model == "gamma-escape" ~ p))

## compare to MCMC results
x <- readRDS("gamma_fits_agecross/group_hpds_sigma005.rds")
subset(x, age == "E") %>% 
	left_join(subset(tmp2,model == "gamma")) %>% 
	mutate(cross = factor_cross(cross, long_names = TRUE)) %>%
	ggplot(aes(x = estimate, y = nu_0)) + 
	geom_abline(intercept = 0, slope = 1, lty = "dashed", colour = "grey70") +
	geom_smooth(method = "lm", fullrange = TRUE, colour = "grey50") +
	geom_segment(aes(x = conf_lo, xend = conf_hi, y = nu_0, yend = nu_0)) +
	geom_segment(aes(y = nu_lo, yend = nu_hi, x = estimate, xend = estimate)) +
	geom_label(aes(label = cross), size = 2.5) + 
	scale_x_continuous("MCMC estimate (with 95% HPDI)", trans = "log2") +
	scale_y_continuous("ML estimate (with 95% bootstrap CI)", trans = "log2") +
	coord_equal()
ggsave("figures/gamma_mcmc/comparison_to_bootstrap_method.pdf", width = 4, height = 4)

## RMSD between the two methods
subset(x, age == "E") %>% 
	left_join(subset(tmp2,model == "gamma")) %>% 
	group_by("all") %>%
	summarise(rmsd = sum( ((nu_0-estimate)^2)/length(nu_0) ),
			  rmsd_rel = rmsd/min(estimate))

## calculate bootstrap estimates of all pairwise differences
tmp_reps <- subset(tmp, rep_local > 0)
tmp_reps$rep <- droplevels(tmp_reps$rep)
nlevels(tmp_reps$rep)
tmp_reps <- tmp_reps %>%
	group_by(cross) %>%
	mutate(rep = seq_along(rep))
pw_boot <- subset(tmp_reps) %>% 
	group_by(rep) %>%
	do(pairwise_rowwise_diffs(., .$cross))
pw_boot_cis <- pw_boot %>%
	group_by(x1, x2) %>%
	do(make_quantile_intervals_for_diffs(.))

## plot bootstrap estimates of pairwise contrasts
ggplot(pw_boot) +
	ggbeeswarm::geom_quasirandom(aes(x = interaction(x2, x1, sep = " - "), y = nu)) +
	geom_hline(yintercept = 0, lty = "dashed", col = "red") +
	theme_slanty_x()

## plot boostrap CIs of pairwise contrasts
ggplot(pw_boot_cis, aes(x = reorder(interaction(x2, x1, sep = " - "), nu0_med))) +
	geom_linerange(aes(ymin = nu0_lo, ymax = nu0_hi)) +
	geom_linerange(aes(ymin = nu0_lomid, ymax = nu0_himid), lwd = 1.5) +
	geom_point(aes(y = nu0_med), pch = 21, size = 2, fill = "white") +
	geom_hline(yintercept = 0, lty = "dashed", colour = "red") +
	theme_slanty_x()

## combined plot of strength of interference (nu) and escape fraction (p)
p1 <- ggplot(tmp2, aes(x = cross, colour = model)) +
	geom_linerange(aes(ymin = nu_lo, ymax = nu_hi),
				   position = position_dodge(width = 1)) +
	geom_linerange(aes(ymin = nu_lomid, ymax = nu_himid), lwd = 1.5,
				   position = position_dodge(width = 1)) +
	geom_point(aes(y = nu_hat), pch = 21, size = 2, fill = "white",
			   position = position_dodge(width = 1)) +
	scale_color_interference_models(guide = "none") +
	scale_x_crosses() +
	scale_y_continuous(expression(hat(nu)~"with 95% CI"), trans = "log2") +
	theme_slanty_x()
p3 <- ggplot(tmp2, aes(x = cross, colour = model)) +
	geom_linerange(aes(ymin = p_lo, ymax = p_hi),
				   position = position_dodge(width = 1)) +
	geom_linerange(aes(ymin = p_lomid, ymax = p_himid), lwd = 1.5,
				   position = position_dodge(width = 1)) +
	geom_point(aes(y = p_hat), pch = 21, size = 2, fill = "white",
			   position = position_dodge(width = 1)) +
	scale_x_crosses() +
	scale_y_continuous(expression(hat(italic(p))~"with 95% CI")) +
	scale_color_interference_models() +
	theme_slanty_x() +
	theme(legend.position = "bottom",
		  legend.title = element_blank())

## plot nu and p against each other, with bootstrap CIs
## (first get raw correlation estimate)
rho <- subset(tmp2, model == "gamma-escape") %>% with(., cor(p, nu, method = "spearman"))
rho_label <- sprintf("atop(\"Spearman's\",rho == %0.2f)", rho)
p4 <- subset(tmp2, model == "gamma-escape") %>%
	ggplot(., aes(group = cross, colour = model)) +
	geom_segment(aes(x = nu_lo, xend = nu_hi, y = p, yend = p)) +
	geom_segment(aes(x = nu_lomid, xend = nu_himid, y = p, yend = p), lwd = 2) +
	geom_segment(aes(y = p_lo, yend = p_hi, x = nu, xend = nu)) +
	geom_segment(aes(y = p_lomid, yend = p_himid, x = nu, xend = nu), lwd = 2) +
	geom_label(aes(x = nu, y = p, label = CROSSES[as.character(cross)]), size = 2) +
	annotate("text", y = 0.23, x = 64, label = rho_label, parse = TRUE) +
	scale_x_continuous(expression(hat(nu)), trans = "log2") +
	scale_y_continuous(expression(hat(italic(p)))) +
	scale_color_interference_models(guide = "none") +
	theme(plot.margin = unit(c(6,3,6,0), "lines"))
print(p4)

## final combo plots
pp1 <- cowplot::plot_grid(p1 + theme(axis.text.x = element_blank()), p3, nrow = 2, align = "v", labels = "AUTO")
print(pp1)
ggsave("figures/stahl_fits_nu_p_bycross.pdf", width = 4, height = 4)

cowplot::plot_grid(pp1, p4, labels = c("", "C"))
ggsave("figures/stahl_fits_combo.pdf", width = 8, height = 6)

tmp <- mutate(tmp, cross = factor_cross(as.character(cross), long_names = TRUE)) 
tmp_rho <- tmp %>%
	group_by(cross) %>%
	summarise(rho = cor(nu, p, method = "spearman"))
tmp %>% 
	ggplot() + 
	geom_point(aes(x = nu, y = p), colour = "darkblue", alpha = 0.25) + 
	geom_text(data = tmp_rho,
			  aes(x = 64, y = 0.05, label = sprintf("rho == %0.02f", rho)), parse = TRUE) +
	scale_x_continuous(expression("strength of interference"~(nu)), trans = "log2") + 
	scale_y_continuous(expression("rate of interference escape"~(italic(p))), trans = "sqrt") + 
	facet_wrap(~ cross, ncol = 2)
ggsave("figures/stahl_model_param_correlations.pdf", width = 6, height = 6)
