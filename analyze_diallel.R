setwd("~/Dropbox/pmdvlab/xoi/")

library(tidyverse)
library(BayesDiallel)

devtools::load_all("~/Dropbox/pmdvlab/mouser/")
devtools::load_all("~/Dropbox/pmdvlab/popcorn/")

source("helpers/bayesdiallal_helpers.R")
old_theme <- theme_set(theme_minimal())

nxo_all <- readRDS("xos_inferred/xo_counts.rds")

obj <- DiallelAnalyzer(data = mutate(nxo_all, is_female = (sex == 2), age = (age == "E")),
					   mother.strain = "mat", father.strain = "pat",
					   Models = "amv", burnin = 1000, lengthChains = 10000, 
					   #Models = "Babmdw",
					   phenotype = "n_xo", is.female = "is_female",
					   LogTransform = FALSE,
					   FixedEffects = "age")

rez1 <- diallel_hpdi(obj, prob = c(0.5, 0.9))
pred1 <- predict_posteriors(obj)

saveRDS(obj, "xos_inferred/bayesdiallel_result.rds")

plot_diallel(obj, label = TRUE) + 
	#scale_fill_distiller("mean COs", palette = "RdBu", na.value = NA) +
	viridis::scale_fill_viridis("mean COs") +
	coord_equal() +
	theme(axis.line = element_blank()) +
	ggtitle("COs per meiosis")

## show effect sizes for additive, maternal and fixed (=age) effect terms
tmp <- tidy_chains(obj, strains = c("F","G","H")) %>% 
	subset(group %in% c("additive","maternal","fixed effects")) %>% 
	mutate(strain = case_when(original_col == "FixedEffect:1" ~ "EY",
							  TRUE ~ strain),
		   group = case_when(original_col == "FixedEffect:1" ~ "age",
		   				  TRUE ~ group))
tmp_labs <- c("H" = "WSB", "G" = "PWK", "F" = "CAST", "EY" = "old vs young")
tmp_cols <- c( setNames(unname(mus_colors()[ c("dom","mus","cas") ]), c("H","G","F")),
			   "EY" = "black")
tmp$group <- factor(tmp$group, levels = c("additive","maternal","age"))
tmp$strain <- factor(tmp$strain, names(tmp_labs))
		   
p1 <- tmp %>%
	ggplot() + 
	tidybayes::stat_dist_halfeye(aes(x = strain, y = value, colour = strain)) +
	geom_hline(aes(yintercept = 0), lty = "dashed", colour = "grey70") +
	facet_grid(. ~ group, scales = "free_x", space = "free_x") +
	scale_y_continuous("effect size (COs per meiosis)\n") +
	scale_x_discrete(labels = tmp_labs) +
	scale_colour_manual(values = tmp_cols, guide = "none") +
	coord_cartesian(ylim = c(-3,3)) +
	theme_slanty_x()

## show heritability decompsition (looks weird)
p2 <- get_hsq(obj, use_fixed_effects = FALSE) %>% 
	tidy_hsq(flip_parents = TRUE) %>%
	ggplot() +
	tidybayes::stat_dist_halfeye(aes(x = term, y = value)) +
	scale_y_continuous("\nproportion of variance", trans = "sqrt") +
	coord_flip() +
	theme(axis.title.y = element_blank())

## combo plot for manuscript
cowplot::plot_grid(p1, p2, nrow = 1, rel_widths = c(1, 1),
				   labels = "AUTO")
ggsave("figures/diallel_result_summary.pdf", width = 8, height = 4)

## print summary stats from model output, without making column names pretty
get_chains(obj) %>% summary()
get_hsq(obj, use_fixed_effects = TRUE) %>% summary()

## posterior contrasts
posterior_effect_contrast(obj, group = "additive", focal = c("F"), others = "all") %>% summary()
