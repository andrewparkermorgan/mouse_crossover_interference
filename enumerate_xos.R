## infer crossover events from genotypes, merging across platforms
setwd("~/Dropbox/pmdvlab/xoi")

library(tidyverse)

devtools::load_all("~/Dropbox/pmdvlab/mouser/")
devtools::load_all("~/Dropbox/pmdvlab/popcorn/")
source("helpers/geno_helpers.R")
source("helpers/xo_helpers.R")

## set global ggplot2 theme
old_theme <- theme_set(theme_minimal())

mm_mini <- read_map("array_annot/mini_mm39.map")
mm_mega <- read_map("array_annot/mega_mm39.map")

# bootstrap_map_length <- function(x, n = 100L, hmm_error_prob = 0.01, ...) {
# 	
# 	## run HMM to get genotype probs, then get Viterbi solution and use that for analysis
# 	message("Estimating genotype probabilities ... ")
# 	geno_smooth <- qtl2::viterbi(x, error_prob = hmm_error_prob, lowmem = TRUE, quiet = TRUE)
# 	
# 	## get crossover count (ind x chr matrix)
# 	xos <- qtl2::count_xo(geno_smooth)
# 	n_meioses <- qtl2::n_ind(x)
# 	
# 	## do bootstrapping
# 	rez <- vector("numeric", n)
# 	for (ii in seq(0, n)) {
# 		
# 		keep_iids <- sample(1:nrow(xos), replace = TRUE)
# 		n_xo <- sum(xos[ keep_iids, ])
# 		rez[ii] <- n_xo
# 		
# 	}
# 	
# 	return( tibble(rep = 1:n, n_xo = rez, n_meioses = n_meioses, cM = 100*n_xo/n_meioses) )
# 	
# }

## go from genotypes to output of xoi::convertxoloc
genotypes_to_xoloc <- function(x, hmm_error_prob = 0.01, return_df = TRUE, ...) {

	message("Estimating genotype probabilities ... ")
	breaks <- find_breaks(x, x$gmap, error_prob = hmm_error_prob, ...)
	message("Enumerating crossover locations ... ")
	
	if (return_df) {
		return( convert_xo_loc(breaks) )
	}
	else {
		return(breaks)
	}
		
}

summarise_cross <- function(x, annot = NULL, ...) {
	
	n_ind <- qtl2::n_ind(x)
	this_map <- map_to_df(x, annot) %>% subset(in_map)
	rez <- group_by(this_map, chr) %>%
		summarise(n_markers = n(),
				  pos_min = min(pos),
				  pos_max = max(pos),
				  pos_span = (pos_max - pos_min)/1e6,
				  cM_min = min(cM_for_qtl),
				  cM_max = max(cM_for_qtl),
				  cM_span = cM_max - cM_min) %>%
		mutate(n_ind = n_ind)
	return(rez)
	
}

## load cleaned genotypes
crosses_mini <- readRDS("geno/crosses_mini_clean.rqtl.rds")
crosses_mega <- readRDS("geno/crosses_mega_clean.rqtl.rds")

define_chr_bounds <- function(df, ...) {

	.def_bounds <- function(x) {
		c( max(x$cM_min), min(x$cM_max) )
	}
		
	split(df, df$chr) %>%
	lapply(., .def_bounds)
	
}

## summary of genome coverage
summ_mini <- lapply(crosses_mini, summarise_cross, mm_mini) %>% bind_rows(.id = "cross")
summ_mega <- lapply(crosses_mega, summarise_cross, mm_mega) %>% bind_rows(.id = "cross")
summ_all <- bind_rows( list("mini" = summ_mini, "mega" = summ_mega),
					   .id = "platform" )
summ_all$chr <- factor_chrom( as.character(summ_all$chr), short = TRUE )
## these are chromosome bounds to be used for clipping below
inner_chr_bounds <- split(summ_all, summ_all$cross) %>%
	lapply(., define_chr_bounds)

## Next, need to reduce maps from Mega to Mini
## First clip off left and right ends of chromosome to match array with shortest span
clip_chrom_ends <- function(x, bounds, ...) {
	
	gmap <- x$gmap
	new_map <- vector("list", length(gmap))
	names(new_map) <- names(gmap)
	for (cc in names(gmap)) {
		keeps <- gmap[[cc]] >= bounds[[cc]][1] & gmap[[cc]] <= bounds[[cc]][2]
		new_map[[cc]] <- gmap[[cc]][ unname(keeps) ]
	}
	keeps <- Reduce(c, new_map) %>% names()
	qtl2::pull_markers(x, keeps)
	
}

crosses_mini_clipped <- vector("list", length(crosses_mini))
names(crosses_mini_clipped) <- names(crosses_mini)
for (cc in names(crosses_mini)) {
	crosses_mini_clipped[[cc]] <- clip_chrom_ends(crosses_mini[[cc]], inner_chr_bounds[[cc]]) 
}

crosses_mega_clipped <- vector("list", length(crosses_mega))
names(crosses_mega_clipped) <- names(crosses_mega)
for (cc in names(crosses_mega)) {
	crosses_mega_clipped[[cc]] <- clip_chrom_ends(crosses_mega[[cc]], inner_chr_bounds[[cc]])
	if (cc %in% names(crosses_mini_clipped)) {
		crosses_mega_clipped[[cc]] <- prune_to_nearest_marker(crosses_mega_clipped[[cc]],
															  crosses_mini_clipped[[cc]]$gmap)
	}
}

## how many markers are left?
sapply(crosses_mega_clipped, qtl2::tot_mar)
sapply(crosses_mini_clipped, qtl2::tot_mar)

## summary of genome coverage again, this time after map trimming
summ_mini <- lapply(crosses_mini_clipped, summarise_cross, mm_mini) %>% bind_rows(.id = "cross")
summ_mega <- lapply(crosses_mega_clipped, summarise_cross, mm_mega) %>% bind_rows(.id = "cross")
summ_all <- bind_rows( list("mini" = summ_mini, "mega" = summ_mega),
					   .id = "platform" )
saveRDS(summ_all, "xos_inferred/map_summaries.rds")

## calculate missingness by individual
miss_mini <- lapply(crosses_mini_clipped, missing_by_ind) %>% 
	lapply(enframe) %>% 
	bind_rows() %>%
	rename(iid = name, p_miss = value)
miss_mega <- lapply(crosses_mega_clipped, missing_by_ind) %>% 
	lapply(enframe) %>% 
	bind_rows() %>%
	rename(iid = name, p_miss = value)
miss_all <- bind_rows(miss_mini, miss_mega)
# overall genotyping rate
1-mean(miss_all$p_miss)

## calculate missingness by marker
miss_mini <- lapply(crosses_mini_clipped, missing_by_ind) %>% 
	lapply(enframe) %>% 
	bind_rows() %>%
	rename(iid = name, p_miss = value)
miss_mega <- lapply(crosses_mega_clipped, missing_by_ind) %>% 
	lapply(enframe) %>% 
	bind_rows() %>%
	rename(iid = name, p_miss = value)
miss_all <- bind_rows(miss_mini, miss_mega)

## save the results
saveRDS(crosses_mega_clipped, "geno/crosses_mega_clean_thinned.rqtl.rds")
saveRDS(crosses_mini_clipped, "geno/crosses_mini_clean_thinned.rqtl.rds")

## process genotypes and just count crossovers
nxo_mini <- lapply(crosses_mini_clipped, just_count_crossovers) %>%
	lapply(., rowSums) %>%
	lapply(., enframe, name = "iid", value = "n_xo") %>%
	bind_rows()
nxo_mega <- lapply(crosses_mega_clipped, just_count_crossovers) %>%
	lapply(., rowSums) %>%
	lapply(., enframe, name = "iid", value = "n_xo") %>%
	bind_rows()
nxo_all <- bind_rows( list("mini" = nxo_mini, "mega" = nxo_mega),
					  .id = "platform" )

## re-do crossover analysis with pruned data
xo_mega <- lapply(crosses_mega_clipped, genotypes_to_crossovers)
xo_mega <- bind_rows(xo_mega, .id = "cross")

xo_mini <- lapply(crosses_mini_clipped, genotypes_to_crossovers)
xo_mini <- bind_rows(xo_mini, .id = "cross")

## aggregate crossovers across platforms
xo_all <- bind_rows( list("mini" = xo_mini, "mega" = xo_mega),
					 .id = "platform" )

## NB: some of the samples were initially assigned to wrong (reciprocal) cross due to breeding errors
## master sample sheet has the correct information
samples <- readxl::read_xlsx("~/Dropbox/PAE/sample_info/genotyped.xlsx")
xo_all$cross <- NULL
xo_all <- left_join(xo_all, samples)

nxo_all <- group_by(xo_all, iid, chr) %>%
	summarise(n_xo = n()-1) %>%
	summarise(n_xo = sum(n_xo))
nxo_all <- left_join(nxo_all, samples)
nxo_all <- mutate(nxo_all,
				  mat = substr(cross, 1, 1),
				  pat = substr(cross, 2, 2))

## save results
write_csv(xo_all, "xos_inferred/xos_detail.csv")
saveRDS(xo_all, "xos_inferred/xo_detail.rds")
write_csv(nxo_all, "xos_inferred/xo_counts.csv")

## bootstrap estimates of map length by cross
bootstrap_map_length <- function(df, n = 100L, ...) {
	
	n_meioses <- nrow(df)
	rez <- vector("numeric", n)
	for (ii in seq(1, n)) {
		
		idx <- sample.int(n_meioses, replace = TRUE)
		rez[ii] <- sum(df$n_xo[idx])
		
	}
	
	tibble(n_xo = rez, cM = 100*rez/n_meioses, n_meioses = n_meioses)
	
}

nxo_all$cross <- factor_cross(nxo_all$cross)
nxo_all$age <- factor_age(nxo_all$age)
saveRDS(nxo_all, "xos_inferred/xo_counts.rds")

xo_boot <- group_by(nxo_all, cross, age) %>%
	do(bootstrap_map_length(., n = 1000))
xo_boot_noage <- group_by(nxo_all, cross) %>%
	do(bootstrap_map_length(., n = 1000))
xo_boot_nocross <- group_by(nxo_all, age) %>%
	do(bootstrap_map_length(., n = 1000))

## summarise per-cross estimates, save for later
x <- group_by(xo_boot_noage, cross) %>%
	summarise(mean = mean(cM), conf_lo = quantile(cM, 0.05), conf_hi = quantile(cM, 0.975), se = sd(cM))
write_csv(x, "xos_inferred/maplen_by_cross.csv")

## plot crossovers per meiosis, individual level
nxo_all %>%
	ggplot() +
	tidybayes::stat_dist_dotsinterval(aes(x = cross, y = n_xo), dotsize = 1.2) +
	scale_y_continuous("number of COs\n") +
	scale_x_crosses() +
	theme_slanty_x()

p1 <- nxo_all %>%
	ggplot() +
	tidybayes::stat_dist_dotsinterval(aes(x = cross, y = n_xo, colour = age), dotsize = 1.2,
									  position = position_dodge(width = 0.5)) +
	MetBrewer::scale_color_met_d("Demuth", name = "paternal age", labels = AGES) +
	scale_y_continuous("number of COs\n") +
	scale_x_crosses() +
	theme_slanty_x() +
	theme(legend.position = "top")
print(p1)
ggsave("figures/nxo_by_age_cross.pdf", width = 5, height = 4)

## plot bootstrap map length estimates by cross and age
p2 <- ggplot(xo_boot) +
	tidybayes::stat_halfeye(aes(x = cross, colour = age, y = cM),
							position = position_dodge(width = 0.5)) +
	MetBrewer::scale_color_met_d("Demuth", name = "paternal age", labels = AGES) +
	scale_y_continuous("map length (cM)\n", limits = c(650,1800)) +
	scale_x_crosses() +
	theme_slanty_x()
print(p2)
ggsave("figures/map_length_boostrap_by_cross_and_age.pdf", width = 6, height = 3)

## plot bootstrap map length estimates by cross only
ggplot(xo_boot_noage) +
	tidybayes::stat_halfeye(aes(x = cross, y = cM)) +
	scale_y_continuous("map length (cM)\n", limits = c(650,1800)) +
	scale_x_crosses() +
	theme_slanty_x()
ggsave("figures/map_length_boostrap_by_cross.pdf", width = 4, height = 3)

## plot bootstrap map length estimates by age only
ggplot(xo_boot_nocross) +
	tidybayes::stat_halfeye(aes(x = age, y = cM)) +
	scale_y_continuous("map length (cM)\n") +
	scale_x_discrete(labels = AGES) +
	theme_slanty_x()

## combo plot showing both actual crossover counts and bootstrapped map lengths
cowplot::plot_grid(p1 + theme(axis.text.x = element_blank()), p2 + guides(colour = "none"),
				   rel_heights = c(1.25, 1), nrow = 2, align = "v", labels = "AUTO")
ggsave("figures/map_length_combo_by_age_cross.pdf", width = 6, height = 6)

## how many chromosomes have multiple crossovers?
n_doubles <- xo_all %>%
	mutate(iid = factor(iid),
		   chr = factor_chrom(chr, short = TRUE)) %>%
	group_by(iid, chr, .drop = FALSE) %>%
	summarise(is_double = any(censor == 0)) %>%
	summarise(n_double = sum(is_double)) %>%
	left_join(samples)
n_doubles <- mutate(n_doubles,
					mat = substr(cross, 1, 1),
					pat = substr(cross, 2, 2))

n_doubles %>% 
	ggplot() +
	tidybayes::stat_dist_dotsinterval(aes(x = cross, y = n_double, colour = age), dotsize = 1.2,
									  position = position_dodge(width = 0.5)) +
	MetBrewer::scale_color_met_d("Demuth", name = "paternal age", labels = AGES) +
	scale_y_continuous("chromosomes with multiple COs\n") +
	scale_x_crosses() +
	theme_slanty_x() +
	theme(legend.position = "top")

## show breakdown of chromosomes by number of COs
nxo_by_chr <- xo_all %>%
	mutate(iid = factor(iid),
					 chr = factor_chrom(chr, short = TRUE)) %>%
	group_by(iid, chr) %>%
	summarise(n_xo = n()-1) %>%
	left_join(samples) %>%
	mutate(cross = factor_cross(cross, long_names = TRUE),
		   age = factor_age(age, long_names = TRUE))

p1 <- nxo_by_chr %>%
	ggplot() +
	geom_bar(aes(x = cross:age, fill = factor(n_xo)), position = "fill") +
	scale_y_continuous("proportion of meioses\n") +
	scale_x_discrete("\ngenotype : age combinations") +
	MetBrewer::scale_fill_met_d("Demuth", name = "COs/chromsome") +
	theme(axis.text.x = element_text(angle = 90, hjust = 1))

## show relationship between number of COs and number of multiple-CO chromosomes
m3 <- glm(n_xo ~ cross*age, nxo_all, family = poisson())
m3_mm <- lsmeans::lsmeans(m3, c("cross","age"), type = "response") %>% broom::tidy(conf.int = TRUE)
m4 <- glm(n_double ~ cross*age, n_doubles, family = poisson())
m4_mm <- lsmeans::lsmeans(m4, c("cross","age"), type = "response") %>% broom::tidy(conf.int = TRUE)

tmp <- left_join( 
	rename(m3_mm, n_xo_est = rate, n_xo_lo = conf.low, n_xo_hi = conf.high) %>% 
		select(-c(std.error,df,null,statistic,p.value)),
	rename(m4_mm, n_multi_est = rate, n_multi_lo = conf.low, n_multi_hi = conf.high) %>% 
		select(-c(std.error,df,null,statistic,p.value))
)

p2 <- tmp %>%
	mutate(cross = factor_cross(cross, long_names = TRUE),
		   age = factor_age(age, long_names = TRUE)) %>%
	ggplot(., aes(x = n_xo_est, y = n_multi_est)) +
	geom_smooth(method = "lm", fullrange = FALSE, colour = "grey50") +
	geom_segment(aes(x = n_xo_lo, xend = n_xo_hi, y = n_multi_est, yend = n_multi_est)) +
	geom_segment(aes(x = n_xo_est, xend = n_xo_est, y = n_multi_lo, yend = n_multi_hi)) +
	geom_label(aes(label = paste0(cross, "\n", age)), size = 2) +
	scale_x_continuous("total COs") +
	scale_y_continuous("count of chromosomes with multiple COs")

# combo plot, showing breakdown of chromosomes by count of COs
cowplot::plot_grid(p1, p2, nrow = 1, labels = "AUTO")
ggsave("figures/total_cos_vs_multi_cos.pdf", width = 10, height = 6)



## generate non-parametric summary statistics
alpha <- 0.05
group_by(xo_boot, cross, age) %>%
	summarise(cM_median = median(cM),
			  cM_lo = quantile(cM, alpha/2),
			  cM_hi = quantile(cM, 1-alpha/2))

## summarize and plot number of crossovers per chromosome
xos_per_chr <- subset(xo_all, is_autosome(chr)) %>%
	group_by(platform, cross, iid, chr) %>%
	summarise(nxo = n()-1) %>%
	ungroup() %>% 
	group_by(platform, cross) %>%
	mutate(n_ind = length(unique(iid))) %>%
	group_by(platform, cross, chr, nxo, n_ind) %>%
	summarize(n_class = n())

ggplot(xos_per_chr) +
	geom_bar(aes(x = factor_chrom(chr, short = TRUE), fill = factor(nxo),
				 y = n_class/n_ind), stat = "identity") +
	scale_fill_discrete("COs/chrom") +
	scale_x_discrete("\nchromosome") +
	facet_grid(cross ~ platform) +
	theme_bw()

## summarize and plot number of crossovers per individual
xos_summ <- subset(xo_all, is_autosome(chr)) %>%
	group_by(platform, cross, iid, chr) %>%
	summarise(nxo = n()-1) %>%
	summarise(nxo = sum(nxo)) %>%
	left_join(samples) %>% 
	mutate(X_chrom = case_when(grepl("^G", cross) ~ "mus",
							   TRUE ~ "not mus"))

ggplot(xos_summ) +
	ggbeeswarm::geom_quasirandom(aes(x = cross, y = nxo), width = 0.25) +
	facet_grid(. ~ platform)

ggplot(xos_summ) +
	ggbeeswarm::geom_quasirandom(aes(x = cross, y = nxo, colour = age),
								 width = 0.4, dodge.width = 0.2) +
	facet_grid(. ~ platform)

## plot inter-crossover distance for double-crossovers
dummy <- lapply(names(CROSSES), function(f) tibble(width = 100*rexp(1000, 1), cross = f, group = "no interference") ) %>% 
	bind_rows()

subset(xo_all, censor == 0) %>%
	mutate(group = "observed") %>%
	bind_rows(dummy) %>%
	mutate(cross = factor(cross, levels = names(CROSSES), labels = CROSSES)) %>%
	ggplot() +
	stat_ecdf(aes(x = width, colour = group)) +
	#stat_function(fun = pexp, args = list(rate = 1), colour = "grey70") +
	scale_x_continuous("distance between COs (cM)", breaks = seq(0,100,50)) +
	scale_y_continuous("cumulative proportion", breaks = seq(0,1,0.5)) +
	scale_colour_manual("", values = c("grey70","darkblue")) +
	coord_cartesian(xlim = c(0, 105)) +
	facet_wrap(~ cross, nrow = 3) +
	theme(legend.position = "top")
ggsave("figures/interco_ecdfs_by_cross.pdf", width = 4, height = 4)

## summary stats for distribution of inter-CO distances
subset(xo_all, censor == 0) %>%
	mutate(cross = factor_cross(cross, long_names = TRUE)) %>%
	group_by(cross) %>%
	summarize(mu = mean(width),
			  med = median(width),
			  lolo = quantile(width, 0.025),
			  lo = quantile(width, 0.25),
			  hi = quantile(width, 0.75),
			  hihi = quantile(width, 0.975)) %>%
	knitr::kable(format = "latex", digits = 1)

## plot raw inter-xo distances by chromosome
p1 <- subset(xo_all, censor == 0) %>%
	left_join( chromsizes_mm39() %>% enframe(name = "chr", value = "chrlen") ) %>%
	mutate(chr = factor_chrom(chr, short = TRUE)) %>%
	ggplot(aes(x = chr, y = width)) +
	tidybayes::stat_pointinterval() +
	scale_y_continuous("inter-CO distance (cM)") +
	scale_x_discrete("\nchromosome")

chrlen <- chromsizes_cM() %>%
	enframe(name = "chr", value = "chrlen") %>%
	mutate(chr = factor_chrom(gsub("^chr","", chr), short = TRUE))
p2 <- subset(xo_all, censor == 0) %>%
	mutate(chr = factor_chrom(chr, short = TRUE)) %>%
	left_join( chrlen ) %>%
	ggplot(aes(x = chrlen, y = width)) +
	geom_smooth(method = "lm", fullrange = FALSE, colour = "grey50") +
	tidybayes::stat_pointinterval() +
	scale_y_continuous("inter-CO distance (cM)") +
	scale_x_continuous("\n chromosome length (cM)")

cowplot::plot_grid(p1, p2, nrow = 2, align = "v")
ggsave("figures/interco_dist_by_chrom.pdf", width = 4, height = 6)

assign_marker_physical_pos <- function(x, themap) {
	ii <- match(x, themap$marker)
	return( themap$pos[ii] )
}

as_tibble.GRanges <- function(x, ...) {
	rez <- as.data.frame(x) %>% as_tibble()
	colnames(rez)[1] <- "chr"
	return(rez)
}

xo_mini <- lapply( crosses_mini_clipped, genotypes_to_crossovers, return_breaks_table = TRUE) %>% bind_rows()
xo_mini$pos_left <- assign_marker_physical_pos(xo_mini$left, mm_mini)
xo_mini$pos_right <- assign_marker_physical_pos(xo_mini$right, mm_mini)
xo_mini$pos <- with(xo_mini, (pos_left + pos_right)/2)
xo_mega <- lapply( crosses_mega_clipped, genotypes_to_crossovers, return_breaks_table = TRUE) %>% bind_rows()
xo_mega$pos_left <- assign_marker_physical_pos(xo_mega$left, mm_mega)
xo_mega$pos_right <- assign_marker_physical_pos(xo_mega$right, mm_mega)
xo_mega$pos <- with(xo_mega, (pos_left + pos_right)/2)

make_genome_windows <- function(window_size = 1e6, step_size = window_size/2, ...) {
	mouser::chromsizes_mm39(as.seqinfo = TRUE) %>%
		GenomicRanges::GRanges() %>%
		GenomicRanges::slidingWindows(width = window_size, step = step_size) %>%
		unlist()
}

calculate_recomb_rate_on_grid <- function(xos, window_size = 1e6, step_size = window_size/2) {
	
	windows <- make_genome_windows(window_size, step_size)
	xo_gr <- with(xos, GenomicRanges::GRanges(paste0(chr, ":", pos)))
	windows_df <- as_tibble(windows)
	windows_df$n_xo <- GenomicRanges::countOverlaps(windows, xo_gr)
	return(windows_df)
	
}

## scaling factors -- number of meioses by agexcross combos
samples$genotyped <- samples$iid %in% nxo_all$iid
n_meioses <- subset(samples, genotyped) %>%
	group_by(cross, age) %>%
	summarise(n_meioses = n())
n_meioses_cross <- subset(samples, genotyped) %>%
	group_by(cross) %>%
	summarise(n_meioses = n())
n_meioses_age <- subset(samples, genotyped) %>%
	group_by(age) %>%
	summarise(n_meioses = n())

n_meioses_tbl <- n_meioses_cross %>%
	mutate(mat = factor_strains(substr(cross, 1, 1)),
		   pat = factor_strains(substr(cross, 2, 2))) %>%
	pivot_wider(names_from = "pat", values_from = "n_meioses")


## gather up all crossovers
xo_all <- bind_rows(xo_mega, xo_mini)
xo_all <- left_join(xo_all, subset(samples, genotyped))

## generate sliding-window estimates of recomb rate by cross
xo_rate <- xo_all %>%
	group_by(age) %>%
	do(calculate_recomb_rate_on_grid(., 5e6, 1e6)) %>%
	left_join(n_meioses_age) %>%
	mutate(cM_Mb = (100*n_xo/n_meioses)/(width/1e6))
xo_rate$age <- factor(xo_rate$age, names(AGES), labels = AGES)

## plot it along genome, Manhattan-plot style
xo_rate <- linearize_genome(xo_rate, chromsizes_mm39())
subset(xo_rate, is_autosome(chr)) %>%
	ggplot() +
	geom_line(aes(x = .start, y = cM_Mb, colour = age:.colour, group = age)) +
	facet_grid(age ~ .) +
	scale_y_continuous("cM/Mb\n") +
	scale_x_linearized_genome("", chrlen = chromsizes_mm39()[1:19]) +
	scale_colour_manual(values = rev(MetBrewer::met.brewer("Demuth", 4)), guide = "none")
ggsave("figures/windowed_recomb_rate_by_age.pdf", width = 6, height = 4)

## example plot of a single genome, for figure for TAGC talk
tmp <- genotypes_to_crossovers(crosses_mini_clipped$FG)
subset(tmp, iid == "PAE_036_FG_E_028") %>% 
	mutate(chr = factor_chrom(chr, short = TRUE),
		   geno = factor(geno)) %>% 
	ggplot() + 
	geom_segment(aes(x = lo, xend = hi, y = chr, yend = chr, colour = geno),
				 size = 2) +
	scale_y_discrete("", limits = rev) +
	scale_x_continuous("position (cM)") +
	scale_colour_manual(values = c("grey70","darkblue"), guide = "none") +
	theme_classic() +
	theme(axis.line.y = element_blank(),
		  axis.ticks.y = element_blank())
ggsave("figures/individual_haps_example.pdf", width = 2.5, height = 2.75)

## === ##

g <- as_tibble(crosses_mini$FG, mm_mini)
subset(g, chr == "chr6")

plot_coeff_quick <- function(m, ...) {
	
	x <- lsmeans::lsmeans(m, attr(m$terms, "term.labels")[1],
						  type = "response")
	print( plot(x) )
	
}

m0_glm <- glm(n_xo ~ 0+age, nxo_all, family = poisson())
m1_glm <- glm(n_xo ~ cross, nxo_all, family = poisson())
m2_glm <- glm(n_xo ~ age+cross, nxo_all, family = poisson())
m3_glm <- glm(n_xo ~ cross*age, nxo_all, family = poisson())
m4_glm <- glm(n_xo ~ age + mat + pat, nxo_all, family = poisson())
m0 <- lm(n_xo ~ 0+age, nxo_all)
m1 <- lm(n_xo ~ cross, nxo_all)
m2 <- lm(n_xo ~ age+cross, nxo_all)
m3 <- lm(n_xo ~ cross*age, nxo_all)
m4 <- lm(n_xo ~ mat+pat+age, nxo_all)
# m0 <- lme4::lmer(n_xo ~ 0+age + (1|sire), nxo_all)
# m1 <- lme4::lmer(n_xo ~ cross + (1|sire), nxo_all)
# m2 <- lme4::lmer(n_xo ~ age+cross + (1|sire), nxo_all)
# m3 <- lme4::lmer(n_xo ~ cross*age + (1|sire), nxo_all)
# m4 <- lme4::lmer(n_xo ~ age + mat + pat + (1|sire), nxo_all)

ALPHA <- 0.05 #0.05/length(CROSSES) # multiple testing correction
age_effs <- list(glm = marginaleffects::avg_comparisons(m3_glm, type = "response",
														variables = list(age = "pairwise"),
														by = "cross", conf_level = 1-ALPHA),
				 lm = marginaleffects::avg_comparisons(m3, type = "response",
				 									  variables = list(age = "pairwise"),
				 									  by = "cross", conf_level = 1-ALPHA))
age_effs <- lapply(age_effs, as_tibble) %>% bind_rows(.id = "model")

ggplot(age_effs) +
	geom_pointrange(aes(x = cross, y = estimate, ymin = conf.low, ymax = conf.high, colour = model),
					position = position_dodge(width = 0.3)) +
	geom_hline(yintercept = 0, lty = "dashed", colour = "grey70") +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
	ggtitle("paternal age effect by cross")

cross_effs <- list(glm = marginaleffects::avg_comparisons(m3_glm, type = "response",
														variables = list(cross = "pairwise")),
				 lm = marginaleffects::avg_comparisons(m3, type = "response",
				 									  variables = list(cross = "pairwise")))
cross_effs <- lapply(cross_effs, as_tibble) %>% bind_rows(.id = "model")
cross_effs <- left_join(cross_effs, rename(CROSSES_PAIRWISE_TYPES, contrast = code))

ggplot(cross_effs) +
	geom_pointrange(aes(x = name, y = estimate, ymin = conf.low, ymax = conf.high, colour = model),
					position = position_dodge(width = 0.3)) +
	geom_hline(yintercept = 0, lty = "dashed", colour = "grey70") +
	facet_grid(. ~ group, scale = "free_x", space = "free_x") +
	theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
	ggtitle("pairwise comparisons between crosses")

marginaleffects::avg_comparisons(m4, type = "response",
								 variables = list(mat = "pairwise",
								 				  pat = "pairwise")) %>%
	ggplot() +
	geom_pointrange(aes(x = contrast, y = estimate, ymin = conf.low, ymax = conf.high),
					position = position_dodge(width = 0.3)) +
	facet_grid(term ~ .) +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(nxo_all) +
	geom_histogram(aes(x = n_xo, y = after_stat(density)),
				   binwidth = 1) +
	facet_grid(pat ~ mat)

#m3_bayes <- brms::brm(n_xo ~ age*cross, data = nxo_all, family = poisson())

## for crossover interference analyses, use adapted versions of stuff from Karl's xoi package
xoloc <- c( lapply(crosses_mini_clipped, genotypes_to_crossovers),
			lapply(crosses_mega_clipped, genotypes_to_crossovers) )

## to fit the Stalh model, need data that is a list with one entry per meiotic product (= chromosome)
## first make nested list (iid>chrom), can extract invididuals by name later
breaks_all <- c( lapply(crosses_mini_clipped, prep_breaks_by_individual),
				 lapply(crosses_mega_clipped, prep_breaks_by_individual)
				 ) %>% Reduce("c", .)
## TODO: figure out how to carry over cross-specific chromosome lengths?
## for now just use common lengths
L <- c( lapply(crosses_mini_clipped, function(f) sapply(f$gmap, max)),
		lapply(crosses_mega_clipped, function(f) sapply(f$gmap, max)) ) %>%
	simplify2array() %>%
	apply(., 1, max)

## once crossover data is sanitized, can re-fit interference models at will
xoloc_all <- bind_rows(xoloc) %>% left_join(samples)
write_tsv(xoloc_all, "xos_inferred/xoloc_all.txt", col_names = FALSE)
for (cc in unique(xoloc_all$cross)) {
	fn <- paste0("xos_inferred/xoloc_", cc,".txt")
	message(cc, " --> ", fn)
	subset(xoloc_all, cross == cc) %>%
		select(d, censor, iid) %>%
		write_tsv(., fn, col_names = FALSE)
}

## save breaks in format that can eventually be used by xoi::fitStahl()
saveRDS(breaks_all, "xos_inferred/breaks_all.rds")
saveRDS(L, "xos_inferred/chromlen.rds")

x0 <- fit_xoi_gamma(xoloc_all)
x <- xoloc_all %>%
	group_by(age) %>%
	do(fit_xoi_gamma(.))
xx <- xoloc_all %>%
	group_by(cross) %>%
	do(fit_xoi_gamma(.))
xxx <- xoloc_all %>%
	group_by(cross,age) %>%
	do(fit_xoi_gamma(.))

compare_gamma_models(x0, x)
compare_gamma_models(x0, xx)
compare_gamma_models(xx, xxx)
compare_gamma_models(x, xxx)

xxx %>% 
	ggplot() + 
	geom_pointrange(aes(x = cross, y = nu, ymin = nu-1*se, ymax = nu+1*se,
						colour = age),
					position = position_dodge(width = 0.3))


stahl_fits <- samples %>%
	subset(iid %in% xo_all$iid) %>%
	group_by(cross) %>%
	do(get_breaks_and_fit_stahl(.$iid, breaks_all, L, nu_init = c(1,42)))

