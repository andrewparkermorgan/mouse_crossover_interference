## convert raw MegaMUGA genotypes to R/qtl2 cross objects for backcross offspring from PAE experiment
setwd("~/Dropbox/pmdvlab/xoi")

library(tidyverse)

devtools::load_all("~/Dropbox/pmdvlab/mouser/")
devtools::load_all("~/Dropbox/pmdvlab/popcorn/")
source("geno_helpers.R")

## Part 1: MegaMUGA

## new MM map
mm_mega <- read_map("array_annot/mega_mm39.map")

## apply updates to marker map for one argyle::genotypes object
do_updates <- function(f, new_map, ...) {
	old_map <- attr(f, "map")
	new_map <- update_map_positions(old_map, new_map)
	new_map$chr <- shorten_chrom_names(new_map$chr)
	if (!is.character(f))
		rez <- recode_genotypes_to_character(f[ new_map$marker, ])
	else
		rez <- f[ new_map$marker, ]
	rez <- update_genotypes_map(rez, new_map)
	return(rez)
}

## raw MM data for backcross offspring plus FVB/NJ strain; we only need this for FVB/NJ
all_geno_raw <- readRDS("~/Dropbox/PAE/genotypes/mega.rds")
all_geno_raw <- do_updates(all_geno_raw, mm_mega)
x <- all_geno_raw[ ,"FVB/NJm35257" ]

## consensus MM genotypes for CC founder strains
founders <- readRDS("array_annot/mega_CC_founders_consensus_mm39.rds")
founders <- do_updates(founders, mm_mega)

## add FVB/NJ to founder strains
founders <- merge_genotypes_unsafely(founders, x)
founders <- rename_individuals(founders, c(LETTERS[1:8],"V"))
saveRDS(founders, "geno/mega_founders_preclean.rds")

## raw MM data for backcross offspring, split by cross
crosses <- readRDS("~/Dropbox/PAE/genotypes/allcrosses.raw.rds")
crosses <- lapply(crosses, do_updates, attr(founders, "map"))
crosses <- lapply(crosses, function(f) f[ rownames(founders), ] )
saveRDS(crosses, "geno/crosses_mega_preclean.rds")

## check final dimensions
nrow(founders)
sapply(crosses, dim)

## global genotype QC
prop_missing <- p_missing(all_geno_raw)
#prop_maf <- calc_maf(all_geno_raw)

# note: we only want autosomes for crossover analysis
mk_non_autos <- subset(mm_mega, !is_autosome(chr))$marker

## list of markers to drop
mk_bad <- union(
	names( which(prop_missing > 0.1) ),
	mk_non_autos
)

## save list of markers to drop
writeLines(mk_bad, "qc/drop_markers_mega.txt")

clean_cross_genotypes <- function(dam, sire) {
	
	f1_fwd <- paste0(dam,sire, collapse = "")
	f1_rev <- paste0(sire,dam, collapse = "")
	message("Performing genotype QC and conversion on crosses ", f1_fwd, " and ", f1_rev, ".")
	
	## recode genotypes for backcrosses, pooling across reciprocal F1 parents
	xx <- merge_genotypes_unsafely( subset_genotypes(crosses[[f1_fwd]], is_autosome(chr)),
									subset_genotypes(crosses[[f1_rev]], is_autosome(chr))) %>%
		recode_as_backcross(.,
							subset_genotypes(founders, is_autosome(chr)),
							dam,sire,"V")
	
	message("Before genotype QC, there are ", nrow(xx), " markers.")
	
	## QC on genotype proportions; expected proportion is 0.50
	geno_counts <- tally_genotype_calls(xx)
	geno_props <- geno_counts/rowSums(geno_counts)
	geno_max_prop <- apply(geno_props, 1, max)

	## set some thresholds to declare a marker suspicious, by manual inspection
	sus_mk <- geno_max_prop < 0.2 | geno_max_prop > 0.8 | geno_props[,"N"] > 0.1
	## more markers to drop
	sus_mk <- union( names(which(sus_mk)), mk_bad )
	writeLines(sus_mk, paste0("qc/drop_per_cross_mega_",dam,sire,".txt"))
	
	mk_keep <- setdiff( rownames(xx),sus_mk )
	message("After genotype QC, there will be ", length(mk_keep), " markers.")
	
	## now go back and recode each cross separately
	crosses_recoded <- vector("list", 2L)
	names(crosses_recoded) <- c(f1_fwd, f1_rev)
	
	message("Final recoding genotypes for ", f1_fwd, " (", ncol(crosses[[f1_fwd]]), " individuals) ...")
	crosses_recoded[[f1_fwd]] <- recode_as_backcross( crosses[[f1_fwd]][ mk_keep, ],
													  founders[ mk_keep, ],
													  dam, sire, "V" )
	
	message("Final recoding genotypes for ", f1_rev, " (", ncol(crosses[[f1_rev]]), " individuals) ...")
	crosses_recoded[[f1_rev]] <- recode_as_backcross( crosses[[f1_rev]][ mk_keep, ],
													  founders[ mk_keep, ],
													  sire, dam, "V" )

	## convert to R/qtl2 format 	
	crosses_final <- lapply(crosses_recoded, convert_backcross_to_qtl2)
	
	return(crosses_final)
	
}

x1 <- clean_cross_genotypes("F","G")
x2 <- clean_cross_genotypes("G","H")
saveRDS( c(x1,x2), "geno/crosses_mega_clean.rqtl.rds" )


