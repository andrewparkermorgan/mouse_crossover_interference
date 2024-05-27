## helper functions for getting crossover locations from R/qtl2

## some constants
CROSSES <- c("FG" = "CASTxPWK",
			 "GF" = "PWKxCAST",
			 "FH" = "CASTxWSB",
			 "HF" = "WSBxCAST",
			 "GH" = "PWKxWSB",
			 "HG"=  "WSBxPWK")
AGES <- c("Y" = "young", "E" = "old")

CROSSES_PAIRWISE <- c(
	"GF - FG" = "PWKxCAST - CASTxPWK",
	"HF - FH" = "WSBxCAST - CASTxWSB",
	"HG - GH" = "WSBxPWK - PWKxWSB",
	"FH - FG" = "CASTxWSB - CASTxPWK",
	"FH - GF" = "CASTxWSB - PWKxCAST",
	"GH - FG" = "PWKxWSB - CASTxPWK",
	"GH - FH" = "PWKxWSB - CASTxWSB",
	"GH - GF" = "PWKxWSB - PWKxCAST",
	"GH - HF" = "PWKxWSB - WSBxCAST",
	"HF - FG" = "WSBxCAST - CASTxPWK",
	"HF - GF" = "WSBxCAST - PWKxCAST",
	"HG - FG" = "WSBxPWK - CASTxPWK",
	"HG - FH" = "WSBxPWK - CASTxWSB",
	"HG - GF" = "WSBxPWK - PWKxCAST",
	"HG - HF" = "WSBxPWK - WSBxCAST"
)

CROSSES_PAIRWISE_TYPES <- tibble(code = names(CROSSES_PAIRWISE),
								 name = unname(CROSSES_PAIRWISE),
								 group = ifelse(seq_along(CROSSES_PAIRWISE) <= 3, "reciprocals","crosswise"))

factor_cross <- function(x, long_names = FALSE, ...) {
	if (long_names)
		factor(x, levels = names(CROSSES), labels = CROSSES)
	else
		factor(x, levels = names(CROSSES))
}

factor_age <- function(x, long_names = FALSE, ...) {
	if (long_names)
		factor(x, levels = names(AGES), labels = AGES)
	else
		factor(x, levels = names(AGES))
}

factor_strains <- function(x, long_names = FALSE, ...) {
	if (long_names)
		factor(x, levels = names(CC_STRAINS()[c(8,7,6)]), CC_STRAINS()[c(8,7,6)])
	else
		factor(x, levels = names(CC_STRAINS()[c(8,7,6)]))
}

scale_x_crosses <- function(...) {
	scale_x_discrete(labels = CROSSES)
}

scale_x_strains<- function(...) {
	scale_x_discrete(labels = CC_STRAINS()[c("F","G","H")])
}

## this function modified from R/qtl2
geno_to_ivls <- function (geno, map) 
{
	if (all(is.na(geno))) 
		return(NULL)
	stopifnot(length(geno) == length(map))
	map <- map[!is.na(geno)]
	geno <- geno[!is.na(geno)]
	d <- diff(geno)
	rl <- rle(geno)
	xo_int <- which(d != 0)
	data.frame(lo = map[ c(1, xo_int + 1) ],
			   hi = map[ c(xo_int, length(map)) ],
			   geno = geno[ c(xo_int, length(map)) ],
			   n_markers = rl$lengths,
			   left_marker = names(map)[ c(1, xo_int + 1) ],
			   right_marker = names(map)[ c(xo_int, length(map)) ]
			   )
}

## define genotype blocks from Vitberbi genotypes
define_xos <- function(geno, gmap, return_breaks = FALSE, return_breaks_table = FALSE, ...) {
	
	## make table of genotype intervals
	.find_xo_per_chrom_ind <- function(X, this_map) {
		ivl <- geno_to_ivls(X, this_map)
		ivl_df <- tibble::as_tibble(ivl)
		ivl_df$censor <- 3
		if (nrow(ivl_df) > 1) {
			ivl_df$censor[1] <- 1
			ivl_df$censor[ nrow(ivl_df) ] <- 2
			if (nrow(ivl_df) > 2) {
				ivl_df$censor[ 2:(nrow(ivl_df)-1) ] <- 0
			}
		}
		return(ivl_df)
	}
	
	## make list of breakpoint locations between intervals
	## (this is for xoi::fitStahl and similar)
	.make_break_list_per_chrom_ind <- function(X, this_map) {
		ivl <- geno_to_ivls(X, this_map)
		if (nrow(ivl) == 1) {
			return(numeric(0))
		}
		else {
			his <- ivl$hi[ -nrow(ivl) ]
			los <- ivl$lo[ -1 ]
			mids <- (los+his)/2
			#print(ivl)
			#print(mids)
			return(mids)
		}
	}
	
	## make table of crossovers (just breakpoints, not intervals with censoring info)
	.make_break_table <- function(X, this_map) {
		ivl <- geno_to_ivls(X, this_map)
		ivl_df <- tibble::as_tibble(ivl)
		if (nrow(ivl_df) > 1) {
			mk_left <- ivl_df$right_marker[ -nrow(ivl_df) ]
			mk_right <- ivl_df$left_marker[ -1 ]
			cM_mid <- with(ivl_df, (hi[ -nrow(ivl_df) ] + lo[ -1 ])/2)
			return( tibble(left = mk_left, right = mk_right, cM_mid = cM_mid) )
		}
		else {
			return(NULL)
		}
	}
		
	## input is a list of genotype matrices, one per chrom
	rez <- vector("list", length = length(geno))
	for (ii in seq_along(geno)) {
		if (return_breaks) {
			this_chrom <- apply(geno[[ii]], 1, .make_break_list_per_chrom_ind, gmap[[ii]])
			rez[[ii]] <- this_chrom
		}
		else if (return_breaks_table) {
			this_chrom <- apply(geno[[ii]], 1, .make_break_table, gmap[[ii]])
			rez[[ii]] <- dplyr::bind_rows(this_chrom, .id = "iid")			
		}
		else {
			this_chrom <- apply(geno[[ii]], 1, .find_xo_per_chrom_ind, gmap[[ii]])
			rez[[ii]] <- dplyr::bind_rows(this_chrom, .id = "iid")			
		}

	}

	if (!return_breaks & !return_breaks_table) {
		rez <- dplyr::bind_rows(rez, .id = "chr")
		rez$width <- with(rez, hi-lo)
	}
	else if (return_breaks_table) {
		rez <- dplyr::bind_rows(rez, .id = "chr")
	}
	
	return(rez)
	
}

genotypes_to_crossovers <- function(x, hmm_error_prob = 0.01, return_breaks = FALSE, return_breaks_table = FALSE, ...) {
	
	## run HMM to get genotype probs, then get Viterbi solution and use that for analysis
	message("Estimating genotype probabilities ... ")
	geno_smooth <- qtl2::viterbi(x, error_prob = hmm_error_prob, lowmem = TRUE, quiet = TRUE)
	
	## convert genotypes to haplotype intervals; implicitly defines crossovers
	message("Enumerating crossovers ... ")
	xos <- define_xos(geno_smooth, x$gmap, return_breaks = return_breaks, return_breaks_table = return_breaks_table)
	
	return(xos)
	
}

just_count_crossovers <- function(x, hmm_error_prob = 0.01, ...) {
	
	## run HMM to get genotype probs, then get Viterbi solution and use that for analysis
	message("Estimating genotype probabilities ... ")
	geno_smooth <- qtl2::viterbi(x, error_prob = hmm_error_prob, lowmem = TRUE, quiet = TRUE)
	
	## get crossover count (ind x chr matrix)
	message("Counting crossovers ...")
	xos <- qtl2::count_xo(geno_smooth)
	
	return(xos)
	
}

## intersect marker maps
intersect_qtl2_maps <- function(m1, m2) {
	
	chroms <- intersect(names(m1), names(m2))
	rez <- vector("list", length(chroms))
	names(rez) <- chroms
	for (cc in chroms) {
		nm <- intersect( names(m1[[cc]]), names(m2[[cc]]) )
		loc1 <- m1[[cc]][nm]
		loc2 <- m2[[cc]][nm]
		keeps <- unname(abs(loc1-loc2) < 1e-5)
		nm <- nm[ keeps ]
		rez[[cc]] <- sort(loc1[nm])
	}
	
	return(rez)
	
}

## compare marker maps; keep those in m1 but not m2
complement_qtl2_maps <- function(m1, m2) {
	
	chroms <- intersect(names(m1), names(m2))
	rez <- vector("list", length(chroms))
	names(rez) <- chroms
	for (cc in chroms) {
		nm <- setdiff( names(m1[[cc]]), names(m2[[cc]]) )
		loc1 <- m1[[cc]][nm]
		rez[[cc]] <- sort(loc1[nm])
	}
	
	return(rez)
	
}

subset_cross <- function(X, mm, gmap = NULL, ...) {
	
	mm <- sapply(mm, names) %>% Reduce(c, .)
	
	if ("viterbi" %in% class(X)) {
		
		for (cc in names(X)) {
			this_mm <- colnames(X[[cc]])
			keeps <- unname(this_mm %in% mm)
			X[[cc]] <- X[[cc]][ ,keeps,drop=FALSE ]
			gmap[[cc]] <- gmap[[cc]][ colnames(X[[cc]]) ]
		}
	
		attr(X, "gmap") <- gmap
		return(X)
			
	}
	
}

map_to_df <- function(X, annot = NULL, ...) {

	.map_to_df <- function(m) {
		m <- tibble::enframe(m) 
		colnames(m) <- c("marker","cM_for_qtl")
		m$in_map <- TRUE
		return(m)
	}

	rez <- lapply(X$gmap, .map_to_df) %>% bind_rows(.id = "chr")
	if (!is.null(annot))
		rez <- left_join(annot, rez[,c("marker","cM_for_qtl","in_map")])
	return(rez)
	
}

## prune genotypes in DENSE cross to nearest marker in SPARSE map
prune_to_nearest_marker <- function(X, sparse_map, ...) {
	
	message("Starting with ", qtl2::tot_mar(X) ," markers ...")
	
	mm <- vector("character", 0L)
	for (cc in names(X$geno)) {
		chrom <- rep(cc, length(sparse_map[[cc]]))
		mm <- c(mm, qtl2::find_marker(X$gmap, chrom, sparse_map[[cc]] ))
	}
	
	mm <- unique(mm)
	message("Trying to extract ", length(mm), " markers ...")

	rez <- qtl2::pull_markers(X, mm)
	message("Ending with ", qtl2::tot_mar(rez), " markers.")
		
	return(rez)
	
}

## calculate distance in DENSE map to nearest marker in SPARSE map
dist_to_nearest_marker <- function(X, sparse_map, ...) {
	
	rez <- vector("list", length(X$gmap))
	names(rez) <- names(X$gmap)
	
	for (cc in names(X$gmap)) {
		chrom <- rep(cc, length(sparse_map[[cc]]))
		mm <- qtl2::find_marker(X$gmap, chrom, sparse_map[[cc]])
		mm_pos <- qtl2::find_markerpos(X, unique(mm))
		rez[[cc]] <- tibble(chr = chrom, query = sparse_map[[cc]],
							marker = mm, pos = mm_pos[mm,2])
		rez[[cc]]$gap <- with(rez[[cc]], abs(query-pos))
	}
	
	return( bind_rows(rez) )
	
}


## replacement for xoi::find.breaks that works with R/qtl2
## this provides input for xoi::fitStahl
## NB: this assumes backcross data for now
## result is nested list: [ chroms ][ iids ]  
find_breaks <- function(geno, map = NULL, crosstype = "bc", is_x_chr = FALSE, ...) {
	
	if (is.null(map))
		map <- geno$gmap
	
	gprobs <- qtl2::viterbi(geno, ..., lowmem = TRUE)
	
	nchr <- length(gprobs)
	rez <- vector("list", nchr)
	names(rez) <- names(gprobs)
	
	for (ii in seq_along(gprobs)) {
		
		breaks <- qtl2:::.locate_xo(t(gprobs[[ii]]), map[[ii]], crosstype, is_x_chr)
		names(breaks) <- rownames(gprobs[[ii]])
		rez[[ii]] <- breaks
		attr(rez[[ii]], "L") <- max(map[[ii]])
		
	}
	
	attr(rez, "L") <- sapply(map, max)
	return(rez)
	
}

## adapted from xoi::convertxoloc, but keep track of individual IDs
convert_xo_loc <- function(breaks) {
	
	f <- function(x, L) {
		if (length(x) == 0) 
			return(tibble(d = L, censor = 3L))
		else {
			d <- diff(c(0, x, L))
			cen <- c(2L, rep(0L, length(x) - 1), 1L)
			return(tibble(d = d, censor = cen))
		}
	}
	if (is.list(breaks[[1]])) {
		v <- vector("list", length(breaks))
		names(v) <- names(breaks)
		for (i in 1:length(breaks)) {
			v[[i]] <- lapply(breaks[[i]], f, attr(breaks[[i]], "L"))
			rez <- bind_rows(v[[i]], .id = "iid")
			v[[i]] <- rez
		}
		v <- bind_rows(v, .id = "chr")
		v <- v[ ,c("d","censor","chr","iid") ]
	}
	else {
		v <- lapply(breaks, f, attr(breaks, "L"))
		v <- matrix(unlist(v), ncol = 2, byrow = TRUE)
	}
	#v <- as.data.frame(v)
	#names(v) <- c("distance", "censor")
	return(v)
	
}

convert_simulated_xo_loc <- function(breaks) {
	
	f <- function(x, L) {
		if (length(x) == 0) 
			return(tibble(d = L, censor = 3L))
		else {
			d <- diff(c(0, x, L))
			cen <- c(2L, rep(0L, length(x) - 1), 1L)
			return(tibble(d = d, censor = cen))
		}
	}
	
	v <- lapply(breaks, f, attr(breaks, "L"))
	v <- bind_rows(v, .id = "iid")
	v <- v[ ,c("d","censor","iid") ]
	
	#v <- as.data.frame(v)
	#names(v) <- c("distance", "censor")
	return(v)
	
}

## take results from convert_xo_loc (or similar) and fit gamma model
fit_xoi_gamma <- function(df, lo = 1, hi = 20, return_se = TRUE, return_interval = FALSE, ...) {
	
	message("Obtaining rough estimate of nu (range ", lo, ", ... , ", hi, " ) ...")
	rough <- xoi::fitGamma(d = df$d, censor = df$censor,
						   nu = seq(lo, hi, by = 2))
	nu0 <- rough[ which.max(rough[,2]),1 ]
	message("    nu0 = ", nu0)
	message("    loglik = ", max(rough[,2]))
	message("Refining estimate of nu ...")
	fine <- xoi::fitGamma(d = df$d, censor = df$censor,
						   lo = max(nu0 - 2,0), hi = nu0 + 2,
						   se = return_se, supint = return_interval,
						   ...)
	return( as_tibble(fine) )
	
}

## fit Stahl model (gamma + escape) 
fit_xoi_stahl <- function(breaks, L = NULL, nu_init = c(1, 24), p_init = 0.1, ...) {
	
	rez <- xoi::fitStahl(breaks, chrlen = L, nu = nu_init, p = p_init, ...)
	names(rez) <- c("nu","p","loglik","nu_0","loglik_0","lnLR")
	rez <- as.list(rez)
	return( as_tibble(rez) )
	
}

compare_gamma_models <- function(m0, m1, ...) {
	
	if (nrow(m0) >= nrow(m1))
		stop("Models may be swapped: should be m0 = null/smaller, m1 = alternative/larger")
	
	lambda <- -2*(sum(m0$loglik) - sum(m1$loglik))
	names(lambda) <- "lambda"
	ndf <- nrow(m1)-nrow(m0)
	
	rez <- list(statistic = lambda,
				parameter = NULL,
				p.value = pchisq(lambda, ndf, lower.tail = FALSE),
				alternative = "one-sided",
				null.value = c("chisq", 1),
				method = "Likelihood ratio rest (chisq)",
				data.name = "m0 and m1")
	class(rez) <- c("htest", class(rez))
	return(rez)
	
}

compare_stahl_models <- function(m0, m1, ...) {
	
	if (nrow(m0) >= nrow(m1))
		stop("Models may be swapped: should be m0 = null/smaller, m1 = alternative/larger")
	
	## LRT for full models
	lambda <- -2*(sum(m0$loglik) - sum(m1$loglik))
	ndf <- 2*(nrow(m1)-nrow(m0))
	
	## LRT for the gamma sub-models
	lambda_0 <- -2*(sum(m0$loglik_0) - sum(m1$loglik_0))
	ndf_0 <- 1*(nrow(m1)-nrow(m0))
	
	rez <- tibble( statistic = c(lambda, lambda_0),
				   df = c(ndf, ndf_0),
				   p.value = pchisq(c(lambda, lambda_0), c(ndf, ndf_0), lower.tail = FALSE),
				   method = "LRT",
				   models = c("gamma-escape","gamma")
				   )
	
	return(rez)
	
}


sim_gamma_model <- function(nu, nchrom = 1000, L = 100, verbose = FALSE, obligate_chiasma = FALSE, ...) {
	
	breaks <- xoi::simStahl(nchrom, nu, L, p = 0, obligate_chiasma = obligate_chiasma)
	if (verbose)
		message("Simulated ", length(breaks), " meioses with nu = ", nu, " and L = ", L)
	breaks <- convert_simulated_xo_loc(breaks)
	return(breaks)
	
}

prep_breaks_by_individual <- function(geno, ...) {

	#breaks <- find_breaks(geno, ...)
	breaks <- genotypes_to_crossovers(geno, ..., return_breaks = TRUE)
	
	L <- attr(breaks, "L")
	n_chrom <- length(breaks)
	n_ind <- length(breaks[[1]])
	
	by_ind <- vector("list", length = n_ind)
	names(by_ind) <- names(breaks[[1]])
	for (ii in 1:n_ind) {
		by_ind[[ii]] <- vector("list", length = n_chrom)
	}
	for (ii in 1:n_chrom) {
		for (jj in 1:n_ind) {
			by_ind[[jj]][[ii]] <- breaks[[ii]][[jj]]
		}
	}
	
	#print(by_ind)
	attr(by_ind, "L") <- L
	return(by_ind)
	
}

lookup_breaks_by_individual <- function(breaks, iids, L = NULL, chroms = NULL, ...) {
	
	if (is.null(L))
		L <- attr(breaks, "L")
	
	if ( !all(iids %in% names(breaks)) )
		stop("Some of the supplied iids not found in breaks")
	
	rez <- Reduce("c", breaks[iids])
	lens <- rep(L, length(iids))
	return( list(pts = rez, L = lens) )
	
}

## try fitting the Stahl model
get_breaks_and_fit_stahl <- function(iids, breaks, L = NULL, ...) {
	
	message("Retreiving crossover events for ", length(iids), " individuals ...")
	b <- lookup_breaks_by_individual(breaks, iids, L)
	message("Fitting Stahl model from ", length(b$pts), " meiotic products ...")
	rez <- fit_xoi_stahl(b$pts, b$L, ...)
	rez$n_ind <- length(iids)
	rez$n_chrom <- length(b$pts)
	return(rez)
	
}

## fit Stahl model on bootstrap replicates
## first replicate (i == 0) is the real data
bootstrap_stahl <- function(df, n = 10, cores = 8, ...) {
	
	.do_boot <- function(k) {
		if (k == 0) {
			iids <- df$iid
		}
		else {
			ii <- sample.int(nrow(df), replace = TRUE)
			#print(ii)
			iids <- df$iid[ii]
		}
		message(" ... bootstrap replicate ", k, " of ", n)
		get_breaks_and_fit_stahl(iids, ...)
	}
	
	lapply(0:n, function(i) .do_boot(i)) %>%
		bind_rows(.id = "rep") %>%
		mutate(rep = as.integer(rep)-1)
	
}

missing_by_ind <- function(X, ...) {
	
	.n_miss <- function(g) {
		rowSums(is.na(g))
	}
	
	nn <- sapply(X$geno, .n_miss) %>% rowSums()
	nm <- qtl2::tot_mar(X)
	return( nn/nm )
	
}

missing_by_marker <- function(X, ...) {
	
	.n_miss <- function(g) {
		colSums(is.na(g))
	}
	
	nn <- lapply(X$geno, .n_miss) %>% Reduce("c", .)
	nm <- qtl2::n_ind(X)
	return( nn/nm )
	
}

calc_hpdi <- function(x, burnin = NULL, prob = 0.95, ...) {
	
	if (!is.null(burnin)) {
		xx <- coda::mcmc(x[ seq(burnin+1, nrow(x)),,drop=FALSE ])
		message("Original chain length ", nrow(x), "; after dropping ", burnin, " iterations, final chain length is ", nrow(xx))
	}
	else {
		xx <- coda::mcmc(x)
	}
		
	rez <- summary(xx)
	hpd <- coda::HPDinterval(xx, prob = prob, ...)
	
	rez <- tibble( term = rownames(rez$statistics),
			estimate = rez$statistics[,"Mean"],
			conf_lo = hpd[,1],
			conf_hi = hpd[,2])
	return(rez)
	
}

calc_hpdi_from_vector <- function(x, prob = 0.95, ...) {
	
	x <- as.vector(x)
	x <- matrix(x, ncol = 1)
	x <- coda::as.mcmc(x)
	hpd <- coda::HPDinterval(x, prob = prob, ...)
	rez <- tibble( estimate = mean(x, na.rm = TRUE),
				   conf_lo = hpd[,1],
				   conf_hi = hpd[,2])
	return(rez)
	
}

## from <https://github.com/walmes/wzRfun/blob/master/R/pairwise.R>
all_pairwise_comparisons_matrix <- function(lev, ...) {
	
	combos <- utils::combn(seq_along(lev), 2)
	M <- matrix(0, nrow = ncol(combos), ncol = length(lev))
	colnames(M) <- lev
	rownames(M) <- paste(lev[combos[1, ]], lev[combos[2, ]], sep = "-")
	
	for (ii in 1:ncol(combos)) {
		
		M[ ii,combos[1,ii] ] <- +1
		M[ ii,combos[2,ii] ] <- -1
		
	}
	
	return(M)
}
