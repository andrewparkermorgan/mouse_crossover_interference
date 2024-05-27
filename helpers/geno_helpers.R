## helper functions for genotypes objects
## this will eventually replace the `argyle` packaage, which I am in proces of rewriting

.is_mutable_genotypes <- function(G, ...) {
	
	rez <- TRUE
	
	if (!is.null(attr(G,"mutable")))
		rez <- as.logical(attr(G,"mutable"))
	
	return(rez)
	
}

.is_same_genotype_encoding <- function(G1, G2) {
	
	if (!.is_mutable_genotypes(G1) | !.is_mutable_genotypes(G2))
		return(FALSE)
	
	if (is.character(G1) & is.character(G2))
		return(TRUE)
	else if (is.numeric(G1) & is.numeric(G2))
		return(TRUE)
	else
		return(FALSE)
	
}

tally_genotype_calls <- function(G, ...) {
	
	G <- unclass(G)
	
	pb <- txtProgressBar(min = 0, max = nrow(G), style = 3)
	
	if (is.character(G)) {
	
		counts <- matrix(NA_integer_, nrow = nrow(G), ncol = 6)
			
		for (ii in seq_len(nrow(G))) {
			As <- sum(G[ii,] == "A", na.rm = TRUE)
			Cs <- sum(G[ii,] == "C", na.rm = TRUE)
			Gs <- sum(G[ii,] == "G", na.rm = TRUE)
			Ts <- sum(G[ii,] == "T", na.rm = TRUE)
			Hs <- sum(G[ii,] == "H", na.rm = TRUE)
			Ns <- sum(G[ii,] == "N" | is.na(G[ii,]),  na.rm = TRUE)
			this_row <- c(As,Cs,Gs,Ts,Hs,Ns)
			counts[ii,] <- this_row
			setTxtProgressBar(pb, ii)
		}
		
		rownames(counts) <- rownames(G)
		colnames(counts) <- c("A","C","G","T","H","N")
		
	}
	else if (is.numeric(G)) {
	
		counts <- matrix(NA_integer_, nrow = nrow(G), ncol = 4)
			
		for (ii in seq_len(nrow(G))) {
			hom_ref <- sum(G[ii,] == 0, na.rm = TRUE)
			het <- sum(G[ii,] == 1, na.rm = TRUE)
			hom_alt <- sum(G[ii,] == 2, na.rm = TRUE)
			Ns <- sum(is.na(G[ii,]),  na.rm = TRUE)
			this_row <- c(hom_ref,het,hom_alt,Ns)
			counts[ii,] <- this_row
			setTxtProgressBar(pb, ii)
		}
		
		rownames(counts) <- rownames(G)
		colnames(counts) <- c("0","1","2","N")
		
	}
	
	message("\n")
	return(counts)
	
}

guess_alleles <- function(G, ...) {
	
	NT <- c("A","C","G","T")
	
	nt_counts <- tally_genotype_calls(G)
	nt_counts <- nt_counts[ ,1:4,drop=FALSE ]
	alleles <- matrix(NA_character_, nrow = nrow(G), ncol = 2)
	for (ii in seq_len(nrow(G))) {
		nz <- which(nt_counts[ii,] > 0)
		na <- length(nz)
		if (na > 2) {
			alleles[ii,] <- c("X","X")
		}
		else if (na == 1) {
			alleles[ii,1] <- c(NT[nz])
		}
		else if (na == 2) {
			alleles[ii,] <- c(NT[nz])
		}
	}
	
	rownames(alleles) <- rownames(G)
	return(alleles)
	
}

`[.genotypes` <- function(G, i, j, ..., drop = FALSE) {
	
	if (missing(i))
		i <- rep(TRUE, nrow(G))
	if (missing(j))
		j <- rep(TRUE, ncol(G))
	
	r <- NextMethod("[", drop = drop)
	## handle "array indexing"
	#if (is.matrix(i)) {
	#	ii <- i
	#	i <- ii[,1, drop = TRUE]
	#	j <- ii[,2, drop = TRUE]
	#}
	r <- .reconcile_metadata(r, attr(G,"map"), attr(G, "ped"))
	return(r)
	
}

.reconcile_metadata <- function(G, mm = NULL, ped = NULL, ...) {
	
	ii <- rownames(G)
	jj <- colnames(G)
	
	if (is.null(mm))
		mm <- attr(G, "map")
	if (is.null(ped))
		ped <- attr(G, "ped")
	
	if (!all(ii %in% mm$marker))
		stop("Not all markers in genotypes are present in marker map.")
	mm <- mm[ match(ii, mm$marker),,drop=FALSE ]
	if (!tibble::is_tibble(mm))
		mm <- tibble::as_tibble(mm)
	
	if (!all(jj %in% ped$iid))
		stop("Not all markers in genotypes are present in marker map.")
	ped <- ped[ match(jj, ped$iid),,drop=FALSE ]
	if (!tibble::is_tibble(ped))
		ped <- tibble::as_tibble(ped)
	
	attr(G, "map") <- mm
	attr(G, "ped") <- ped
	
	if (!("genotypes" %in% class(G)))
		class(G) <- c("genotypes", class(G))
	
	return(G)
	
}

## find markers that are homozygous different genotypes between two focal individuals
find_informative_markers <- function(G, f1, f2) {
	
	G <- unclass(G)
	
	f1 <- G[ ,f1 ]
	f2 <- G[ ,f2 ]
	
	infm <- f1 != f2 & f1 != "H" & f1 != "H" & !is.na(f1) & !is.na(f2)
	return(infm)
	
}

recode_as_backcross <- function(G, G_ref, p1, p2, tester, ...) {
	
	G <- unclass(G)
	
	if (!all(rownames(G) == rownames(G_ref)))
		stop("Marker order betwen focal and reference samples does not match.")
	if (!.is_mutable_genotypes(G))
		stop("Genotypes cannot be recoded.")
	
	## find markers p1 != p2
	infm0 <- find_informative_markers(G_ref, p1, p2)
	G_infm <- G[ infm0,,drop=FALSE ]
	G_ref_infm <- G_ref[ infm0,,drop=FALSE ]
	
	## case 1: P1 != P2, P1 != tester
	infm1 <- find_informative_markers(G_ref_infm, p1, tester)
	
	## case 2: P1 != P2, P2 != tester
	infm2 <- find_informative_markers(G_ref_infm, p2, tester)
	
	## now create new genotypes
	G_new <- matrix(NA_integer_, nrow = nrow(G_ref_infm), ncol = ncol(G))
	rownames(G_new) <- rownames(G_ref_infm)
	colnames(G_new) <- colnames(G)
	for (ii in seq_len(ncol(G))) {
		
		is_het <- G_infm[,ii] == "H"
		
		## case 1
		G_new[ is_het & infm1,ii ] <- 1L
		G_new[ !is_het & infm1,ii ] <- 2L
		
		## case 2
		G_new[ is_het & infm2,ii ] <- 2L
		G_new[ !is_het & infm2,ii ] <- 1L
		
	}
	
	old_map <- attr(G, "map")
	attr(G_new, "map") <- old_map[ match(rownames(G_new), old_map$marker), ]
	attr(G_new, "ped") <- attr(G, "ped")
	class(G_new) <- c("genotypes", class(G_new))
	attr(G_new, "mutable") <- FALSE
	
	return(G_new)
	
}

convert_backcross_to_qtl2 <- function(G, ...) {
	
	G <- unclass(G)
	
	.process_chrom <- function(df) {
		
		m <- df$marker
		geno <- t( G[ m,,drop=FALSE ] )
		gmap <- setNames( df$cM, df$marker )
		pmap <- setNames( df$pos, df$marker )
		
		return( list(geno = geno, gmap = gmap, pmap = pmap) )
		
	}
	
	the_map <- attr(G, "map")
	by_chrom <- split(the_map, the_map$chr, drop = TRUE)
	nchrom <- length(by_chrom)
	
	rez <- list( crosstype = "bc", 
				 geno = vector("list", nchrom),
				 gmap = vector("list", nchrom),
				 pmap = vector("list", nchrom),
				 is_x_chr = rep(FALSE, nchrom), # fake
				 is_female = rep(FALSE, ncol(G)), # fake
				 cross_info = matrix(0L, ncol = 1, nrow = ncol(G), dimnames = list(colnames(G))),
				 alleles = c("A","B") )
	
	for (ii in seq_along(by_chrom)) {
		
		converted <- .process_chrom( by_chrom[[ii]] )
		rez$geno[[ii]] <- converted$geno
		rez$gmap[[ii]] <- converted$gmap
		rez$pmap[[ii]] <- converted$pmap
		
	}
	
	names(rez$geno) <- names(by_chrom)
	names(rez$gmap) <- names(by_chrom)
	names(rez$pmap) <- names(by_chrom)
	
	class(rez) <- c("cross2", class(rez))
	
	return(rez)
	
}

## subset genotypes by markers; could be by logical, index or name
slice_by_markers <- function(G, m, ...) {
	
	if (is.logical(m) & length(m) == nrow(G))
		m <- rownames(G)[ which(m) ]
	else if (is.numeric(m) & !is.logical(m)) {
		m_conf <- m
	}
	else if (is.character(m)) {
		not_in_geno <- !(m %in% rownames(G))
		m_conf <- m[ !not_in_geno ]
	}
	
	g <- unclass(G)[ m_conf,,drop=FALSE ]
	m_new <- rownames(g)
	
	map <- attr(G, "map")[ match(m_new, attr(G,"map")$marker),,drop=FALSE ]
	
	attr(g, "map") <- as_tibble(map)
	attr(g, "ped") <- as_tibble( attr(G, "ped") )
	class(g) <- c("genotypes", class(g))
	return(g)
	
}

## subset genotypes by individuals by name only
slice_by_ind <- function(G, ii, ...) {
	
	if (!is.numeric(ii) & !is.logical(ii))
		ii <- which(colnames(G) %in% ii)
	
	g <- unclass(G)[ ,ii,drop=FALSE ]
	
	attr(g, "map") <- as_tibble( attr(G, "map") )
	attr(g, "ped") <- as_tibble( attr(G, "ped")[ ii,,drop=FALSE ] )
	class(g) <- c("genotypes", class(g))
	return(g)
	
}

## convert genotypes from numeric (0=hom_ref,1=het,2=hom_alt) to character representation
recode_genotypes_to_character <- function(G, ...) {
	
	if (!is.numeric(G))
		stop("Genotypes are not encoded as numeric.")
	
	if (!.is_mutable_genotypes(G))
		stop("Genotypes cannot be recoded.")
	
	message("Recoding genotypes as character ...")
	
	G <- unclass(G)
	
	.recode <- function(g) {
		hom_ref <- which(g == 0)
		het <- which(g == 1)
		hom_alt <- which(g == 2)
		g[ hom_ref ] <- attr(G, "map")$A1[ hom_ref ]
		g[ het ] <- "H"
		g[ hom_alt ] <- attr(G, "map")$A2[ hom_alt ]
		return(g)
	}
	
	G_new <- matrix(NA_character_, nrow = nrow(G), ncol = ncol(G))
	rownames(G_new) <- rownames(G)
	colnames(G_new) <- colnames(G)
	
	for (ii in seq_len(ncol(G))) {
		G_new[ ,ii ] <- .recode(G[ ,ii ])
	}
	
	attr(G_new, "map") <- attr(G, "map")
	attr(G_new, "ped") <- attr(G, "ped")
	class(G_new) <- c("genotypes", class(G_new))
	return(G_new)
	
}

## update an old marker map with new genetic and physical positions; leave alleles the same
update_map_positions <- function(old_map, new_map, ...) {
	
	message("Old map has ", nrow(old_map), " markers.")
	
	# position of markers in both maps, on old map
	m_both <- which(old_map$marker %in% new_map$marker)
	message("New map will have ", length(m_both), " markers.")
	
	# position of markers in both maps, on new map
	ii <- match(old_map$marker[m_both], new_map$marker)
	
	# new map, in order of NEW map
	rez <- new_map[ ii, ]
	
	# still use OLD alleles
	rez$A1 <- old_map$A1[m_both]
	rez$A2 <- old_map$A2[m_both]
	
	# sort in order of NEW map
	rez$chr <- factor_chrom(rez$chr)
	rez <- arrange(rez, chr, cM, pos)
	
	return(rez)
	
}

## update marker map for a genotypes object
update_genotypes_map <- function(G, new_map, ...) {
	
	m_new <- new_map$marker
	if(!all(m_new %in% rownames(G)))
		stop("Not all markers in new map are in genotype matrix.")
	
	G_new <- G[ m_new,,drop=FALSE ]
	attr(G_new, "map") <- new_map
	attr(G_new, "ped") <- attr(G, "ped")
	class(G_new) <- c("genotypes", class(G_new))
	
	return(G_new)
	
}

add_individuals <- function(G, G_add, ...) {
	
	message("This will blindly concatenate genotype matrices column-wise and will NOT check that alleles are consistent.")
	
	if (!is.character(G) | !is.character(G_add))
		stop("Can only concatenate genotypes encoded as character (will allow to catch allele mix-ups).")
	
	m <- rownames(G)
	G_new <- matrix(NA_character_, nrow = nrow(G), ncol = ncol(G_add))
	rownames(G_new) <- rownames(G)
	colnames(G_new) <- colnames(G_add)
	G_add <- unclass(G_add) # get around S3 method dispatch for `[` operator on genotypes objects
	
	#print(head(G_add))
	#print(head(G))
	#print(head(G_new))
	
	## find intersecting markers
	m_idx <- match(m, rownames(G_add))
	in_both <- !is.na(m_idx)
	m_idx_both <- m_idx[in_both]
	
	for (jj in colnames(G_add)) {
		G_new[ which(in_both), jj ] <- G_add[ m_idx_both,jj ]
	}
	
	mm <- attr(G, "map")
	ped_new <- rbind( attr(G, "ped"), attr(G_add, "ped") )
	G_new <- cbind( unclass(G), G_new )
	attr(G_new, "map") <- mm
	attr(G_new, "ped") <- ped_new
	class(G_new) <- c("genotypes", class(G_new))
	
	return(G_new)
	
}

rename_individuals <- function(G, new_names, ...) {
	
	if (length(new_names) != ncol(G))
		stop("New names must match length of old names.")
	
	if (!length(unique(new_names)) == length(new_names))
		stop("Individual names must be unique.")

	colnames(G) <- as.character(new_names)
	attr(G, "ped")$iid <- new_names
	return(G)
		
}


## find markers in first but not second object
disjoint_markers <- function(G1, G2) {
	setdiff( rownames(G1), rownames(G2) )
}

## find markers that are present in ALL objects provided
## return them in order they appear in FIRST object
common_markers <- function(...) {
	
	ll <- list(...)
	if (length(ll) < 1)
		stop("Must supply at least one genotypes object.")
	
	mm <- lapply(ll, rownames)
	in_all <- Reduce(intersect, mm)
	
	m1 <- rownames(ll[[1]])
	return(m1[ m1 %in% in_all ])
	
}

## merge genotypes objects using markers in common, without checking alleles
merge_genotypes_unsafely <- function(G1, G2) {
	
	if (!.is_same_genotype_encoding(G1,G2))
		stop("Genotypes encoding may not match; cannot merge.")
	
	mm <- common_markers(G1, G2)
	
	G1_new <- G1[mm,]
	G2_new <- G2[mm,]
	
	return( add_individuals(G1_new, G2_new) )
	
}

subset_genotypes <- function(G, expr, ...) {
	
	e <- substitute(expr)
	r <- eval(e, attr(G, "map"), parent.frame())
	r <- r & !is.na(r)
	return( G[ r,,drop=FALSE ] )
	#return( `[.genotypes`(G, r, drop = FALSE) )
	
}

recode_genotypes_to_numeric <- function(G, ...) {
	
	if (!is.character(G))
		stop("Genotypes are not encoded as character.")
	if (!.is_mutable_genotypes(G))
		stop("Genotypes cannot be recoded.")
	
	message("Recoding genotypes as numeric ...")
	
	G <- unclass(G)
	G <- toupper(G)
	
	A1 <- toupper( attr(G,"map")$A1 )
	A2 <- toupper( attr(G,"map")$A2 )
	
	rez <- matrix(NA_integer_, nrow = nrow(G), ncol = ncol(G),
				  dimnames = list(rownames(G), colnames(G)))
	
	for (jj in seq_len(ncol(G))) {
		
		hom_ref <- G[ ,jj ] == A1
		het <- G[ ,jj ] == "H"
		hom_alt <- G[ ,jj ] == A2
		rez[ hom_ref,jj ] <- 0
		rez[ het,jj ] <- 1
		rez[ hom_alt,jj ] <- 2
		
	}
	
	attr(rez,"map") <- attr(G,"map")
	attr(rez,"ped") <- attr(G,"ped")
	class(rez) <- c("genotypes", class(rez))
	
	return(rez)
	
}

.fudge_missing_genotypes <- function(G, ...) {
	
	if (!is.numeric(G))
		stop("Genotypes must be encoded as numeric.")
	
	G <- unclass(G)
	
	mafs <- rowMeans(G, na.rm = TRUE)
	for (ii in seq_len(nrow(G))) {
		nas <- is.na(G[ii,])
		G[ii,nas] <- mafs[ii]
	}
	
	return(G)
	
}

## calculate proportion missing genotypes
p_missing <- function(G, by = c("markers","individuals"), ...) {
	
	G <- unclass(G)
	by_what <- match.arg(by)
	
	nas <- is.na(G)
	
	if (by_what == "markers") {
		return( rowMeans(nas) )
	}
	else if (by_what == "individuals") {
		return( colMeans(nas) )
	}
	
}

## calculate proportion heterozygous genotypes
p_het <- function(G, by = c("markers","individuals"), ...) {
	
	G <- unclass(G)
	by_what <- match.arg(by)
	
	if (is.character(G)) {
		hs <- G == "H"
	}
	else if (is.numeric(G)) {
		hs <- G == 1
	}
	
	if (by_what == "markers") {
		return( rowMeans(hs, na.rm = TRUE) )
	}
	else if (by_what == "individuals") {
		return( colMeans(hs, na.rm = TRUE) )
	}
	
}

## calculate minor allele frequency
calc_maf <- function(G, ...) {
	
	if (!is.numeric(G))
		G <- recode_genotypes_to_numeric(G)
	
	x <- tally_genotype_calls(G)[ ,1:3 ]
	afs <- (x[ ,1 ] + x[ ,2 ])/rowSums(x)
	return( pmin(afs, 1-afs) )
	
}

## do PCA on genotypes
pca_genotypes <- function(G, max_missing = 0.1, K = 10, ...) {
	
	message("Starting with genotype matrix of size ", nrow(G), " markers x ", ncol(G), " individuals.")
	
	prop_miss <- p_missing(G)
	G2 <- G[ prop_miss < max_missing, ]
	message("Dropped ", nrow(G) - nrow(G2), " markers with proportion > ", max_missing, " missing calls.")
	
	message("Imputing missing values ...")
	G2_imp <- .fudge_missing_genotypes(G2)
	
	message("performing PCA on final genotype matrix of size ", nrow(G2), " markers x ", ncol(G2), " individuals ...")
	pc <- prcomp(G2_imp, center = TRUE, scale. = TRUE, rank. = K)
	eigvals <- pc$sdev^2
	varexp <- eigvals <- sum(eigvals)
	rez <- tibble::tibble(iid = colnames(G2), tibble::as_tibble(pc$rotation))
	
	attr(rez,"eigenval") <- eigvals
	attr(rez,"eigvals") <- eigvals
	attr(rez,"explained") <- varexp
	attr(rez,"result") <- pc
	class(rez) <- c("pca_result", class(rez))
	
	message("Done.")
	return(rez)
	
}

drop_intensity <- function(G, ...) {
	
	if (!is.null(attr(G,"intensity"))) {
		message("Has intensity data; deleting ...")
		attr(G,"intensity") <- NULL
	}
	else {
		message("No intensity data; will do nothing.")
	}
		
	return(G)
	
}

## convert raw R/qtl2 genotype matrix to tibble, adding marker map if available
as_tibble.cross2 <- function(x, mm = NULL, ...) {
	
	.geno_to_tibble <- function(g) {
		as_tibble(g, rownames = "iid") %>%
			pivot_longer(., -iid, names_to = "marker")
	}
	
	rez <- lapply(x$geno, .geno_to_tibble) %>% bind_rows()
	if (tibble::is_tibble(mm))
		rez <- left_join(rez, mm)
	return(rez)
	
}

## convert raw R/qtl2 genotype matrix to "wide"-format tibble (markers = rows, individuals = cols)
as_tibble_wide <- function(x, mm = NULL, ...) {

	.geno_to_tibble <- function(g) {
		t(g) %>% as_tibble(rownames = "marker")
	}
	
	rez <- lapply(x$geno, .geno_to_tibble) %>% bind_rows()
	
	if (!is.null(mm)) {
		this_m <- rez[,"marker"]
		this_m <- left_join(this_m, mm)
		rez <- bind_cols( this_m, select(rez, -marker) )
	}
	
	return(rez)
		
}

## convert argyle::genotypes object to "wide"-format tibble (markers = rows, individuals = cols)
as_tibble.genotypes <- function(x, ...) {
	mm <- attr(x, "map")
	rez <- as_tibble(unclass(x), rownames = "marker")
	this_m <- rez[ ,"marker" ]
	this_m <- left_join(this_m, mm)
	rez <- bind_cols( this_m, select(rez, -marker) )
	return(rez)
}

## take R/qtl2 genotpyes, converted to tibble, and plot them like a grid
plot_geno_matrix <- function(df, ...) {
	
	mm <- unique( df[ ,c("chr","marker","pos") ] )
	
	mm$chr <- factor_chrom(mm$chr)
	mm <- arrange(mm, chr, pos)
	print(head(mm))
	df$marker <- factor(df$marker, levels = rev(mm$marker))
	
	
	p <- ggplot(df) +
		geom_tile(aes(x = iid, y = marker, fill = factor(value))) +
		theme(axis.text.x = element_blank()) +
		theme(axis.text.y = element_text(size = 6))
	return(p)
	
}