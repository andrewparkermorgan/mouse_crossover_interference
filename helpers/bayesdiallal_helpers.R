merge_chains <- function(chains, burnin = 0, ...) {
	niter <- min(sapply(chains, nrow))
	rez <- lapply(chains, function(f) f[ (burnin+1):niter,, drop = FALSE ] )
	coda::mcmc( do.call(rbind, rez) )
}

get_chains <- function(bd, ...) {
	merge_chains(bd$AllDiallelObs[[1]]$cent.chains, bd$burnin)
}

diallel_hpdi <- function(bd, prob = c(0.5, 0.95), chains = NULL, ...) {
	
	burnin <- bd$burnin
	if (!is.null(chains) && inherits(bd, "mcmc")) {
		# use pre-computed chains (ie. to accomodate post-computed variables)
		chains <- window(chains, start = burnin+1)
	}
	else {
		chains <- merge_chains(bd$AllDiallelObs[[1]]$cent.chains, burnin)
	}
	hpd <- lapply(prob, function(p) cbind(HPDinterval(chains, prob = p), summary(chains, 0.5)[[2]], p))
	rez <- lapply(hpd, rename_effects, strains = bd$strain.map[,1]) %>% bind_rows()
	rez$prob_class <- factor(rez$prob, sort(unique(stats::na.omit(rez$prob))))
	return(rez)
	
}

rename_effect_cols <- function(thenames, strains, ...) {
	
	n <- thenames
	#print(n)
	d <- grepl("^dominancej", n)
	a <- grepl("^aj", n)
	m <- grepl("^motherj", n)
	v <- grepl("^SymCrossjk", n)
	w <- grepl("^ASymCrossjk", n)
	i <- grepl("BetaInbred", n)
	ga <- grepl("^Gender\\:aj", n)
	gm <- grepl("^Gender\\:motherj", n)
	fe <- grepl("^Fixed", n)
	re <- grepl("^RandomEffect", n)
	
	df <- dplyr::tibble(original_col = n,
						group = character(length(n)),
						strain = character(length(n)),
						strain2 = character(length(n)))
	
	df$group[i] <- "inbreeding"
	df$strain[i] <- "all"
	df$group[d] <- "dominance"
	df$strain[d] <- strains
	df$group[a] <- "additive"
	df$strain[a] <- strains
	df$group[m] <- "maternal"
	df$strain[m] <- strains
	df$group[ga] <- "sex x additive"
	df$strain[ga] <- strains
	df$group[gm] <- "sex x maternal"
	df$strain[gm] <- strains
	df$group[fe] <- "fixed effects"
	df$group[re] <- "random effects"
	df$strain[fe] <- "all"
	
	cbo <- combn(length(strains), 2)
	mat <- strains[cbo[2,]]
	pat <- strains[cbo[1,]]
	df$group[w] <- "cross-specific"
	df$strain[w] <- mat
	df$strain2[w] <- pat
	df$group[v] <- "cross-specific"
	df$strain[v] <- pat
	df$strain2[v] <- mat
	df$term[fe] <- gsub("^FixedEffect\\:[0-9]+\\:*", "", n[fe])
	
	return(df)
	
}
	
rename_effects <- function(x, strains, ...) {
	
	n <- rownames(x)
	df <- rename_effect_cols(n, strains, ...)
	
	df$lower <- x[,1]
	df$upper <- x[,2]
	df$estimate <- x[,3]
	#df$prob <- attr(x, "Probability")
	df$prob <- x[,4]
	
	#df$term <- with(df, paste0(group, " :: ", strain))
	#df$term[ is.na(df$strain) ] <- df$group[ is.na(df$strain) ]
	#df$term[ df$strain2 != "" ] <- with(df, paste0(group, " :: (", strain, ";", strain2, ")"))
	df <- subset(df, group != "")
	return(df)
	
}

factor_diallel_strains <- function(x, bd, ...) {
	y <- bd$strain.map[ x,1 ]
	factor(y, bd$strain.map[,1])
}

get_observed_values <- function(bd, ...) {
	
	x <- tibble::as_tibble(bd$Listjk)
	colnames(x) <- c("mat","pat")
	x$mat <- factor_diallel_strains(x$mat, bd)
	x$pat <- factor_diallel_strains(x$pat, bd)
	x$obs <- bd$Y
	return(x)
	
}


get_hsq <- function(bd, method = c("hsq","dsq"), use_fixed_effects = FALSE, ...) {
	
	method <- match.arg(method)
	if (method == "hsq") {
		#rez <- bd$AllDiallelObs[[1]]$PosteriorHSq()
		rez <- BayesDiallel:::PosteriorHSq.DiallelOb(bd$AllDiallelObs[[1]], AFD = bd, UseFixed = use_fixed_effects)
	}
	else if (method == "dsq") {
		#rez <- bd$AllDiallelObs[[1]]$PosteriorDSq()
	}
	else {
		stop("Must choose Hsq or Dsq method for decomposing variance.")
	}
	return( rez[[1]] %>% merge_chains(burnin = bd$burnin) )
	
}

clean_hsq_terms <- function(f, flip_parents = FALSE, ...) {
	
	replacer <- c("aj" = "additive",
				  "dominancej" = "dominance",
				  "BetaInbred:Av" = "inbreeding",
				  "motherj" = if (flip_parents) "maternal" else "paternal",
				  "SymCrossjk" = "cross-specific",
				  "AsymCrossjk" = "cross-specific (asymmetric)",
				  "FixedEffect:1" = "fixed effect")
	last <- c("total.explained" = "total explained",
			  "Noise" = "residual")
	effterms <- f
	effterms <- c(replacer, last)[ effterms ]
	effterms[ is.na(effterms) ] <- f[ is.na(effterms) ]
	f <- effterms
	
	levs <- unique(f)
	newlevs <- setdiff(levs, c(replacer, last))
	f <- factor(f, levels = c(replacer, newlevs, last))
	
	return(f)
	
}

tidy_hsq <- function(x, ...) {
	
	x %>%
		tidybayes::gather_draws(`^[a-zA-Z].+`, regex = TRUE) %>%
		ungroup() %>%
		mutate(term = clean_hsq_terms(.variable, ...)) %>%
		rename(value = .value)
	
}

tidy_mcmc <- function(x, prob = 0.95, ...) {
	
	xx <- summary(x)
	ss <- xx$statistics 
	if (is.null(dim(ss))) {
		## mcmc chain only had one variable
		rez <- as.matrix(ss) %>% t()
		rownames(rez) <- colnames(x)
		rez <- as_tibble(rez, rownames = "term")
	}
	else {
		## mcmc chain had >1 variable; a matrix returns, no issue
		rez <- as_tibble(ss, rownames = "term")
	}
	
	colnames(rez)[2] <- "estimate"
	hpd <- coda::HPDinterval(x, prob = prob) %>% as_tibble(rownames = "term")
	return( left_join( rez[,c("term","estimate")], hpd ) )
	
}

summarise_hsq <- function(bd, prob = 0.95, ...) {
	
	hsq <- get_hsq(bd)
	rez <- tidy_mcmc(hsq)
	#print(rez)
	rez$term <- clean_hsq_terms(rez$term)
	return(rez)
	
}

summarise_observed <- function(bd, ...) {
	
	x <- get_observed_values(bd)
	x %>% group_by(mat, pat) %>%
		summarise(mu_obs = mean(obs, na.rm = TRUE),
				  med_obs = median(obs, na.rm = TRUE),
				  se = sd(obs, na.rm = TRUE)/sqrt(sum(!is.na(obs))))
	
}

predict_posteriors <- function(bd, ...) {
	
	## to get posterior samples rather than summaries, do: getPosteriorPredictionsMeanChains()
	## here we use the summaries straight away
	#post <- tibble::as_tibble(bd$AllDiallelObs[[1]]$getPosteriorPredictionsMeanSummary())
	post <- BayesDiallel:::PosteriorPredSummary.DiallelOb(bd$AllDiallelObs[[1]], AFD = bd)
	colnames(post) <- c("mat","pat","mu","med","lo","hi")
	strains <- bd$strain.map
	post$mat <- factor_diallel_strains(post$mat, bd)
	post$pat <- factor_diallel_strains(post$pat, bd)
	post$term <- attr(post, "row.names")
	return( as_tibble(post) )
	
}

tidy_chains <- function(bd, strains, ...) {
	
	todrop <- c( colnames(bd$RandomEffects) )
	chains <- bd$AllDiallelObs[[1]]$cent.chains
	tokeep <- setdiff(colnames(chains[[1]]), todrop)
	chains <- chains[ ,tokeep ]
	names(chains) <- LETTERS[ 1:length(chains) ]
	
	rez <- lapply(chains, as_tibble, rownames = ".iter") %>%
		bind_rows(.id = ".chain")
	nn <- rename_effect_cols(colnames(rez), strains)
	rez <- pivot_longer(rez, -c(.iter, .chain)) %>% rename(original_col = name)
	rez <- left_join(rez, nn)
	
	return(rez)
	
}

all_pairs <- function(x, by, ...) {
	
	combiner <- function(df) {
		inner_join( subset(x, keycol == df$x[1]),
					subset(x, keycol == df$y[1]),
					by = c("keycol", quote(...)) ) %>%
			select(-keycol)
	}
	
	vals <- select(x, quo(by))[[1]]
	if (is.factor(vals))
		levs <- levels(vals)
	else
		levs <- unique(vals)
	x$keycol <- vals
	pp <- tidyr::crossing(x = levs, y = levs)
	print(pp)
	group_by(pp, x, y) %>%
		do(combiner(.))
	
}


posterior_effect_contrast <- function(bd, group = c("additive","dominance","maternal"), focal = NULL, others = NULL, ...) {

	# get MCMC samples, clipping off the burn-in iterations
	chains <- get_chains(bd)
	
	# get indices of 'focal' and 'background' strains, if specified
	strains <- bd$strain.map[,1]
	idx <- as.integer(bd$strain.map[,2])
	if (is.null(focal)) {
		focal <- 1
	}
	else {
		focal <- match(focal, strains)
	}
	
	if (is.null(others) || "all" %in% others) {
		others <- setdiff(idx, focal)
	}
	else {
		others <- match(others, strains)
	}
	
	if (any(is.na(others), is.na(focal)))
		return(NULL)
	
	# figure out what group of parameters we're after, in order to retrieve them from the sample
	group <- match.arg(group)
	s <- c("additive" = "a", "maternal" = "mother", "dominance" = "dominance")[group]
	# helper function to construct column names for indexing into sample matrix
	.make_colnames <- function(ss, ff) {
		template <- "${effect}j:${focal} - mean(${effect}j)"
		sapply(ff, function(f) stringr::str_interp(template, list(effect = ss, focal = f)))
	}
	
	# construct the contrast vector
	s.focal <- .make_colnames(s, focal)
	s.other <- .make_colnames(s, others)
	if (!all(c(s.focal, s.other) %in% colnames(chains)))
		return(NULL)
	
	nf <- length(s.focal)
	no <- length(s.other)
	beta <- c( rep(1/nf, nf) , rep(-1/no, no) )
	
	# compute the actual contrast
	X <- chains[ ,unname(c(s.focal, s.other)), drop = FALSE ]
	Y <- as.mcmc(X %*% beta)
	colnames(Y) <- "contrast"
	
	# done
	return(Y)
	
}

plot_diallel <- function(x, y = NULL, circles = FALSE, labels = FALSE, flip_parents = TRUE, strain_order = NULL, ...) {

	if (inherits(x, c("FullDiallelAnalyze","FullDiallel"))) {
		pred <- predict_posteriors(x)
		obs <- get_observed_values(x)
	}
	else {
		pred <- x
		obs <- y
	}
	
	obs2 <- group_by(obs, mat, pat) %>%
		summarise(mu = mean(obs, na.rm = TRUE), sigma = sd(obs, na.rm = TRUE)/sqrt(sum(!is.na(obs))))
	obs2 <- ungroup(obs2)
	
	labs <- if (!flip_parents) c("maternal","paternal") else c("paternal","maternal")
	strains <- levels(pred$mat)
	
	toplot <- bind_rows(
		mutate(pred, what = "model posterior"),
		mutate(tidyr::complete(obs2, mat, pat), what = "observed")
	)
	
	toplot$what <- relevel(factor(toplot$what), "observed")
	toplot$pat <- factor(toplot$pat, levels = rev(levels(toplot$mat)))
	toplot$prettylab <- as.character(round(toplot$mu, 1))
	toplot$prettylab[ is.na(toplot$prettylab) ] <- "-"
	toplot <- unique(toplot)
	
	if (!is.null(strain_order)) {
		toplot$mat <- factor(toplot$mat, strain_order)
		toplot$pat <- factor(toplot$pat, levels = rev(strain_order))
	}
	
	if (circles) {
		
		ns <- nlevels(toplot$mat)
		toplot$x0 <- as.numeric(toplot$mat)
		toplot$y0 <- as.numeric(toplot$pat)
		maxr <- max(toplot$mu, na.rm = TRUE)
		minr <- min(toplot$mu, na.rm = TRUE)
		toplot$r <- (toplot$mu - minr)/(2*(maxr-minr))
		toplot$r[ is.na(toplot$mu) ] <- 0.1
		
		p <- ggplot(toplot) +
			ggforce::geom_circle(aes(x0 = x0, y0 = y0, r = r, fill = r), colour = "grey90") +
			scale_x_continuous(labs[1], breaks = seq(ns), labels = levels(toplot$mat)) +
			scale_y_continuous(labs[2], breaks = seq(ns), labels = levels(toplot$pat)) +
			facet_grid(. ~ what)
		
	}
	else {
		
		p <- ggplot(toplot) +
			geom_tile(aes(x = mat, y = pat, fill = mu)) +
			scale_x_discrete(labs[1]) +
			scale_y_discrete(labs[2]) +
			facet_grid(. ~ what)
		
	}
	
	if (labels) {
		p <- p +
			geom_text(aes(x = mat, y = pat, label = prettylab, colour = (prettylab == "-"))) +
			scale_color_manual(values = c("white","grey80"), guide = FALSE)
	}
	
	return(p)
	
}

estimate_parent_effect <- function(fitted, focal = "G", others = c("all","F","H","Z","R"), ...) {
	tmp <- tidyr::crossing("effect" = c("additive","maternal","dominance"), "focal" = focal, "comparing" = others)
	get_contrast <- function(e, o) {
		#print(list(e, o))
		eff <- posterior_effect_contrast(fitted, e, focal, o)
		if (!is.null(eff)) {
			rez <- tidy_mcmc(eff)
			rez$effect <- e
			rez$focal <- paste0(focal, collapse = ",")
			rez$comparator <- paste0(o, collapse = ",")
		}
		else {
			rez <- tibble()
		}
		return(rez)
	}
	betas <- vector("list", length = nrow(tmp))
	for (ii in seq_len(nrow(tmp))) {
		betas[[ii]] <- get_contrast(unlist(tmp$effect[ii]), unlist(tmp$comparing[ii]))
	}
	betas <- bind_rows(betas)
	return(betas)
}