## calculate log-likelihood for data given nu under gamma model from `xoi`
loglik <- function(nu, data, max_conv = 10, integr_tol = 1e-6, max_subd = 100, min_subd = 10) {
	
	d <- data$d/100 # cM --> M
	censor <- data$censor
	len_d <- as.integer(length(d))
	
	result <- .C("GammaS",
				 len_d,
				 as.double(d), 
				 as.integer(censor),
				 1L,
				 nu = as.double(nu), 
				 loglik = as.double(0.0),
				 as.integer(max_conv), 
				 0L,
				 as.double(integr_tol),
				 as.integer(max_subd), 
				 as.integer(min_subd),
				 PACKAGE = "xoi")
	
	return(result$loglik)
	
}

## (log)prior prob at current values of theta
logprior <- function(theta, mu = NULL, sigma = NULL, ...) {
	
	if (is.null(mu))
		mu <- rep(0, length(theta))
	if (is.null(sigma))
		sigma <- diag(length(theta))
	
	mvtnorm::dmvnorm(theta, mean = mu, sigma = sigma, log = TRUE)
	
}

## proposal distribution
## multivariate gaussian, centered on current values of theta
proposal <- function(theta, Sigma, ...) {
	
	#Sigma <- sigma*diag(length(theta))
	mvtnorm::rmvnorm(1, mean = theta, sigma = Sigma)
	
}

#' Set up variables for hierarchical model using formula notation
.prep_data_splits <- function(formula, data, ...) {
	
	M <- model.matrix(formula, data)
	nn <- colnames(M)
	grouper <- apply(M, 1, paste0, collapse = "")
	k <- ncol(M)
	rez <- list( data = split.data.frame(data, grouper),
				 design = split.data.frame(M, grouper),
				 coefs = nn,
				 k = k )
	return(rez)
	
}

#' Calculate log-likelihood in each split (under current and proposed `nu`); sum them up
.splitwise_loglik <- function(obj, theta, ...) {
	
	ll <- rep(0.0, length(obj$data))
	for (ii in seq_along(obj$data)) {
		nu_curr <- exp(obj$design[[ii]][1,] %*% t(theta))
		#message("nu_curr = ", nu_curr, " at index = ", ii)
		ll_this <- loglik(nu_curr, obj$data[[ii]])
		ll[ii] <- ll_this
	}
	lltot <- sum(ll)
	return(lltot)
	
}

# Metropolis-Hastings MCMC function
do_mcmc <- function(data, num_iterations, intercept_mu = 2.5, proposal_sigma = 0.5, prior_sigma = 0.5, null_model = FALSE) {
	
	# prepare data
	fm <- formula(~ cross + age)
	#fm <- formula(~ age)
	#fm <- formula(~ cross)
	data_split <- .prep_data_splits(fm, data)
	k <- data_split$k
	
	# prepare containers for MCMC trace
	chain <- matrix(0.0, ncol = k, nrow = num_iterations) # To store the sampled values
	ll <- numeric(num_iterations) # save log-likelihoods
	
	# define parameter covariance matrix
	# under "null" (intercept-only) model, only the intercept term varies and others are fixed at zero
	# under "full" model, all parameters can vary
	if (null_model) {
		Sigma_prior <- 0*diag(k)
		Sigma_prior[1,1] <- prior_sigma
		Sigma_prop <- 0*diag(k)
		Sigma_prop[1,1] <- proposal_sigma
		message("Fitting NULL MODEL")
	}
	else {
		Sigma_prior <- prior_sigma*diag(k)
		Sigma_prop <- proposal_sigma*diag(k)
	}
	#print(Sigma)

	# draw initial values for theta
	# here we allow the intercept term to have different mean that reflects knowledge that nu ~= 12 ~= exp(2.5)
	mu_prior <- c(intercept_mu, rep(0, k-1))
	initial_value <- mvtnorm::rmvnorm(1, mu_prior, Sigma_prior)
	
	current <- initial_value # Initial state
	ll_curr <- .splitwise_loglik(data_split, current) # current log-likelihood

	message("Prior mean = ", intercept_mu)
	message("Prior covariance = ")
	print(Sigma_prior)
	message("")
	message("Proposal covariance = ")
	print(Sigma_prop)
	message("")

	message("Initial theta = [", paste0(round(current, 2), collapse = " "), "] , loglik = ", ll_curr)
	
	did_accept <- 0
	t_last <- Sys.time()
	ii_last <- 0
	for (ii in 1:num_iterations) {
		
		if (ii %% 10 == 0) {
			
			t_now <- Sys.time()
			t_diff <- difftime(t_last, t_now, units = "secs")/(ii_last-ii)
			message(" ... iter ", ii, " : theta = [", paste0(round(current, 2), collapse = " "), "] , loglik = ", ll_curr)
			message("      ", round(t_diff, 3), " secs per iteration ; cumulative accept = ", did_accept/ii)
			
			ii_last <- ii
			t_last <- t_now
			
		}
		
		# Generate a proposal from the proposal distribution
		proposed <- proposal(current, Sigma_prop)
		
		# log-likelihood of proposed value
		ll_prop <- .splitwise_loglik(data_split, proposed)
		
		# Calculate acceptance probability
		acceptance_prob <- min(1.0, exp((ll_prop + logprior(proposed, mu_prior, Sigma_prior))
										- (ll_curr + logprior(current, mu_prior, Sigma_prior))))
		if (is.na(acceptance_prob))
			acceptance_prob <- 0
		
		# Accept or reject the proposed value
		if (runif(1) < acceptance_prob) {
			chain[ii,] <- proposed
			ll[ii] <- ll_prop
			current <- proposed
			ll_curr <- ll_prop
			did_accept <- did_accept + 1
		} else {
			chain[ii,] <- current
			ll[ii] <- ll_curr
		}
		
	}
	
	message("Acceptance prob = ", did_accept / num_iterations)
	rez <- list(
		intercept_mu = intercept_mu,
		prior_sigma = Sigma_prior,
		proposal_sigma = Sigma_prop,
		trace = chain,
		loglik = ll,
		acceptance_prop = did_accept/num_iterations,
		data = data_split )
	
	return(rez)
	
}
