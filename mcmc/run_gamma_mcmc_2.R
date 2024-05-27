#! /usr/bin/env Rscript

library(optparse)
library(tidyverse)
library(xoi)

source("gamma_mcmc_2.R")

opts <- list(
	make_option("--prop_sigma", type="double", default = 0.5, 
				help="std dev of proposal distribution [default %default]",
				metavar="proposal std dev"),
	make_option("--prior_sigma", type="double", default = 0.5, 
				help="std dev by which coovariance matrix of prior distribution is scaled [default %default]",
				metavar="prior std dev"),
	make_option("--intercept_mu", type="double", default = 2.5, 
				help="prior mean for intercept term [default %default]",
				metavar="prior mean for intercept only"),
	make_option("--null_model", action="store_true", default=FALSE,
		    		help="fit null (intercept-only) model"),
	make_option(c("-n","--niter"), type="integer", default = 10000, 
				help="number of MCMC iterations [default %default]",
				metavar="num iterations"),
	make_option(c("-o","--output"), type="character", default = "out.rds", 
				  help="name of output file [default %default]",
				  metavar="output file"),
	make_option(c("-i","--input"), type="character", default = "input.txt", 
				help="input file of interval widths (cM) and censoring codes [default %default]",
				metavar="input file")
	)

## parse command-line args
args <- parse_args( OptionParser(option_list = opts) )

## read input data
message("Reading input data from ", args$input, " ...")
#data <- read_table(args$input, col_names = FALSE)
#colnames(data) <- c("d","censor")
data <- readRDS(args$input)

## run the MCMC
result <- do_mcmc(data, args$niter,
				  intercept_mu = args$intercept_mu,
				  proposal_sigma = args$prop_sigma,
				  prior_sigma = args$prior_sigma,
				  null_model = args$null_model)

result$intercept_mu <- args$intercept_mu
result$proposal_sigma <- args$prop_sigma
result$prior_sigma <- args$prior_sigma

saveRDS(result, args$output)
