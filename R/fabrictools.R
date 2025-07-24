library(phytools)
library(MASS)
library(dplyr)

calc_pred_beta <- function(tree, sig2) {
  branch_lengths <- sapply(tree$edge.length, function(x) x)[-1]
  pred_beta <- 2 * sqrt(3 * sig2 * branch_lengths)
  return(list(branch_lengths = branch_lengths, pred_beta = pred_beta))
}

test_zero_branch_length <- function(tree) {
  for (node in tree$edge) {
    if (node != tree$root && tree$edge.length[node] == 0) {
      cat("tree has a zero length branch for node ", node, "\n")
      quit(status = 1)
    }
  }
}

get_sig2 <- function(trees, argv2) {
  if (is.numeric(argv2)) {
    return(as.numeric(argv2))
  } else {
    data <- read.csv(argv2)
    bind_data(trees, data)
    
    model <- contrasts(trees)
    ML <- fit_all_maximum_likelihood(model)
    
    return(ML$Sigma2[1])
  }
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  cat(paste(args[1], "takes a tree file and a data file.\n"))
  quit(status = 1)
}

trees <- read.nexus(args[1])
test_zero_branch_length(trees[[1]])

sig2 <- get_sig2(trees, args[2])

result <- calc_pred_beta(trees[[1]], sig2)
branch_lengths <- result$branch_lengths
pred_beta <- result$pred_beta

df <- data.frame(branch_length = branch_lengths, pred_beta = pred_beta)
write.csv(df, paste0(args[1], ".csv"), row.names = FALSE)

cat("sigma2:", sig2, "\n")

header <- c("distribution", "lh", "parameter 1", "parameter 2")
cat(header, sep = "\t", "\n")

# parameter order is not consistent. They are the order BayesTraits prior settings expect them. 

fit_gamma <- fitdistr(pred_beta, "gamma", start = list(shape = 1, rate = 1))
lh_gamma <- dgamma(pred_beta, shape = fit_gamma$estimate[1], rate = fit_gamma$estimate[2])
cat("gamma", sum(log(lh_gamma)), fit_gamma$estimate[1], fit_gamma$estimate[2], sep = '\t', "\n")

fit_weibull <- fitdistr(pred_beta, "weibull", start = list(shape = 1, scale = 1))
lh_weibull <- dweibull(pred_beta, shape = fit_weibull$estimate[1], scale = fit_weibull$estimate[2])
cat("weibull", sum(log(lh_weibull)), fit_weibull$estimate[2], fit_weibull$estimate[1], sep = '\t', "\n")

fit_lognorm <- fitdistr(pred_beta, "lognormal", start = list(meanlog = 0, sdlog = 1))
lh_lognorm <- dlnorm(pred_beta, meanlog = fit_lognorm$estimate[1], sdlog = fit_lognorm$estimate[2])
cat("lognorm", sum(log(lh_lognorm)), fit_lognorm$estimate[2], fit_lognorm$estimate[1], sep = '\t', "\n")