#' Selecting pathway case with causal forests
#'
#' This function performs a causal forest analysis using the 'grf' package to estimate
#' the Conditional Average Treatment Effects (CATE) on the provided data. It identifies
#' the observation with the highest estimated CATE, and returns its value along with a
#' confidence interval.
#'
#' @param data A data frame containing the variables for analysis. It must include
#' both the treatment and outcome variables, as well as any other covariates.
#' @param treatment The name of the treatment variable in the data frame.
#' @param outcome The name of the outcome variable in the data frame.
#' @param num.trees Number of trees to grow in the causal forest model. Defaults to 2000.
#' @param min.node.size Minimum size of the nodes in the causal forest. Defaults to 1.
#' @param honesty A logical value indicating whether to use honest estimation.
#' Defaults to TRUE.
#' @param honesty.fraction Fraction of the data to be used for estimation of
#' treatment effects in honest estimation. Defaults to 0.5.
#' @param seed An optional integer value for random seed to ensure reproducibility.
#' If NULL, the seed is not set. Defaults to NULL.
#' @param ci.group.size Size of the groups for confidence interval estimation.
#' Defaults to 2.
#' @param alpha Significance level for confidence interval calculation.
#' Defaults to 0.05.
#' @param imbalance.penalty A penalty term for imbalance in the splitting rule.
#' Defaults to 0.
#' @param clusters An optional vector of cluster IDs to be used for clustered sampling.
#' If NULL, no clustering is used. Defaults to NULL.
#' @param sample.fraction Fraction of the data to be used for building each tree.
#' Defaults to 0.5.
#' @param mtry Number of variables randomly sampled as candidates at each split.
#' Defaults to the floor of the square root of the number of variables in 'data' minus 2.
#'
#' @return A character string describing the observation number with the highest CATE,
#' the estimated CATE value for this observation, and its confidence interval.
#'
#' @examples
#' # Assuming 'my_data' is a data frame with variables 'treatment_var' and 'outcome_var'
#' result <- pathway_case(my_data, "treatment_var", "outcome_var")
#' print(result)
#'
#' @export
#' @import grf
pathway_case <- function(data, treatment, outcome, num.trees = 2000, min.node.size = 1,
                         honesty = TRUE, honesty.fraction = 0.5, seed = NULL,
                         ci.group.size = 2, alpha = 0.05, imbalance.penalty = 0,
                         clusters = NULL, sample.fraction = 0.5, mtry = floor(sqrt(ncol(data) - 2))) {
  if (!(treatment %in% names(data))) {
    stop("Treatment variable not found in the data.")
  }
  if (!(outcome %in% names(data))) {
    stop("Outcome variable not found in the data.")
  }
  Y <- data[[outcome]]
  W <- data[[treatment]]
  X <- data[, !(names(data) %in% c(treatment, outcome))]
  if (!is.null(seed)) {
    set.seed(seed)
  }
  cforest <- causal_forest(X, Y, W, num.trees = num.trees, min.node.size = min.node.size,
                           honesty = honesty, honesty.fraction = honesty.fraction,
                           ci.group.size = ci.group.size, alpha = alpha,
                           imbalance.penalty = imbalance.penalty,
                           clusters = clusters, sample.fraction = sample.fraction,
                           mtry = mtry)
  cate_estimates <- predict(cforest)$predictions
  ci <- predict(cforest, estimate.variance = TRUE)$variance.estimates
  ci_lower <- cate_estimates - qnorm(1 - alpha / 2) * sqrt(ci)
  ci_upper <- cate_estimates + qnorm(1 - alpha / 2) * sqrt(ci)

  highest_cate_index <- which.max(cate_estimates)
  highest_cate <- cate_estimates[highest_cate_index]
  highest_cate_ci_lower <- ci_lower[highest_cate_index]
  highest_cate_ci_upper <- ci_upper[highest_cate_index]

  return(paste("The case with the highest CATE is observation number", highest_cate_index,
               "which has the CATE =", highest_cate,
               "with a confidence interval of [", highest_cate_ci_lower, ",", highest_cate_ci_upper, "]"))
}






