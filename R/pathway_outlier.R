#' Outlier-Adjusted Causal Forest Analysis
#'
#' This function performs a causal forest analysis after adjusting for outliers in
#' the treatment and outcome variables. It identifies the observation with the
#' highest estimated Conditional Average Treatment Effect (CATE) among the non-outliers.
#'
#' The function first excludes outliers in the treatment and outcome variables
#' based on the specified threshold. It then performs a causal forest analysis
#' using the 'grf' package on the remaining data to estimate CATEs.
#'
#' @param data A data frame containing the variables for analysis, including
#' treatment and outcome variables and other covariates.
#' @param treatment The name of the treatment variable in the data frame.
#' @param outcome The name of the outcome variable in the data frame.
#' @param outlier_threshold A numeric value representing the threshold for outlier
#' detection. Observations with treatment or outcome values beyond this threshold
#' will be considered outliers and excluded from the analysis.
#' @param num.trees Number of trees to grow in the causal forest model.
#' Defaults to 2000.
#' @param min.node.size Minimum size of the nodes in the causal forest.
#' Defaults to 1.
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
#' @return A character string describing the observation number with the highest CATE
#' among non-outliers, the estimated CATE value for this observation, and a statement
#' that this observation is among the non-outliers. The index of the observation is
#' with respect to the original dataset.
#'
#' @export
#' @importFrom grf causal_forest
#' @importFrom stats quantile
pathway_outlier <- function(data, treatment, outcome, outlier_threshold,
                            num.trees = 2000, min.node.size = 1,
                            honesty = TRUE, honesty.fraction = 0.5,
                            seed = NULL, ci.group.size = 2,
                            alpha = 0.05, imbalance.penalty = 0,
                            clusters = NULL, sample.fraction = 0.5,
                            mtry = floor(sqrt(ncol(data) - 2))) {

  # Ensure that treatment and outcome are in the data
  if (!(treatment %in% names(data))) {
    stop("Treatment variable not found in the data.")
  }
  if (!(outcome %in% names(data))) {
    stop("Outcome variable not found in the data.")
  }

  # Define the lower and upper bounds for outlier detection
  lower_bound <- quantile(data[[treatment]], probs = outlier_threshold / 2, na.rm = TRUE)
  upper_bound <- quantile(data[[treatment]], probs = 1 - outlier_threshold / 2, na.rm = TRUE)

  # Identify non-outlier observations based on treatment
  non_outlier_indices_treatment <- which(data[[treatment]] >= lower_bound & data[[treatment]] <= upper_bound)

  # Repeat for outcome variable
  lower_bound <- quantile(data[[outcome]], probs = outlier_threshold / 2, na.rm = TRUE)
  upper_bound <- quantile(data[[outcome]], probs = 1 - outlier_threshold / 2, na.rm = TRUE)

  # Identify non-outlier observations based on outcome
  non_outlier_indices_outcome <- which(data[[outcome]] >= lower_bound & data[[outcome]] <= upper_bound)

  # Combine both indices
  non_outlier_indices <- intersect(non_outlier_indices_treatment, non_outlier_indices_outcome)

  # Prepare the data for causal forest
  X <- data[non_outlier_indices, !(names(data) %in% c(treatment, outcome))] # Covariates
  Y <- data[non_outlier_indices, outcome] # Outcome variable
  W <- data[non_outlier_indices, treatment] # Treatment variable

  # Estimate the causal forest
  cforest <- causal_forest(X, Y, W, num.trees = num.trees, min.node.size = min.node.size,
                           honesty = honesty, honesty.fraction = honesty.fraction,
                           ci.group.size = ci.group.size, alpha = alpha,
                           imbalance.penalty = imbalance.penalty,
                           clusters = clusters, sample.fraction = sample.fraction,
                           mtry = mtry)

  # Estimate CATEs
  cate_estimates <- predict(cforest)$predictions

  # Find the observation with the highest CATE among the non-outliers
  highest_cate_index <- which.max(cate_estimates)
  highest_cate <- cate_estimates[highest_cate_index]

  # Return a descriptive statement
  return(paste("The case with the highest CATE among the non-outliers (excluding extremes in treatment and outcome) is observation number",
               non_outlier_indices[highest_cate_index], "in the original dataset, which has the CATE =", highest_cate))
}


