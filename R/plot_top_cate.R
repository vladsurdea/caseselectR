#' Plot Top N Cases Based on CATE
#'
#' This function plots the top N cases based on Conditional Average Treatment Effects (CATE)
#' estimated using a causal forest model. It shows the CATE values along with their
#' 95% confidence intervals for the selected top cases.
#'
#' The function first fits a causal forest model to the data, then calculates the CATE
#' and its confidence intervals for each observation. It selects the top N cases based on
#' the highest CATE values and generates a plot displaying these values and their
#' confidence intervals.
#'
#' @param data A data frame containing the variables for analysis, including
#' the treatment and outcome variables and other covariates.
#' @param treatment The name of the treatment variable in the data frame.
#' @param outcome The name of the outcome variable in the data frame.
#' @param N The number of top cases to display based on CATE. Defaults to 10.
#' @param num.trees Number of trees to grow in the causal forest model.
#' Defaults to 2000.
#' @param alpha Significance level for confidence interval calculation.
#' Defaults to 0.05.
#'
#' @return A ggplot object displaying the top N cases based on CATE, each represented by a point
#' with a corresponding confidence interval. Points are ordered by the magnitude of the CATE.
#'
#' @examples
#' # Assuming 'my_data' is a data frame with variables 'treatment_var' and 'outcome_var'
#' plot_result <- plot_top_cate(my_data, "treatment_var", "outcome_var")
#' print(plot_result)
#'
#' @export
#' @importFrom grf causal_forest
#' @importFrom ggplot2 ggplot aes geom_point geom_errorbar geom_label scale_x_discrete labs
#' @importFrom ggpubr theme_pubr
#' @importFrom stats qnorm
plot_top_cate <- function(data, treatment, outcome, N = 10, num.trees = 2000, alpha = 0.05) {
  # Ensure that treatment and outcome are in the data
  if (!(treatment %in% names(data))) {
    stop("Treatment variable not found in the data.")
  }
  if (!(outcome %in% names(data))) {
    stop("Outcome variable not found in the data.")
  }

  # Prepare data for causal forest
  Y <- data[[outcome]]
  W <- data[[treatment]]
  X <- data[, !(names(data) %in% c(treatment, outcome))]

  # Estimate the causal forest and CATEs
  cforest <- causal_forest(X, Y, W, num.trees = num.trees)
  cate_estimates <- predict(cforest, estimate.variance = TRUE)
  cate_values <- cate_estimates$predictions
  cate_se <- sqrt(cate_estimates$variance.estimates)

  # Calculate confidence intervals
  ci_lower <- cate_values - qnorm(1 - alpha / 2) * cate_se
  ci_upper <- cate_values + qnorm(1 - alpha / 2) * cate_se

  # Data frame with CATEs, CIs, and observation numbers
  top_cate <- data.frame(Observation = 1:length(cate_values),
                         CATE = cate_values,
                         CI_lower = ci_lower,
                         CI_upper = ci_upper)

  # Select top N cases
  top_cate <- top_cate[order(-top_cate$CATE), ][1:N, ]

  # Adding an offset to position the label above the CI_upper
  label_position <- top_cate$CI_upper + max(top_cate$CI_upper - top_cate$CI_lower) * 0.1

  # Create the plot
  ggplot(top_cate, aes(x = reorder(Observation, -CATE), y = CATE)) +
    geom_point(shape = 17, size = 3.5) +  # Layer with empty triangles
    geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.1) +
    geom_label(aes(label = paste("CATE=", round(CATE, 2)), y = label_position),
               vjust = 0,
               fill = "white",
               color = "black") +
    scale_x_discrete(labels = top_cate$Observation) +
    labs(title = paste("Top", N, "Cases Based on CATE"),
         x = "Case",
         y = "CATE with 95% Confidence Interval") +
    theme_pubr() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

