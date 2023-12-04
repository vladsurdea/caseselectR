#' Time Period Analysis of CATE with Causal Forest
#'
#' This function performs a causal forest analysis to estimate Conditional Average Treatment Effects
#' (CATE) over different time periods for each unique identifier in the dataset. It returns a plot
#' showing the evolution of CATEs over time for the ID with the highest overall CATE.
#'
#' The function estimates CATEs separately for each unique identifier in the data, considering
#' the specified time period. It then identifies the identifier with the highest overall CATE and
#' generates a plot of CATEs across different time periods for this identifier.
#'
#' @param data A data frame containing the variables for analysis, including identifiers, time periods,
#' treatment and outcome variables, and other covariates.
#' @param id The name of the identifier variable in the data frame.
#' @param time_period The name of the time period variable in the data frame.
#' @param treatment The name of the treatment variable in the data frame.
#' @param outcome The name of the outcome variable in the data frame.
#' @param num.trees Number of trees to grow in the causal forest model. Defaults to 2000.
#'
#' @return A list containing two elements: 'plot', a ggplot object depicting the evolution of CATEs over
#' time for the ID with the highest overall CATE, and 'pred', a vector of the predicted CATE values for
#' this ID.
#'
#' @examples
#' # Assuming 'my_data' is a data frame with variables 'ID', 'Time_Period', 'treatment_var', and 'outcome_var'
#' result <- pathway_time(my_data, "ID", "Time_Period", "treatment_var", "outcome_var")
#' print(result$plot)
#'
#' @export
#' @importFrom grf causal_forest
#' @importFrom ggplot2 ggplot aes geom_line geom_point labs
#' @importFrom ggpubr theme_pubr
pathway_time <- function(data, id, time_period, treatment, outcome, num.trees = 2000) {
  # Ensure required columns are in the data
  for (col in c(id, time_period, treatment, outcome)) {
    if (!(col %in% names(data))) {
      stop(paste("Column", col, "not found in the data."))
    }
  }

  # Initialize a data frame to store CATEs for each time period and id
  cate_data <- data.frame(ID = integer(), Time_Period = integer(), CATE = numeric())

  # Calculate CATEs for each id
  for (unique_id in unique(data[[id]])) {
    id_data <- subset(data, data[[id]] == unique_id)
    Y <- id_data[[outcome]]
    W <- id_data[[treatment]]
    X <- id_data[, !(names(id_data) %in% c(id, treatment, outcome))]

    # Estimate the causal forest and predict CATEs
    cforest <- causal_forest(X, Y, W, num.trees = num.trees)
    cate_estimates <- predict(cforest)$predictions

    # Combine with time periods
    cate_data <- rbind(cate_data, data.frame(ID = unique_id, Time_Period = id_data[[time_period]], CATE = cate_estimates))
  }

  # Find the ID with the highest overall CATE
  max_cate_id <- cate_data$ID[which.max(cate_data$CATE)]

  # Filter data for this ID
  plot_data <- cate_data[cate_data$ID == max_cate_id, ]

  # Create the plot
  plot <- ggplot(plot_data, aes(x = Time_Period, y = CATE)) +
    geom_line(linetype="dashed") +
    geom_point(shape = 17, size = 3.5) +  # Layer with empty triangles
    labs(title = paste("Evolution of CATEs Over Time for ID", max_cate_id),
         x = "Time Period",
         y = "CATE") +
    theme_pubr()

  # Print the plot
  print(plot)

  # Return the predictions as a vector
  pred <- plot_data$CATE
  return(list(plot = plot, pred = pred))
}

