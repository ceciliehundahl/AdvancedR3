#' Calculate descriptive statistics of each metabolite.
#'
#' @param data The lipidomics dataset.
#'
#' @return A data.frame/tibble.
#'
create_table_descriptive_stats <- function(data) {
  data |>
    dplyr::group_by(metabolite) |>
    dplyr::summarise(dplyr::across(value, list(
      mean = mean,
      sd = sd,
      median = median,
      iqr = IQR
    ))) |>
    dplyr::mutate(dplyr::across(
      tidyselect::where(is.numeric),
      \(x) round(x, digits = 1)
    )) |>
    dplyr::mutate(
      MeanSD = glue::glue("{value_mean} ({value_sd})"),
      medianIQR = glue::glue("{value_median} ({value_iqr})")
    ) |>
    dplyr::select(
      Metabolite = metabolite,
      `Mean (SD)` = MeanSD,
      `Median (IQR)` = medianIQR
    )
}

#' Plot for basic distribution of metabolite data.
#'
#' @param data The lipidomics dataset.
#'
#' @return A plot object.
#'
create_plot_distributions <- function(data) {
  data |>
    ggplot2::ggplot(ggplot2::aes(x = value)) +
    ggplot2::geom_histogram() +
    ggplot2::facet_wrap(ggplot2::vars(metabolite), scales = "free") +
    ggplot2::theme_classic() +
    ggplot2::labs(x = "Concentration", y = "Count")
}

#' Do some cleaning to fix issues in the data.
#'
#' @param data The lipidomics data.
#'
#' @returns A data frame.
#'
clean <- function(data) {
  data |>
    dplyr::group_by(dplyr::pick(-value)) |>
    dplyr::summarise(value = mean(value), .groups = "keep") |>
    dplyr::ungroup()
}

#' Preprocess the data.
#'
#' @param data The lipidomics data.
#'
#' @returns A data frame.
#'
preprocess <- function(data) {
  data |>
    mutate(
      class = as.factor(class),
      value = scale(value)
    )
}

#' Fit the model to the data and get the results.
#'
#' @param data The data to fit.
#' @param model The formula.
#'
#' @returns A data frame of the results.
#'
fit_model <- function(data, model) {
  glm(
    formula = model,
    data = data,
    family = binomial
  ) |>
    broom::tidy(exponentiate = TRUE) |>
    dplyr::mutate(
      metabolite = unique(data$metabolite),
      model = format(model),
      .before = tidyselect::everything()
    )
}

#' Create model results for report.
#'
#' @param data The lipidomics data.
#'
#' @returns A data frame of model results.
#'
create_model_results <- function(data) {
  data |>
    dplyr::filter(metabolite == "Cholesterol") |>
    preprocess() |>
    fit_model(class ~ value)
}
