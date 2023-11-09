#' Descriptive function
#'
#' @param data Liipidomics data set
#'
#' @return tibble with descriptive statistics
#'
descriptive_stats <- function(data) {
  data |>
    dplyr::group_by(metabolite) |>
    dplyr::summarise(dplyr::across(value, list(mean = mean, sd = sd, iqr = IQR))) |>
    dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), ~ round(.x, digits = 1)))
}


#' Simple geom_histogram plot
#'
#' @param data Lipidomics data set
#'
#' @return ggplot object
plot_distributions <- function(data) {
  data |> ggplot2::ggplot(ggplot2::aes(x = value)) +
    ggplot2::geom_histogram() +
    ggplot2::facet_wrap(vars(metabolite), scales = "free")
}


#' Convert columns to snakecase
#'
#' @param data Data set
#' @param cols Vector of columns. Tidyverse compatible
#'
#' @return tibble
column_values_to_snake_case <- function(data, cols) {
  data |> dplyr::mutate(across({{ cols }}, snakecase::to_snake_case))
}

#' Pivot lipidomics metabolite to wider
#'
#' @param data Lipidomics data set
#'
#' @return tibble
metabolite_to_wider <- function(data) {
  data |>
    column_values_to_snake_case(metabolite) |>
    tidyr::pivot_wider(
      names_from = metabolite,
      values_from = value,
      values_fn = mean,
      names_prefix = "metabolite_"
    )
}

#' Specifying recipes
#'
#' @param data Lipidomics data set
#' @param metabolite_variable matabolite variable
#'
#' @return recipe
create_recipe_spec <- function(data, metabolite_variable) {
  recipes::recipe(data) %>%
    recipes::update_role({{ metabolite_variable }}, age, gender, new_role = "predictor") %>%
    recipes::update_role(class, new_role = "outcome") %>%
    recipes::step_normalize(tidyselect::starts_with("metabolite_"))
}
