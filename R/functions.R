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
    recipes::update_role({{ metabolite_variable }}, age, gender, new_role = "predictor") |>
    recipes::update_role(class, new_role = "outcome") |>
    recipes::step_normalize(tidyselect::starts_with("metabolite_"))
}


#' Create workflow object
#'
#' @param model_specs Parsnip model specs
#' @param recipe_specs Recipe specs
#'
#' @return workflows object
create_model_workflow <- function(model_specs, recipe_specs) {
  workflows::workflow() |>
    workflows::add_model(model_specs) |>
    workflows::add_recipe(recipe_specs)
}

#' Tidy model output of workflow fitted model
#'
#' @param workflow_fitted_model workflow fitted model
#'
#' @return tibble
tidy_model_output <- function(workflow_fitted_model) {
  workflow_fitted_model |>
    workflows::extract_fit_parsnip() |>
    broom::tidy(exponentiate = TRUE)
}


#' Split Lipidomics data set to list by metabolites
#'
#' @param data Lipidomics data set long format
#'
#' @return list
split_by_metabolite <- function(data, col) {
  data %>%
    column_values_to_snake_case(metabolite) |>
    dplyr::group_split(metabolite) |>
    purrr::map(metabolite_to_wider)
}

#' Generate model results
#'
#' @param data Lipidomics wide format
#'
#' @return list
generate_model_results <- function(data) {
  create_model_workflow(
    parsnip::logistic_reg() |>
      parsnip::set_engine("glm"),
    data |>
      create_recipe_spec(tidyselect::starts_with("metabolite_"))
  ) |>
    fit(data) |>
    tidy_model_output()
}


#' Add original metabolite names to model results
#'
#' @param model_results Model results
#' @param data Original lipidomics
#'
#' @return tibble
add_original_metabolite_names <- function(model_results, data) {
  data %>%
    dplyr::mutate(term = metabolite) %>%
    column_values_to_snake_case(term) %>%
    dplyr::mutate(term = stringr::str_c("metabolite_", term)) %>%
    dplyr::distinct(term, metabolite) |>
    dplyr::right_join(model_results, by = "term")
}

#' Calculate model estimates
#'
#' @param data lipidomics data
#'
#' @return tibble
calculate_estimates <- function(data) {
  data |>
    split_by_metabolite() |>
    map(generate_model_results) |>
    list_rbind() |>
    filter(str_detect(term, "metabolite_")) |>
    add_original_metabolite_names(data)
}


#' Plotting model estimates
#'
#' @param results tibble of estimates
#'
#' @return ggplot object
plot_estimates <- function(results) {
  results |>
    ggplot(aes(
      x = estimate,
      y = metabolite,
      xmin = estimate - std.error,
      xmax = estimate + std.error
    )) +
    geom_pointrange() +
    coord_fixed(xlim = c(0, 5))
}
