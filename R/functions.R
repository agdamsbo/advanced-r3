#' Descriptive function
#'
#' @param data Liipidomics data set
#'
#' @return tibble with descriptive statistics
#'
descriptive_stats <- function(data) {
  data |>
    dplyr::group_by(metabolite) |>
    dplyr::summarise(dplyr::across(value, list(mean = mean, sd = sd))) |>
    dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), ~ round(.x, digits = 1)))
}


#' Simple geom_histogram plot
#'
#' @param data Lipidomics data set
#'
#' @return ggplot object
plot_distributions <- function(data){
  data |> ggplot2::ggplot(ggplot2::aes(x = value)) +
    ggplot2::geom_histogram() +
    ggplot2::facet_wrap(vars(metabolite), scales = "free")
}
