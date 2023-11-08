install.packages("pak")

pak::pak("rostools/r3")

r3::install_packages_advanced()

# There will be a pop-up to type in your name (first and 
# last), as well as your email
r3::setup_git_config()

r3::check_setup()


prodigenr::setup_with_git()
usethis::use_blank_slate("project")
usethis::use_r("functions", open = FALSE)

r3::create_qmd_doc()

usethis::use_data_raw("nmr-omics")

lipidomics <- readr::read_csv(here::here("data/lipidomics.csv"))

lipidomics

usethis::use_git_ignore("data-raw/nmr-omics/")

r3::check_project_setup_advanced()
