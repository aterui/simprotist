
# build package -----------------------------------------------------------

usethis::use_mit_license()
usethis::use_roxygen_md()
usethis::use_package_doc()
devtools::document()
devtools::load_all()
devtools::check(vignettes=FALSE)


# check syntax ------------------------------------------------------------

lintr::lint_package()


# build website -----------------------------------------------------------

# Run once to configure package to use pkgdown
usethis::use_pkgdown()

# Run to build the website
pkgdown::build_site()

usethis::use_pkgdown_github_pages()
