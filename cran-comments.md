## Submission notes bsvarSIGNs v1.0

## R CMD check results

There were no ERRORs, WARNINGs, or NOTEs.

## GitHub R-CMD-check using `usethis::use_github_action_check_standard()`

Passing on all platforms!

## Check at using `devtools::check(manual = TRUE, remote = TRUE, incoming = TRUE)`

This shows:
>  Maintainer: ‘Xiaolei Wang <adamwang15@gmail.com>’
>  Found the following (possibly) invalid URLs:
>    URL: https://www.linkedin.com/in/tomaszwwozniak
>      From: README.md
>      Status: 999
>    URL: https://www.linkedin.com/in/xiaolei-adam-wang/
>      From: README.md
>      Status: 999

This is not a problem with the link, but how LinkedIn responds to automatic checks as documented e.g. [HERE](https://stackoverflow.com/questions/27231113/999-error-code-on-head-request-to-linkedin) and [HERE](https://http.dev/999)

It also shows:
>  Package in Depends/Imports which should probably only be in LinkingTo: ‘RcppArmadillo’

This is not a problem as `RcppArmadillo` is used in the package and is not a problem to be in `Depends` or `Imports`. We tested various scenarios of includeing `RcppArmadillo` in dependencies, and only the current setup compiles and installs the package correctly. It shows this note nevertheless.


## Done some more tests from `usethis::use_release_issue()`

All good here!