## CRAN submission comments for SDALGCP version 0.5.0

This is a resubmission of the SDALGCP package, which was previously archived on CRAN on 2022-05-04 due to its dependency on the now-archived `geoR` package.

### Summary of changes since last CRAN version:

- The package has been fully revised to be CRAN-compliant.
- Support for `sf`, `terra`, and `ggplot2` has been added, replacing legacy spatial dependencies.
- Replaced deprecated or archived functions and improved performance and compatibility with current spatial packages.
- Documentation has been rewritten to meet current CRAN guidelines.
- Examples and internal checks have been revised for compatibility and performance.
- `geoR` is now listed as an optional dependency only for specific model support.
- Improved support for spatio-temporal modeling and visualisation.
- Functions have been refactored to ensure better modularity and maintainability.

### R CMD check results (on Ubuntu 22.04 with R 4.5.1):

- 0 errors ✔
- 0 warnings ✔
- 1 NOTE ✖

The NOTE is:

