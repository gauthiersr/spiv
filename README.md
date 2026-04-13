# spiv: System Projections with Instrumental Variables

The `spiv` package provides an R implementation of the System Projections with Instrumental Variables (SP-IV) methodology using Local Projections.

## Methodology

This package implements the methods developed in: \* Lewis, D. J., & Mertens, K. (2024). *Dynamic Identification Using System Projections and Instrumental Variables*, available at <https://karelmertens.com/research/>. \* Methodology sourced from the Federal Reserve Bank of Dallas Working Paper 2204 : <https://doi.org/10.24149/wp2204>.

**Disclaimer:** This is an unofficial, community-contributed implementation of the SP-IV methodology. It is not affiliated with, nor endorsed by, the original authors or the Federal Reserve Bank of Dallas. The maintainer is solely responsible for any errors in the software implementation.

## Installation

You can install the development version of `spiv` from GitHub with:

```         
devtools::install_github("https://github.com/gauthiersr/spiv")
```

## Computational considerations

The grid search evaluates $L^K$ points (by brute forcing every combination), where $L$ = `grid_length`. When working with $K>2$,use a small grid (set `grid_length` to 20 or 30). AR should be prioritized over KLM if efficiency is a concern. Set `cores` to match your hardware capacity to distribute the load.

## Note on `plot`

Impulse Response Functions for both the outcome and endogenous variables are available in the output of the `spiv` function. The `plot` command only plots IRFs with HAC confidence intervals.
