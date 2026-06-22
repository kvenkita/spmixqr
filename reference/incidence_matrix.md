# Unit-incidence matrix R (n x L), sparse, from an observation-\>unit index.

Unit-incidence matrix R (n x L), sparse, from an observation-\>unit
index.

## Usage

``` r
incidence_matrix(unit_idx, L)
```

## Arguments

- unit_idx:

  Length-n integer index of each observation's unit (1..L).

- L:

  Number of units.

## Value

A sparse `n x L` incidence matrix.
