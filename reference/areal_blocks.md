# Contiguous areal blocks for the spatial-block bootstrap.

Partitions the CAR units into `nblk` contiguous blocks. Disconnected
components are kept as separate blocks; large components are further
split by a coarse spatial k-means on a 2-D MDS/grid embedding of the
adjacency (no igraph dependency). Returns a length-n factor (one block
label per observation) so whole contiguous clusters of units are
resampled together (Lahiri 2003: under-covers with few blocks;
disclosed).

## Usage

``` r
areal_blocks(obj, nblk = 8L)
```
