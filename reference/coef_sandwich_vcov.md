# Classification-conditional smoothed-QR sandwich for one regime's coefficients. Bread-meat-bread form with the density read off the fitted smoothed kernel (never an ALD stand-in). Disclosed as classification-conditional.

On the CAR spatial-error path `Xt`/`Pen` arrive as sparse Matrix
`dgCMatrix` objects; the bread and meat are dense `(p + L')^2` blocks,
so we densify the inputs (`as.matrix`) and work with base linear
algebra. This is acceptable because the sandwich already materialises a
dense `(p + L')^2` inverse; for very large `L'` (areal L up to ~400 is
fine) the bootstrap is the recommended inference path.

## Usage

``` r
coef_sandwich_vcov(Xt, y, tau, w, beta, Pen, h)
```
