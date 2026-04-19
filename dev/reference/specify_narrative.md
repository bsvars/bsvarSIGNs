# vector specifying one narrative restriction

The class narrative specifies a single narrative restriction.

## Usage

``` r
specify_narrative(
  start,
  periods = 1,
  type = "S",
  sign = 1,
  shock = 1,
  var = NA
)
```

## Arguments

- start:

  positive integer - the period in which the narrative starts (greater
  than the number of lags).

- periods:

  positive integer - the number of periods the narrative restriction
  lasts.

- type:

  character - the type of the narrative restriction (one of "S", "A",
  "B"), where: "S" for restrictions on structural shocks; "A" for type A
  restrictions on historical decomposition, i.e. if the absolute value
  of the historical decomposition of shock to var is the greatest/least
  among all shocks; "B" for type B restrictions on historical
  decomposition, i.e. if the absolute value of the historical
  decomposition of shock to var is the greater/less than the sum of all
  other shocks.

- sign:

  integer - the sign of the narrative restriction (1 or -1).

- shock:

  positive integer - the index of the shock to which the narrative
  restriction applies.

- var:

  positive integer - the index of the variable to which the narrative
  restriction applies.

## Value

An object of class `narrative` specifying one narrative restriction.

## Examples

``` r
# specify a narrative restriction
narrative       = specify_narrative(
                    start = 166, 
                    periods = 1, 
                    type = "S", 
                    sign = 1, 
                    shock = 1, 
                    var = 6
                  )
# use it to specify the model
specification   = specify_bsvarSIGN$new(monetary, sign_narrative = list(narrative))
```
