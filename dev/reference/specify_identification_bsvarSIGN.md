# R6 Class Representing IdentificationBSVARSIGN

The class IdentificationBSVARSIGN presents the identifying restrictions
for the Bayesian Structural VAR models with sign and narrative
restrictions.

## Public fields

- `VB`:

  a list of `N` matrices determining the unrestricted elements of matrix
  \\B\\.

- `sign_irf`:

  a `NxNxH` array of sign restrictions on the impulse response
  functions.

- `sign_narrative`:

  a `ANYx6` matrix of narrative sign restrictions.

- `sign_structural`:

  a `NxN` matrix of sign restrictions on contemporaneous relations.

- `max_tries`:

  a positive integer with the maximum number of iterations for finding a
  rotation matrix \\Q\\ that would satisfy sign restrictions.

## Methods

### Public methods

- [`specify_identification_bsvarSIGN$new()`](#method-IdentificationBSVARSIGN-new)

- [`specify_identification_bsvarSIGN$get_identification()`](#method-IdentificationBSVARSIGN-get_identification)

- [`specify_identification_bsvarSIGN$set_identification()`](#method-IdentificationBSVARSIGN-set_identification)

- [`specify_identification_bsvarSIGN$clone()`](#method-IdentificationBSVARSIGN-clone)

------------------------------------------------------------------------

### Method `new()`

Create new identifying restrictions IdentificationBSVARSIGN.

#### Usage

    specify_identification_bsvarSIGN$new(
      N,
      sign_irf,
      sign_narrative,
      sign_structural,
      max_tries = Inf
    )

#### Arguments

- `N`:

  a positive integer - the number of dependent variables in the model.

- `sign_irf`:

  a `NxNxH` array - sign and zero restrictions on the impulse response
  functions, ±1 for positive/negative sign restriction 0 for zero
  restrictions and NA for no restrictions, the `h`-th slice `NxN` matrix
  contains the restrictions on the `h-1` horizon.

- `sign_narrative`:

  a list of objects of class "narrative" - narrative sign restrictions.

- `sign_structural`:

  a `NxN` matrix with entries ±1 or NA - sign restrictions on the
  contemporaneous relations `B` between reduced-form errors `E` and
  structural shocks `U` where `BE=U`.

- `max_tries`:

  a positive integer with the maximum number of iterations for finding a
  rotation matrix \\Q\\ that would satisfy sign restrictions.

#### Returns

Identifying restrictions IdentificationBSVARSIGN.

------------------------------------------------------------------------

### Method `get_identification()`

Returns the elements of the identification pattern
IdentificationBSVARSIGN as a `list`.

#### Usage

    specify_identification_bsvarSIGN$get_identification()

------------------------------------------------------------------------

### Method `set_identification()`

Set new starting values StartingValuesBSVARSIGN.

#### Usage

    specify_identification_bsvarSIGN$set_identification(
      N,
      sign_irf,
      sign_narrative,
      sign_structural
    )

#### Arguments

- `N`:

  a positive integer - the number of dependent variables in the model.

- `sign_irf`:

  a `NxNxH` array - sign and zero restrictions on the impulse response
  functions, ±1 for positive/negative sign restriction 0 for zero
  restrictions and NA for no restrictions, the `h`-th slice `NxN` matrix
  contains the restrictions on the `h-1` horizon.

- `sign_narrative`:

  a list of objects of class "narrative" - narrative sign restrictions.

- `sign_structural`:

  a `NxN` matrix with entries ±1 or NA - sign restrictions on the
  contemporaneous relations `B` between reduced-form errors `E` and
  structural shocks `U` where `BE=U`.

- `max_tries`:

  a positive integer with the maximum number of iterations for finding a
  rotation matrix \\Q\\ that would satisfy sign restrictions

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    specify_identification_bsvarSIGN$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
# recursive specification for a 5-variable system
specify_identification_bsvarSIGN$new(N = 5)
#> <IdentificationBSVARSIGN>
#>   Public:
#>     VB: list
#>     clone: function (deep = FALSE) 
#>     get_identification: function () 
#>     initialize: function (N, sign_irf, sign_narrative, sign_structural, max_tries = Inf) 
#>     max_tries: Inf
#>     set_identification: function (N, sign_irf, sign_narrative, sign_structural) 
#>     sign_irf: NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA ...
#>     sign_narrative: list
#>     sign_structural: 1 NA NA NA NA NA 1 NA NA NA NA NA 1 NA NA NA NA NA 1 NA  ...

# specify sign restrictions of the first shock on the contemporaneous IRF
# + no effect on the first variable
# + positive effect on the second variable
sign_irf = matrix(c(0, 1, rep(NA, 23)), 5, 5)
specify_identification_bsvarSIGN$new(N = 5, sign_irf = sign_irf) 
#> <IdentificationBSVARSIGN>
#>   Public:
#>     VB: list
#>     clone: function (deep = FALSE) 
#>     get_identification: function () 
#>     initialize: function (N, sign_irf, sign_narrative, sign_structural, max_tries = Inf) 
#>     max_tries: Inf
#>     set_identification: function (N, sign_irf, sign_narrative, sign_structural) 
#>     sign_irf: 0 1 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA N ...
#>     sign_narrative: list
#>     sign_structural: NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA ...
```
