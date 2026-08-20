# A 6-variable US monetary policy data, from 1965 Jan to 2007 Aug

A sample data to identify monetary policy shock.

## Usage

``` r
data(monetary)
```

## Format

A matrix and a `ts` object with time series of over two hundred
observations on 5 variables:

- gdpc1:

  monthly real gross domestic product

- gdpdef:

  monthly gross domestic product: implicit price deflator

- cprindex:

  monthly consumer price index

- totresns:

  monthly reserves of depository institutions

- bognonbr:

  monthly non-borrowed reserves of depository institutions

- fedfunds:

  monthly federal funds effective rate

## Source

Replication package,
<https://www.aeaweb.org/articles?id=10.1257/aer.20161852>

## References

Antolín-Díaz & Rubio-Ramírez (2018) Narrative Sign Restrictions for
SVARs, American Economic Review, 108(10), 2802-29,
\<doi:10.1257/aer.20161852\>.
