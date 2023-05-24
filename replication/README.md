# Replication files

# DOTS data

## Covariate data

Covariate data is included in `data_X_1995.csv`
and similarly for the other years.
The columns are as follow:

- `ISO`
  Three-letter code for each country

- `Population`
  The population of each country

- `Population_bracket`
  The decile of each country's population

- `GDP`
  The GDP of each country

- `GDP_bracket`
  The decile of each country's GDP

- `GDP_per_capita`
  The GDP per capita of each country

- `GDP_per_capita_bracket`
  The decile of each country's GDP per capita

- `GATT`
  Whether each country is in the General Agreement on Tariffs and Trade

## Trade data

Trade data is included in `data_W_1995.csv`
and similarly for the other years.
The data is organized as follows.

- The first row is a list of all the ISO country codes.
- The subsequent rows indicate the trade volume between countries,
  forming an upper triangular matrix with missing entries given as -Inf.
