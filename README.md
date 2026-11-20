# NHL Statistical Modeling and Predictive Analytics


An academic study of NHL game outcomes during the 2021–2022 season, developed in the context of an undergraduate internship. The workflow covers data collection, exploratory analysis, statistical modeling and historical out-of-sample evaluation.

## Methods
- HTML data extraction, cleaning and match-level transformation.
- Double Poisson regression with team attack/defense and home advantage.
- Zero-truncated Skellam likelihood, numerical maximum likelihood estimation and uncertainty analysis.
- Regular-season/playoff and chronological holdout evaluation using Brier score, accuracy and log loss.
- Retrospective expected-value analysis against closing market odds; not a recommendation to wager.

## Data
The data were obtained for the academic analysis from a publicly accessible website, not from an internal Kaizen Gaming dataset. The scripts use the [Sportsbook Reviews Online NHL 2021–22 archive](https://www.sportsbookreviewsonline.com/scoresoddsarchives/nhl-odds-2021-22/). The raw data are not redistributed here. Script 01 downloads the HTML if it is not already present. If the source is unavailable, place an authorized copy named `sbr_nhl_2021-22.html` in the repository root. The source parser expects two rows per match. Check the source terms before retrieving or redistributing data. 

## Run
Install R and run these commands from the repository root:

```sh
Rscript install_packages.R
Rscript run_all.R
```

The five numbered scripts can also be run in order. CSVs and PNG figures are generated in the repository root. The workflow uses `rvest`, `moments`, `multcomp` and `ggplot2`.

## Limitations
A single historical season does not establish generalization or profitability. Closing-market benchmarks contain information not available to the score-only models. Independent score and zero-truncation assumptions simplify overtime outcomes. The workflow validates row pairing and odds and stops if output files cannot be written. Full execution in R with the source HTML has not been independently verified. AIC values from models fitted to different responses (paired scores versus goal differences) are not directly comparable; use common held-out probability metrics. The t-tests and bootstrap are exploratory and assume independent outcomes; repeated teams and threshold exploration limit inference.

## Academic context
Undergraduate statistics project. Publish only with the appropriate university and internship permissions. No internship report or internal company data is included. Public accessibility of the source is distinct from permission to redistribute data or publish internship-related work.
