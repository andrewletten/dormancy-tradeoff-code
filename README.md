
<!-- README.md is generated from README.Rmd. Please edit that file -->

## Microbial dormancy supports multi-species coexistence under resource fluctuations

Code and simulation output from:

Letten, Yamamichi, Richardson and Ke. Microbial dormancy supports
multi-species coexistence under resource fluctuations.

#### Core scripts

`dormancy-tradeoff-fig1acd.R` - Simulate competition between a gleaner
and a dormancy strategist (Fig 1a,c&d in ms).

`dormancy-tradeoff-fig1b.R` - Simulate competition between a gleaner and
a dormancy strategist with different parameter combinations for the
dormancy strategist’s growth and dormant phase mortality (Fig 1b in ms).

`dormancy-tradeoff-fig2.R` - Simulate competition between a gleaner, an
opportunist and a dormancy strategist for one resource supply
concentration and pulsing frequency (Fig 2 in ms).

`dormancy-tradeoff-3sppsims.R` - Simulate competition between all
combinations of a gleaner, opportunist and dormancy strategist under
different patterns of resource supply.

`dormancy-tradeoff-3sppsims-chemo.R` - Simulate competition between all
combinations of a gleaner, opportunist and dormancy strategist under
continuous/chemostat resource supply.

`dormancy-tradeoff-fig3.R` - Summarize and plot simulation output from
`dormancy-tradeoff-3sppsims.R` and `dormancy-tradeoff-3sppsims-chemo.R`
(Fig 3 in ms).

#### Other scripts

`dormancy-tradeoff-spp-params.R` - Parameters for core analyses.

`dormancy-tradeoff-functions.R` - Functions for defining ODE, extracting
and wrangling saved sims, plotting etc.

#### Simulation outputs

The runtime for the core simulations in `dormancy-tradeoff-3sppsims.R`
(Fig 3 in ms) is likely to be on the order of days (or at least many
hours). To save having to rerun all simulations, all outputs from
`dormancy-tradeoff-3sppsims.R` and `dormancy-tradeoff-3sppsims-chemo.R`
are provided as individual zipped csv files in the latest
[release](https://github.com/andrewletten/dormancy-tradeoff-code/releases/tag/latest)
associated with this repo. As described below, these can be downloaded
using the `piggyback` package, and distributed into the appropriate
directory structure for importing into `dormancy-tradeoff-fig3.R`. Note,
the download time for each of the four batches of sim files (all 3 spp
and the three pairwise dynamics) is likely to be on the order of hours
(or at least many mins, but not days!). Note also that the code in
`dormancy-tradeoff-fig3.R` can tested on one batch of data at a time
(i.e. there is no need to download all four batches at the same time).
Finally, please note there are 960 zipped csv files in the release with
a combined size of ~15GB.

``` r
library(piggyback)

simnames <- pb_list(repo = "andrewletten/dormancy-tradeoff-code", tag = "latest")
# Make a new directory in the data directory in the cloned repo
dir.create("data/N1N2N3_gz", recursive = TRUE) # all 3 spp
dir.create("data/N1N2_gz", recursive = TRUE) # gleaner and dormancy strategist
dir.create("data/N1N3_gz", recursive = TRUE) # opportunist and dormancy strategist
dir.create("data/N2N3_gz", recursive = TRUE) # gleaner and opportunist

N1N2N3_files <- simnames$file_name[grep("^N1N2N3_", simnames$file_name)]
pb_download(
  file = N1N2N3_files,
  dest = "data/N1N2N3_gz",
  repo = "andrewletten/dormancy-tradeoff-code",
  tag = "latest"
)

N1N2_files <- simnames$file_name[grep("^N1N2_", simnames$file_name)]
pb_download(
  file = N1N2_files,
  dest = "data/N1N2_gz",
  repo = "andrewletten/dormancy-tradeoff-code",
  tag = "latest"
)

N1N3_files <- simnames$file_name[grep("^N1N3_", simnames$file_name)]
pb_download(
  file = N1N3_files,
  dest = "data/N1N3_gz",
  repo = "andrewletten/dormancy-tradeoff-code",
  tag = "latest"
)

N2N3_files <- simnames$file_name[grep("^N2N3_", simnames$file_name)]
pb_download(
  file = N2N3_files,
  dest = "data/N2N3_gz",
  repo = "andrewletten/dormancy-tradeoff-code",
  tag = "latest"
)
```

Once the sim files are downloaded, run `dormancy-tradeoff-fig3.R` to
reproduce the different iterations (species combinations) in the panels
of Fig 3.
