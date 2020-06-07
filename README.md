# covid-19
COVID 19 Case Projection Model for Turkey

## Usage

To fit data using all observations up to latest known day and project 20 days in the future: 
```sh
Rscript ModelFit.R
```

To fit using observations up to a certain date (e.g. May 12th, 2020):

```sh
Rscript ModelFit.R 2020-05-12
```

To fit using observations up to a certain date (e.g. May 12th, 2020), and project e.g. 10 days into the future from that date:

```sh
Rscript ModelFit.R 2020-05-12 10
```

The fit and relevant information will be saved to `./models/`, the plots to `./plots/` and tabulated projections to `./tables/` once run is complete.
