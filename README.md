# Detection of diversity of coastal avian and mammalian fauna

*Oleksii Dubovyk, Ella DiPetto, Chi Wei, Iroshmal Peiris, Eric L. Walters*

Data and analyses for the VAS 2024 annual meeting presentation (and subsequent pubs, hopefully) on observed avian and mammalian diversity of shorelines in coastal Virginia. Data privided by Ella DiPetto.

## Data source

### Wildlife observations
The data collected during the dissertation research of [Ella DiPetto](https://edipetto.weebly.com/).

#### [detections.csv](detections.csv)

The main dataset containing the data on wildlife observations.

#### [deployments.csv](deployments.csv)

Supplementary dataset on duration of each deployment of trail cameras.

### Functional traits

#### [BirdFuncDat.txt](BirdFuncDat.txt) and [MamFuncDat.txt](MamFuncDat.txt)
The dataset used was EltonTraits 1.0[^eltontraits].

### Tides

#### [tides.csv](tides.csv)

The data on tides from the [NOAA](https://tidesandcurrents.noaa.gov/waterlevels.html?id=8638610&units=standard&bdate=20220401&edate=20230201&timezone=GMT&datum=MLLW&interval=h&action=data).

## Contents

### [main.R](main.R)

The main script, all analyses are here.

### [suntime.R](suntime.R)

A function to find a local time of sunrise and sunset for a given date and coordinates. Based on the [USGS calculator](https://gml.noaa.gov/grad/solcalc/calcdetails.html).

**Arguments**

  - `date` - chr, date of interest in "yyyy-mm-dd" format
  - `lat` - num, decimal latitude
  - `lon` - num, decimal longitude
  - `utc_offset` - num, time zone offset relative to the UTC: e.g., EDT is `-4`, EST is `-5`, PST is `-8`.

**Usage**

To calculate the sunrise and sunset time for Norfolk on Apr 2nd 2024, we call

```r
suntime(date = "2024-04-02", lat = 36.8794, lon = -76.2892, utc_offset = -4)
## [1] "06:48:22" "19:28:38"
```

### [bigrarefaction.R](bigrarefaction.R)

A group of functions to build rarefaction curves through interpolation or extrapolation procedures. `interpolation(..., mode = "l")` is built to handle large numbers.

It is assumed that for a community with $S$ species and $N$ individuals, such that there are $N_i$ individuals of species $i$ and, therefore, $\sum \limits_{i=1}^{S} N_i = N$, when $n$ individuals are drawn, the interpolated species richness ($S(x)$ represents the species richness observed when $x$ individuals are drawn) can be estimated as:

$S(n) = S(N) - {\binom{N}{n}}^{-1} \times \sum \limits_{i=1}^{S(N)} \binom{N-N_i}{n}$

and extrapolated values are estimated through Chao1 estimator[^chao], $\hat{f_0} = f_1^2 / 2f_2$, where $f_x$ represents the number of species for which $N_i = x$,

$S(N+m) = S(N) + \hat{f_0}\left[ 1 - \left( 1 - \frac{f_1}{N \hat{f_0} + f_1} \right)^m \right]$.

### [probrar.R](probrar.R)

*All questions regarding this section should be addressed to Oleksii, oadubovyk@gmail.com*

Desperate attempts to go away from the singleton/doubleton-based Chao[^chao] approximations of extrapolated rarefaction curves. Use at your own risk: the approach has not been peer reviewed and mostly relies on thoughts and prayers.

The notation is the following: `N` denotes a vector of values representing abundances of different species within a community. 

`prob_same(N, m)` estimates the expected probability of getting an unobserved before species in a sample drawn from community `N` at `m`th individual with replacement (therefore, inaccurate). If we let $S$ be the observed species richness, $J$ -- overall number of individuals, and $p_i = \frac{N_i}{J}$ -- percentages of species within a community, then the probability of the $m$th individual drawn to represent a new species is

$P(m) = \sum \limits_{i = 1}^{S} \prod \limits_{k = 1}^m \frac{N_i - p_i (k-1)}{J - (k-1)}$

Again, this estimation is inaccurate.

`probs_roll(N)` estimates the expected probabilities when drawing a sequence of individuals from a community simply applying `prob_same(N, m)` to the sequence of individuals $m = \{1, 2, 3, \dotsb, N-2, N-1, N\}$.

## Prerequisites
### R stuff
- The latest [R version](https://cran.r-project.org/bin/windows/base/)
- [Posit/RStudio](https://posit.co/download/rstudio-desktop/)
- Install `tidyverse`, `lubridate`, `data.table`, `caret`:
```r
packages <- c("tidyverse", "lubridate", "data.table", "caret")
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
invisible(lapply(packages, library, character.only = TRUE))
```
### Git Bash (if you only want to get the newest code)
- [Install Git Bash](https://carpentries.github.io/workshop-template/)
- Open the folder you want our repository to be copied to, e.g.,
  ```bash session
  cd /c/Users/username/CoastalDiversity
  ```
- Let the Git know who you are
  ```bash session
  git config --global user.name "Your Name"
  git config --global user.email "youremail@domain.com"
  ```
- Type
  ```bash session
  git clone https://github.com/OleksiiDubovyk/CoastalDiversity
  ```
- Whenever you want to get the newest code, type
  ```bash session
  git pull
  ```
### Setup your GitHub account (if you plan to contribute to coding)
- Create a GitHub account
- Install [GitHub CLI](https://github.com/cli/cli?tab=readme-ov-file#installation)
    - Open Windows PowerShell and run
    ```console
    winget install --id GitHub.cli
    ```
    - Restart Git Bash, navigate to the working directory, an run
    ```bash session
    gh auth login
    ```
    - Follow the prompts: GitHub.com -> HTTPS -> Y -> Login with a web browser
- Open Git Bash, run
  ```bash session
  git remote set-url origin https://{TOKEN}@github.com/OleksiiDubovyk/CoastalDiversity.git/
  ```
- Whenever you want to edit the code, type
  ```bash session
  git pull
  git add filename.extension # specify the file you have just changed
  git commit -m "Your comments on what you've added"
  git push
  ```

# References

[^eltontraits]: Wilman, H., J. Belmaker, J. Simpson, C. de la Rosa, M. M. Rivadeneira, and W. Jetz. 2014. EltonTraits 1.0: species-level foraging attributes of the world’s birds and mammals. Ecology 95:2027–2027. https://doi.org/10.1890/13-1917.1

[^chao]: Chao, A., N. J. Gotelli, T. C. Hsieh, E. L. Sander, K. H. Ma, R. K. Colwell, and A. M. Ellison. 2014. Rarefaction and extrapolation with Hill numbers: a framework for sampling and estimation in species diversity studies. Ecological Monographs 84:45–67. https://doi.org/10.1890/13-0133.1
