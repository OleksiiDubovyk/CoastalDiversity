# Detection of diversity of coastal avian and mammalian fauna

*Oleksii Dubovyk, Ella DiPetto, Chi Wei, Alex Wright, Iroshmal Peiris, Maizer Sparkman, Eric L. Walters*

Data and analyses for the VAS 2024 annual meeting presentation (and subsequent pubs, hopefully) on observed avian and mammalian diversity of shorelines in coastal Virginia. Data privided by Ella DiPetto.

## Data source
The data collected during the dissertation research of [Ella DiPetto](https://edipetto.weebly.com/).

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

## Prerequisites
### R stuff
- The latest [R version](https://cran.r-project.org/bin/windows/base/)
- [Posit/RStudio](https://posit.co/download/rstudio-desktop/)
- Install `tidyverse`, `lubridate`:
```r
packages <- c("tidyverse", "lubridate")
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
    - Go back to Git Bash an run
    ```bash session
    gh auth login
    ```
    - Follow the prompts: GitHub.com -> HTTPS -> Y -> Login with a web browser
- Open Git Bash, run
  ```bash session
  git remote add origin https://{TOKEN}@github.com/OleksiiDubovyk/CoastalDiversity.git/
  ```
- Whenever you want to edit the code, type
  ```bash session
  git pull
  git add filename.extension # specify the file you have just changed
  git commit -m "Your comments on what you've added"
  git push
  ```
