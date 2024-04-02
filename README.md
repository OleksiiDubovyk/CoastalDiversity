# CoastalDiversity

*Oleksii Dubovyk, Ella DiPetto, Chi Wei, Alex Wright, Iroshmal Peiris, Maizer Sparkman, Eric L. Walters*

Data and analyses for the VAS 2024 annual meeting presentation (and subsequent pubs, hopefully) on observed avian and mammalian diversity of shorelines in coastal Virginia. Data privided by Ella DiPetto.

## Data source
The data collected during the dissertation research of [Ella DiPetto](https://edipetto.weebly.com/).

## Prerequisites
### R stuff
- The latest [R version](https://cran.r-project.org/bin/windows/base/)
- [Posit/RStudio](https://posit.co/download/rstudio-desktop/)
- Install `tidyverse`:
```
install.packages("tidyverse")
```
### Git Bash (if you only want to get the newest code)
- [Install Git Bash](https://carpentries.github.io/workshop-template/)
- Open the folder you want our repository to be copied to, e.g.,
  ```
  cd /c/Users/username/CoastalDiversity
  ```
- Type
  ```
  git clone https://github.com/OleksiiDubovyk/CoastalDiversity
  ```
- Whenever you want to get the newest code, type
  ```
  git pull
  ```
### Setup your GitHub account (if you plan to contribute to coding)
- Create a GitHub account
- Install (GitHub CLI)[https://github.com/cli/cli?tab=readme-ov-file#installation]
    - Open Windows PowerShell and run
    ```
    winget install --id GitHub.cli
    ```
    - Go back to Git Bash an run
    ```
    gh auth login
    ```
    - Follow the prompts: GitHub.com -> HTTPS -> Y -> Login with a web browser
- Open Git Bash, run
  ```
  git remote add origin https://{TOKEN}@github.com/OleksiiDubovyk/CoastalDiversity.git/
  ```
  [comment]: <Token is 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>
- Whenever you want to edit the code, type
  ```
  git commit -m "Your comments on what you've added"
  git push
  ```
