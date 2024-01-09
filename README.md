# Water fetching and climate

This repository contains code and data used for in "Temperature and precipitation affect the water fetching time burden in Sub-Saharan Africa" by Abigail Harvey Paulos, David Carroll II, Julie Powers, Jake Campolo, Daehyun Daniel Kim, Avery Cohn, and Amy Pickering.

Preprint: https://www.researchsquare.com/article/rs-3789072/v1

This code was developed on R version 4.3.1

Required R packages:
* data.table
* sandwich
* lm.test
* tidyverse
* doParallel (for parallelization, when appropriate)

Included are:
1. The water fetching - climate combined dataset. Data are split into three separate files to get the size of each file to fall below 100 MBs. To combine, read in each file and use rbind() to combine
2. The codebook for the combined dataset, including all variable definitions
3. Code to run the fixed effects models on the full dataset
4. Code to run the spatial first differences models on the full dataset.
* Note - to run the full dataset through the spatial first differences code will take >12 hours