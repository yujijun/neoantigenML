# This is an prediction model immunogenicity based on the physical and chemical properties of peptides by machine learning method.

## 1. Prepare Environment
```
> sessioninfo()
R version 4.0.2 (2020-06-22)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS  10.16

attached base packages:
  [1] parallel  stats     graphics  grDevices utils     datasets  methods   base

other attached packages:
  [1] NeoantigenML_0.0.0 foreach_1.5.1      snow_0.4-4         data.table_1.14.2
[5] forcats_0.5.1      stringr_1.4.0      dplyr_1.0.7        purrr_0.3.4
[9] readr_2.1.1        tidyr_1.1.4        tibble_3.1.6       ggplot2_3.3.5
[13] tidyverse_1.3.1    corrplot_0.92      devtools_2.4.3     usethis_2.1.5

loaded via a namespace (and not attached):
[1] Rcpp_1.0.7        lubridate_1.8.0   prettyunits_1.1.1 ps_1.6.0          assertthat_0.2.1
[6] rprojroot_2.0.2   utf8_1.2.2        cellranger_1.1.0  R6_2.5.1          backports_1.4.1
[11] reprex_2.0.1      httr_1.4.2        pillar_1.6.4      rlang_0.4.12      readxl_1.3.1
[16] rstudioapi_0.13   callr_3.7.0       desc_1.4.0        munsell_0.5.0     broom_0.7.10
[21] compiler_4.0.2    modelr_0.1.8      pkgconfig_2.0.3   pkgbuild_1.3.1    doSNOW_1.0.19
[26] tidyselect_1.1.1  codetools_0.2-18  fansi_0.5.0       tzdb_0.2.0        crayon_1.4.2
[31] dbplyr_2.1.1      withr_2.4.3       grid_4.0.2        jsonlite_1.7.2    gtable_0.3.0
[36] lifecycle_1.0.1   DBI_1.1.2         magrittr_2.0.1    scales_1.1.1      stringi_1.7.6
[41] cli_3.1.0         cachem_1.0.6      pbapply_1.5-0     fs_1.5.2          remotes_2.4.2
[46] testthat_3.1.1    xml2_1.3.3        ellipsis_0.3.2    generics_0.1.1    vctrs_0.3.8
[51] iterators_1.0.13  tools_4.0.2       glue_1.6.0        hms_1.1.1         processx_3.5.2
[56] pkgload_1.2.4     fastmap_1.1.0     colorspace_2.0-2  sessioninfo_1.2.2 rvest_1.0.2
[61] memoise_2.0.1     haven_2.4.3
```
User could install all related by 00_InstallPackage.R under ./main_running_v3 folder


## 2. Running
Running propress could refer to files under `./main_running_v3/` folder

## 3. Introduction about folder ./R 
#This is an neoantigens filtering repository by machine learning method

** Description about function:**
1. Main function for whole project was in main.r function;
2. All preprocess about datasets were inclueded into dataset.r;
3. All functions were put into *.R;
4. All supplementary functions used in *.R were included into base.r;

