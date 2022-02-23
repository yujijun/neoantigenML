# INSTALL PACKAGES & LOAD LIBRARIES -----------------
cat("INSTALLING PACKAGES & LOADING LIBRARIES... \n\n", sep = "")
packages <- c("Biostrings",
              "car",
              "caret",
              "cvAUC",
              "DescTools",
              "doParallel",
              "doSNOW",
              "foreach",
              "fst",
              "ggpubr",
              "ggsci",
              "gridExtra",
              "igraph",
              "matrixStats",
              "mlr",
              "msa",
              "precrec",
              "pROC",
              "psych",
              "randomForestSRC",
              "rlecuyer",
              "S4Vectors",
              "seqinr",
              "snow",
              "stringdist",
              "survminer",
              "zoo",
              "pbapply",
              "Peptides",
              "VennDiagram") # list of packages to load
n_packages <- length(packages) # count how many packages are required
new.pkg <- packages[!(packages %in% installed.packages())] # determine which packages aren't installed
# install missing packages
if(length(new.pkg)){
  install.packages(new.pkg)
}
cat("ALL PACKAGES HAVE BEEN INSTALLED... \n\n", sep = "")
#--------------Main -----------------------
