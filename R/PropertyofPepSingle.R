# This is a script for property of peptides.
# Update Time:Thu Nov 25 18:38:05 2021
# Author: JijunYu


#### property for whole peptides ######
#' Title Properties of whole peptides
#'
#' @param peptides A vector of peptides with same length
#' @param autoCor A logic value TRUE or FALSE, meaning to calculate autoCovariance value or not.
#' @param lag A value for a lag, the max value is equal to the length of shortest peptide minus one. detail ref 'Peptides' package.
#' @param property A property to use as value to be correlated.
#' @param property1 A property to use as value to evaluate the cross-covariance.
#' @param property2 A property to use as value to evaluate the cross-covariance.
#' @param center  A logic value TRUE or FALSE if the property must be centered
#' @param pH A pH value in charge function.
#' @param pKscale A character string specifying the pKa scale to be used; must be one of "Bjellqvist","Dawson", "EMBOSS", "Lehninger", "Murray", "Rodwell", "Sillero", "Solomon"or "Stryer"
#' @param scale A character string specifying the hydophobicity scale to be used; must be one of "Aboderin", "AbrahamLeo", "Argos", "BlackMould", "BullBreese", "Casari","Chothia", "Cid", "Cowan3.4", "Cowan7.5", "Eisenberg", "Engelman", "Fasman","Fauchere", "Goldsack", "Guy", "HoppWoods", "Janin", "Jones", "Juretic","Kidera", "Kuhn", "KyteDoolittle", "Levitt", "Manavalan", "Miyazawa","Parker", "Ponnuswamy", "Prabhakaran", "Rao", "Rose", "Roseman", "Sweet","Tanford", "Welling", "Wilson", "Wolfenden", "Zimmerman", "interfaceScale_pH8","interfaceScale_pH2", "octanolScale_pH8", "octanolScale_pH2", "oiScale_pH8"or "oiScale_pH2".
#' @param monoisotopic A logical value 'TRUE' or 'FALSE' indicating if monoisotopic weights of aminoacids should be used
#' @param avgScale Set the mass scale to use for average weight only (if ’monoisotopic == FALSE’).Accepts "expasy" (default) or "mascot".
#' @param aaShift Define the mass difference in Dalton of given amino acids as a named vector.Use the amino acid one letter code as names and the mass shift in Dalton as values.
#' @param charge The net charge for which the m/z should be calculated
#' @param pKscale.PI A character string specifying the pK scale in PI function to be used; must be one of "Bjellqvist","EMBOSS", "Murray", "Sillero", "Solomon", "Stryer", "Lehninger", "Dawson"or "Rodwell"
#'
#' @return a dataframe with all properties.
#' @export
#'
#' @examples
#' Neodataset.length <- Neodataset %>%
#' mutate(Length = str_length(wild_Peptide)) %>%
#' filter(Length == 9) %>%
#' unique()
#' peptides <- unique(Neodataset.length$wild_Peptide)
#' PropertyofPep(peptides)
PropertyofPepSingle <- function(peptides,
                          PropAll = c("aaComp","aIndex",
                                      "autoCorrelation","autoCovariance",
                                      "blosumIndices","boman",
                                      "charge","crossCovariance",
                                      "crucianiProperties","fasgaiVectors",
                                      "hmoment","hydrophobicity",
                                      "instaIndex","kideraFactors",
                                      "mswhimScores",
                                      "mw","mz","pl","protFP",
                                      "stScales","tScales","vhseScales",
                                      "zScales"),
                          autoCor = TRUE,
                          lag = 1,
                          property = AAdata$Hydrophobicity$KyteDoolittle,
                          property1 = AAdata$Hydrophobicity$KyteDoolittle,
                          property2 = AAdata$Hydrophobicity$Eisenberg,
                          center = TRUE,
                          pH = 7,
                          pKscale = "Lehninger",
                          scale = "KyteDoolittle",
                          monoisotopic = FALSE,
                          avgScale = "expasy",
                          aaShift = NULL,
                          charge = 2,
                          pKscale.PI = "EMBOSS"){
  library(Peptides)
  library(tidyverse)
  source("./R/base.r")

  data(AAdata)
  if("aaComp" %in%  PropAll) {
    aaComp.pep <- aaComp(seq = peptides)  # A List + matrix: compose of AA
  }
  if("aIndex" %in%  PropAll){
    aIndex.pep <- aIndex(seq = peptides) # A vector: thermal stability of proteins
  }
  if("autoCorrelation" %in% PropAll){
    autoCorrelation.pep <- autoCorrelation(sequence = peptides,
                                           lag = lag,
                                           property = property,
                                           center = center) # a vector
  }
  if("autoCovariance" %in% PropAll){
    autoCovariance.pep <- autoCovariance(sequence = peptides,
                                         lag = lag,
                                         property = property,
                                         center = center) # a vector
  }
  if("crossCovariance" %in% PropAll){
    crossCovariance.pep <- crossCovariance(sequence = peptides,
                                           lag = lag,
                                           property1 = property1,
                                           property2 = property2,
                                           center = center) #
  }
  if("blosumIndices" %in% PropAll){
    blosumIndices.pep <- blosumIndices(seq = peptides) # a list + vector
  }
  if("boman" %in% PropAll){
    boman.pep <- boman(seq = peptides) # a vector:Potential Protein Interaction)
  }
  if("charge" %in% PropAll){
    charge.pep <- charge(seq = peptides,
                         pH = pH,
                         pKscale = pKscale) #  A vector theoretical net charge
  }
  if("crucianiProperties" %in% PropAll){
    crucianiProperties.pep <- crucianiProperties(seq = peptides) # A list + vector
  }
  if("fasgaiVectors" %in% PropAll){
    fasgaiVectors.pep <- fasgaiVectors(seq = peptides) #A list + vector
  }
  if("hmoment" %in% PropAll){
    hmoment.pep <- hmoment(seq = peptides) # A vector Compute the hydrophobic moment of a protein sequence
  }
  if("hydrophobicity" %in% PropAll){
    hydrophobicity.pep <- hydrophobicity(seq = peptides,
                                         scale = scale) # A vector
  }
  if("instaIndex" %in% PropAll){
    instaIndex.pep <- instaIndex(seq =peptides) # A vector: Compute the instability index
  }
  if("kideraFactors" %in% PropAll){
    kideraFactors.pep <- kideraFactors(seq = peptides) #A list + vector
  }
  if("mswhimScores" %in% PropAll){
    mswhimScores.pep <- mswhimScores(seq = peptides) # A list + a vector 36 electrostatic potential properties
  }
  if("mw" %in% PropAll){
    mw.pep <- mw(seq = peptides,
                 monoisotopic = monoisotopic,
                 avgScale = avgScale,
                 label = "none",
                 aaShift = aaShift) # A vector
  }
  if("mz" %in% PropAll){
    mz.pep <- mz(seq = peptides,
                 charge = charge,
                 aaShift = aaShift,
                 cysteins = 57.021464) #A vector
  }
  if("pl" %in% PropAll){
    pI.pep <- pI(seq = peptides,
                 pKscale = pKscale.PI) #A vector
  }
  if("protFP" %in% PropAll){
    protFP.pep <- protFP(seq = peptides) # A list + a vector
  }
  if("stScales" %in% PropAll){
    stScales.pep <- stScales(seq = peptides) # A list + a vector
  }
  if("tScales" %in% PropAll){
    tScales.pep <- tScales(seq = peptides)# A list+ a vector
  }
  if("vhseScales" %in% PropAll){
    vhseScales.pep <- vhseScales(seq = peptides)# A list+ a vector
  }
  if("zScales" %in% PropAll){
    zScales.pep <- zScales(seq = peptides)# A list+ a vector
  }
  #### convert list into dataframe
  allvar <- grep(".pep",ls(),value = T)
  property.list <- list()
  for(var in allvar){
    if(is.list(get(var)) == TRUE){
      # convert list into dataframe
      assign(var, listtodf(get(var),name = peptides))
    }else if(is.vector(get(var)) == TRUE){
      # convert vector into dataframe
      vectortodf <- as.data.frame(get(var))
      colnames(vectortodf) <- var
      assign(var,vectortodf)
    }
    property.list[[var]] <- get(var)
  }

  all.df <- do.call(cbind.data.frame,property.list)
  rownames(all.df) <- peptides
  return(all.df)
}

