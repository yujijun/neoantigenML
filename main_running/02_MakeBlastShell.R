# BASIC INFO ---------------------------
##
## Script name: MakeBlastShell.R
##
## Purpose of script: create blastp script or info
##
## Author: JijunYu
##
## Date Created: 2021-12-08
## Update Date:
##
## Copyright (c) Jijun Yu, 2021
## Email: jijunyuedu@outlook.com
##
## Notes:
##
#--------------Main -----------------------
#### retrieve reffasta from Uniprot seq ####
load("./result/yjj/Tcell_deredundancy_all.distinct.RData")
Uniprot.name <- unique(Tcell_deredundancy_all.distinct$uniprot_5)
# library(Biostrings)
# uniprot.fasta <- readAAStringSet("./data-raw/uniprot_sprot.fasta", format="fasta",
#                                  nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
# for(i in Uniprot.name){
#   i_uniportname <- grep(x = names(uniprot.fasta),pattern = i,value = T)
#   if(is_empty(i_uniportname) == FALSE){
#     i_uniportfasta <- uniprot.fasta[i_uniportname]
#     output_name = paste0(i,".fasta")
#     writeXStringSet(i_uniportfasta,filepath = paste0("./result/iedb_fasta/",output_name))
#   }
# }
# #lots of fasta info couldn't be find from uniprotKB,download them by link in IEDB directly may be a better choose.

allfastawww <- paste0("wget https://www.uniprot.org/uniprot/",
                      Uniprot.name,".fasta")
allfastawww <- as.data.frame(allfastawww)
write_csv(allfastawww,file = "result/iedb_fasta_v1/allfastawww.sh",quote = NULL,col_names = F)

#### create makeblastdb shell code ####
fastafiles <- list.files("./result/iedb_fasta/")
fastafiles.name <- str_split_fixed(string = fastafiles,
                                   pattern = "[.]",
                                   n = Inf)[,1]
makedb <- paste0("makeblastdb -in ",fastafiles," -out ",fastafiles.name, " -dbtype prot")
makedb.df <- as.data.frame(makedb)
write_csv(makedb.df,file = "result/iedb_fasta_db/makedb.df.sh",
          quote = NULL,col_names = F)

#### make query fasta info ####
library(seqinr)
for(i in Uniprot.name){
  df.tmp <- Tcell_deredundancy_all.distinct %>%
    filter(uniprot_5 == i)
  peptides.tmp <- as.list(df.tmp$description)
  names(peptides.tmp) <- peptides.tmp
  fastaname = paste0("query.",i,".fasta")
  write.fasta(sequences = peptides.tmp,names = peptides.tmp,file.out = paste0("./result/iedb_query_fasta/",                                                                           fastaname))
}

#### make blastp shell code ####
blastp <- paste0("blastp -query ", "/mnt/sdc/neoantigenML/result/iedb_query_fasta/",
                 "query.",Uniprot.name,".fasta"," -db ", Uniprot.name,
                 " -out ","/mnt/sdc/neoantigenML/result/iedb_blast/",Uniprot.name, ".result",
                 " -outfmt 6 -num_threads 4")
blastp.df <- as.data.frame(blastp)
write_csv(blastp.df,file = "result/iedb_blast/iedb_blast.sh",
          quote = NULL,col_names = F)

