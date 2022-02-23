
base_path <- "/mnt/sdc/neoantigenML"
setwd(base_path)
input.path = "./result/iedb_blast/"
output.path = "./result/yjj/"

allresult <- list.files(input.path)
allresult <- grep(x = allresult,pattern = "result",value = T)
redundancy <- function(x){
  tryCatch({
    input.tmp <- read.delim(file = paste0(input.path,x),header = F)
  }, error = function(e){
    NULL
  })
  if(!is.null(input.tmp)){
    # 消除内部匹配误差，留下匹配最好的部分
    input.tmp <- input.tmp %>%
      group_by(V1) %>%
      filter(V3 == max(V3)) %>%
      filter(V3 > 75) %>%
      filter(V4 == max(V4)) %>%
      ungroup()

    # 判断交叠的肽段，将交叠信息和原理的信息重合到一起
    if(nrow(input.tmp) != 0){
      library(intervals)
      from = as.matrix(input.tmp[,c(9,10)])
      from = Intervals(from,closed = c(TRUE, TRUE),type = "R")
      from.inter <- interval_overlap(from, from)
      from.inter <- t(do.call(cbind, lapply(lapply(from.inter, unlist), `length<-`, max(lengths(from.inter)))))
      from.inter[is.na(from.inter)] <- "0"
      from.inter  <- as.data.frame(from.inter) %>%
        tidyr::unite("overlap",everything(),sep = ",") %>%
        pull(overlap)
      from.inter <- str_replace_all(string = from.inter,pattern = ",0",replacement = "")
      input.tmp <- input.tmp  %>%
        mutate(overlap = from.inter) %>%
        rownames_to_column()

      # add iedb fasta info
      uniprot_tmp <- str_split(x,pattern = "[.]")[[1]][1]
      input.tmp <- input.tmp %>%
        mutate(uniprotid = uniprot_tmp)
    }else{
      input.tmp = NULL
    }
  }
  return(input.tmp)
}
#redundancy(x = "A0A0A0MRD9.result")
#for(i in list.files(input.path)){
#   print(i)
#   redundancy(i)
# }
resultfinal <- lapply(allresult,redundancy)
resultfinal.df <- do.call(rbind.data.frame,resultfinal)
resultfinal.df.bak <- resultfinal.df
resultfinal.df <- resultfinal.df %>%
  group_by(V1) %>%
  filter(V4 == max(V4)) %>%
  distinct(V1,V4,.keep_all =T) %>%
  ungroup()

write_tsv(resultfinal.df, file = "./result/yjj/deredundancy.blast.tsv",col_names = T,quote = NULL)


#### add the number of IEDB info ####
library(tidyverse)
library(Biostrings)
library(motifStack)
library(ggpubr)
input_folder <- "/Users/yujijun/Documents/work/4_AntigenML/NeoantigenML/input/"
output_folder <- "/Users/yujijun/Documents/work/4_AntigenML/NeoantigenML/result/yjj/"
Tcell_v3 <- read.delim(stringr::str_c(input_folder,"tcell_full_v3.csv"),sep = ",",header = T)
Tcell_v3_clean <- Tcell_v3 %>%
  select(matches("Epitope")|matches("Host")|matches("Assay")|matches("MHC"))
colnames(Tcell_v3_clean) <- Tcell_v3_clean[1,]
Tcell_v3_clean <- Tcell_v3_clean[-1,]
Tcell_v3_clean <- janitor::clean_names(Tcell_v3_clean)

redundancy <- read.delim(file = "./result/yjj/deredundancy.blast.tsv")

#### 添加肽段频率信息 ####
pep.num <- c()
for(i in redundancy$V1){
  i_num <- Tcell_v3_clean %>%
    filter(description == i) %>%
    nrow()
  pep.num <- c(pep.num,i_num)
} ### need some time
redundancy <- redundancy %>%
  mutate(pep.num = pep.num) %>%
  group_by(uniprotid) %>%
  add_count() %>%
  ungroup()
redundancy.single <- redundancy %>%
  filter(n == 1)
redundancy.du <- redundancy %>%
  filter(n != 1)

#### 判断同一个蛋白质肽段中，那些有重叠 ####
is.integer0 <- function(x)
{
  is.integer(x) && length(x) == 0L
}
tmp <- function(overlap1){
  if(is.integer0(which(overlap1 == 1))){
    return(overlap1)
  }else{
    index = seq(2,length(overlap1),1)
    i_all <- c()
    for(i in index){
      if(is.na(overlap1[i] ==1 & overlap1[i-1] != 1) |
         (overlap1[i] ==1 & overlap1[i-1] != 1)){
        i_all <- c(i_all,i,i-1)
      }
    }
    overlap1[i_all] = 1
    return(overlap1)
  }
}

library(groupdata2)
library(plotfunctions)
redundancy.du.test <- redundancy.du %>%
  group_by(uniprotid) %>%
  arrange(V9) %>%
  mutate(V9.1 = move_n_point(V10)) %>%
  mutate(V10_V9 = V9.1 - V9) %>%
  mutate(overlap1 = if_else(V10_V9 >= 3,1,0)) %>%
  mutate(overlap2 = tmp(overlap1)) %>% #将所有有交集的peptides 给找出来了
  ungroup() %>%
  arrange(uniprotid,V9)

#### 对有重叠的部分进行切分 ####
redundancy.du.test.overlap <- redundancy.du.test %>%
  mutate(overlap1 = if_else(is.na(overlap1), 0, overlap1)) %>%
  filter(overlap2 == 1) %>%
  group_split(uniprotid)

redundancy.du.test.nonoverlap <- redundancy.du.test %>%
  mutate(overlap2 = if_else(is.na(overlap2), 0, overlap2)) %>%
  filter(overlap2 != 1)

addgroup <- function(x){
  library(groupdata2)
  group = as.data.frame(group(1:nrow(x),n=which(x$overlap1 == 0),method = "l_starts"))
  group = group %>%
    pull(.groups)
  x <- x %>%
    mutate(group = group)
  return(x)
}
# for(i in 1:2){
#   print(i)
#   x = redundancy.du.test.overlap[[i]] %>%
#     addgroup()
#   print(x[,19:22])
# }
redundancy.du.test.overlap.list <- lapply(redundancy.du.test.overlap,addgroup)
redundancy.du.test.overlap.all <- do.call(rbind.data.frame,redundancy.du.test.overlap.list)

#### 将所有的结果整合到一起 ####
redundancy.du.test.overlap.all <- redundancy.du.test.overlap.all %>%
  group_by(uniprotid,group) %>%
  filter(pep.num == max(pep.num)) %>%
  top_n(n=1,wt = rowname) %>% # 对于肽段数量同样多的情况的处理
  ungroup()

all <- redundancy.du.test.overlap.all %>%
  full_join(redundancy.du.test.nonoverlap) %>% # 将去重后的和同一肽段里面没有重复的情况整合
  full_join(redundancy.single) %>%
  arrange(uniprotid)#和每个蛋白只有一个肽段的情况整合


#### 结果预处理，整合成标准格式 ####
load("./result/yjj/Tcell_deredundancy_all.distinct.RData")
Tcell_deredundancy_all.distinct.blast <- Tcell_deredundancy_all.distinct %>%
  filter(description %in% all$V1)
Tcell_deredundancy_all.distinct.blast.hostMus <- Tcell_deredundancy_all.distinct.blast %>%
  filter(!is.na(str_match(pattern = "Mus musculus",string=name)))
Tcell_deredundancy_all.distinct.blast.hosthuman <- Tcell_deredundancy_all.distinct.blast %>%
  filter(!is.na(str_match(pattern = "Homo sapiens",string=name)))
DeredundancyAlllist <- list(DistinctBeforeBlast = Tcell_deredundancy_all.distinct,
                            DistinctAfterBlast = Tcell_deredundancy_all.distinct.blast,
                            DistinctAfterBlastHostMus = Tcell_deredundancy_all.distinct.blast.hostMus,
                            DistinctAfterBlastHostHuman = Tcell_deredundancy_all.distinct.blast.hosthuman,
                            DistinctAfterBlastDetailInfo = all)
save(DeredundancyAlllist,file = "./result/yjj/DeredundancyAlllist.RData")
blast.all <- Tcell_deredundancy_all.distinct.blast %>%
  select(description,allele_name,judge.final)
blast.all.human <- Tcell_deredundancy_all.distinct.blast.hosthuman %>%
  select(description,allele_name,judge.final)
blast.all.mouse <- Tcell_deredundancy_all.distinct.blast.hostMus %>%
  select(description,allele_name,judge.final)

write_tsv(blast.all,file = "result/yjj/Tcell_deredundancy_all.clean.tsv",
          quote = "none",col_names = T)
write_tsv(blast.all.human,file = "result/yjj/Tcell_deredundancy_all.human.clean.tsv",
          quote = "none",col_names = T)
write_tsv(blast.all.mouse,file = "result/yjj/Tcell_deredundancy_all.mouse.clean.tsv",
          quote = "none",col_names = T)
