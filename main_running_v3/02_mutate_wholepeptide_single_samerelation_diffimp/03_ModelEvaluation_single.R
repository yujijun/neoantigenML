# BASIC INFO ---------------------------
##
## Script name:06_ModelEvaluation.R
##
## Purpose of script: This script is for model evaluation from 00_main.r
##
## Author: JijunYu
##
## Date Created: 2021-12-11
## Update Date:
##
## Copyright (c) Jijun Yu, 2021
## Email: jijunyuedu@outlook.com
##
## Notes:
##
#--------------Main -----------------------
#### 01_load model result ####
output_path <- "./result_v3/03_mutate_wholepeptide_double_diffrelation_result/03_机器学习计算结果/"
Allfiles <- list.files(path = output_path,include.dirs = T,recursive = T,full.names = T)
Alltsv <- grep(Allfiles,pattern = "tsv",value = T)
eval.list <- list()
n =1
for(i in Alltsv){
  df_tmp <- readr::read_delim(i)
  base_name <- basename(i)
  base_name <- str_remove(base_name,pattern = "_model_evaluation.tsv")
  df_tmp$project <- base_name
  eval.list[[n]] <- df_tmp
  n = n+1
}
eval.all <- do.call(rbind.data.frame,eval.list)
eval.all  %>%
  as.data.table() %>%
  mutate(method = stringr::str_remove(outer_task_id,pattern = "neoML[.]")) %>%
  mutate(Pro_method = stringr::str_c(project,"_",method)) %>%
  mutate(project = fct_relevel(project,"mutateML_100_1160_120",after = Inf)) %>%
  mutate(project = as.factor(project)) %>%
  ggplot(aes(x = Pro_method,y = outer_classif.auc,fill = task_id)) +
  geom_boxplot(alpha=0.4) +
  facet_grid(.~project,scales = "free_x") +
  theme_classic() +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x = element_text(angle = 90,hjust = 1,size = 8,face = "bold",family = "ArialMT")) +
  xlab("Methods") +
  ylab("AUC value") -> plot1
ggsave(filename = paste0(output_path,"AllAuc.pdf"),
       width = 45,height = 20,units = "cm")




