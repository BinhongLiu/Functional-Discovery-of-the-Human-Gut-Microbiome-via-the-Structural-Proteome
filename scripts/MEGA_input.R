library(data.table)
library(dplyr)
library(stringr)
library(openxlsx)
library(tibble)
library(RColorBrewer)
library(stats)

options(scipen=10)

TMs <- fread("phage/TMscroes.txt", header = F) %>%
  mutate(V1=gsub("^.*/", "", V1),V2=gsub("^.*/", "", V2)) %>% # split and get the last element
  mutate(V1=str_split_fixed(V1,fixed("."),2)[,1],
         V2=str_split_fixed(V2,fixed("."),2)[,1]) %>%
  dplyr::rename(TM1="V3",TM2="V4",RMSD="V5",ID1="V6",ID2="V7",IDali="V8",L1="V9",L2="V10",Lali="V11") %>%
  mutate(TM=ifelse(L1>L2,round(1-TM1,4) ,round(1-TM2,4))) %>%
  reshape2::dcast(V1 ~ V2,value.var="TM") %>%
  column_to_rownames(var="V1")


TMs[upper.tri(TMs,diag = T)] <- " "

names1 <- rownames(TMs)


lin1 <- paste0("#mega","\n","!Title: Concatenated Files;", "\n","!Format DataType=Distance DataFormat=LowerLeft;","\n", sep = '')

write.table(lin1,"phage/endolysin/tree_new/align.tree.meg", col.names = FALSE, row.names = FALSE,quote = FALSE)
write.table(paste('#',names1, sep = ''),"phage/endolysin/tree_new/align.tree.meg", sep="\t",append = TRUE, col.names = FALSE, row.names = FALSE,quote = FALSE)
write.table("\n","phage/endolysin/tree_new/align.tree.meg", sep="\t",append = TRUE, col.names = FALSE, row.names = FALSE,quote = FALSE)
write.table(TMs,"phage/endolysin/tree_new/align.tree.meg", sep="\t", append = TRUE,col.names = FALSE, row.names = FALSE,quote = FALSE)
options(scipen=0)

