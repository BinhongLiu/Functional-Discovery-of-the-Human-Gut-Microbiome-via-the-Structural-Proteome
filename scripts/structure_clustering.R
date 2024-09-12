library(data.table)
library(dplyr)
library(stringr)
library(openxlsx)
library(tibble)
library(RColorBrewer)
library(stats)

TMs <- fread("phage/TMscroes.txt", header = F) %>%
  mutate(V1=gsub("^.*/", "", V1),V2=gsub("^.*/", "", V2)) %>% # split and get the last element
  mutate(V1=str_split_fixed(V1,fixed("."),2)[,1],
         V2=str_split_fixed(V2,fixed("."),2)[,1]) %>%
  dplyr::rename(TM1="V3",TM2="V4",RMSD="V5",ID1="V6",ID2="V7",IDali="V8",L1="V9",L2="V10",Lali="V11") %>%
  mutate(TM=ifelse(L1>L2,round(1-TM1,4) ,round(1-TM2,4))) %>%
  reshape2::dcast(V1 ~ V2,value.var="TM") %>%
  column_to_rownames(var="V1")

hc <- hclust(as.dist(TMs),method = "average")
clusters <- cutree(hc, h = 0.5)
clusters.info <- data.frame(id=names(clusters),cluster_id=clusters) %>%
  group_by(cluster_id) %>%
  mutate(Spike_cc=sum(str_detect(id,"Spikes#")),query_cc=sum(str_detect(id,"query#"))) %>%
  filter(Spike_cc>0 & query_cc>0)
write.table(clusters.info,"phage/endolysin/clusters.txt",sep="\t", col.names = T, row.names = F, quote=FALSE)

