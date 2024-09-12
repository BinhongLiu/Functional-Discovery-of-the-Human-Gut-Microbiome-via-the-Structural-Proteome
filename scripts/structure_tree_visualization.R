library(dplyr)
library(magrittr)
library(ggnewscale)
library(UniprotR)
library(data.table)
library(tibble)
library(ggplot2)
library(aplot)
library(ggtree)
library(stringr)
library(treeio)
library(openxlsx)
library(rstatix)
library(tidyr)



##
# endolysin rep tree
######
tm <- read.tree("phage/endolysin/tree_new/Newick Export.nwk")
# tm <- read.tree("phage/endolysin/tree_new/align.tree.nwk")
tm$tip.label

endoActives <- read.xlsx("templates.xlsx",sheet = "Endolysin_result") %>%
  filter(activity=="yes") %>%
  select(protein)
tree.domains3 <- read.xlsx("Enzyme.xlsx",sheet = "endolysin_domains") %>%
  filter(is.na(mark)) %>%
  distinct(global_id,domain_info,.keep_all = T) %>%
  arrange(global_id) %>%
  group_by(global_id,domain_type) %>%
  mutate(nn=row_number()) %>%
  mutate(domain_type=paste(domain_type,nn,sep = "_")) %>%
  reshape2::dcast(source+global_id+species+gram~domain_type,value.var= "domain_info") %>%
  mutate(activity=ifelse(global_id %in% endoActives$protein,"yes","other")) %>%
  left_join(fread("phage/metacerberusAll.annotation.txt") %>% select(locus_tag,seqREP_ano),by=c("global_id"="locus_tag")) %>%
  mutate(seqREP_ano=ifelse(source=="Spikes","seqREP_yes",seqREP_ano))
tree.domains3_info <- tree.domains3 %>% column_to_rownames(var = "global_id")  %>%
  mutate(Catalytic_1=ifelse(is.na(Catalytic_1),"other",Catalytic_1),
         Catalytic_2=ifelse(is.na(Catalytic_2),"other",Catalytic_2),
         CBD_1=ifelse(is.na(CBD_1),"other",CBD_1),
         unknown_1=ifelse(is.na(unknown_1),"other",unknown_1))

difDA <- setdiff(tm$tip.label,tree.domains3$global_id)
tm <- tm %>% drop.tip(difDA)

tm %<>% dplyr::left_join(tree.domains3,by=c('label'='global_id'))

p11 <- ggtree(tm, # layout='circular',
              # layout="fan",open.angle = 30,
              size=.2
) +
  # hexpand(.1) +
  hexpand(8) +
  geom_tippoint(aes(color=source,shape=source),size=1.7) +
  # geom_text(aes(label=node)) + # get the node nubmers that need to be flip
  geom_tiplab(aes(x = x+.8),size = 1.5,offset = 0.015) +
  scale_shape_manual(values = c("query"=16,"Spikes"=17)) +
  scale_color_manual(values = c("query"="black","Spikes"="#2E9739")) +
  new_scale_color() +
  geom_tiplab(aes(label=species,color=gram,x = x+.4),fontface="italic",size = 2.2,offset = 0.015) +
  scale_color_manual(values = c("G+"="#E64B35CC","G-"="#4DBBD5CC"))+
  new_scale_color() +
  # scale_color_manual(values = c("positive"="black","negative"="grey"))
  geom_tippoint(aes(color=activity,x = x+.05),size=.5,shape=16) +
  scale_color_manual(values = c("yes"="red","other"="white"))

# p11 <- flip(p11, 102, 70) %>% flip(103, 116) %>% flip(87, 71)%>% flip(88,92)%>% flip(26,94) %>% flip(25,95)%>% flip(23,24) %>% flip(56,126)

p11 <- flip(p11, 60,76) %>% flip(77, 81) %>% flip(92, 94) %>% flip(100,95) %>% flip(26,83) %>% flip(25,84) %>% flip(106,102) %>% flip(49,48) %>% flip(23,24)

p12 <- gheatmap(p11, tree.domains3_info[,"Catalytic_1", drop=F], offset=0.08, width=0.18,legend_title="Catalytic",colnames_angle=0, colnames_offset_y = .25) +
  scale_x_ggtree() +
  scale_fill_manual(name="Catalytic",values = my_col)
# scale_fill_viridis_d(option="D", name="discrete\nvalue")
p13 <- p12 + new_scale_fill()
p14 <- gheatmap(p13, tree.domains3_info[,"Catalytic_2", drop=F], offset=.15, width=0.18,legend_title="Catalytic_2",colnames_angle=0, colnames_offset_y = .25) +
  scale_x_ggtree() +
  scale_fill_d3(name="Catalytic_2")
p15 <- p14 + new_scale_fill()
p16 <- gheatmap(p15, tree.domains3_info[,"CBD_1", drop=F], offset=.25, width=0.18,legend_title="CBD_1",colnames_angle=0, colnames_offset_y = .25) +
  scale_x_ggtree() +
  scale_fill_d3(name="CBD_1")
p17 <- p16 + new_scale_fill()
p18 <- gheatmap(p17, tree.domains3_info[,"unknown_1", drop=F], offset=.32, width=0.18,legend_title="unknown_1",colnames_angle=0, colnames_offset_y = .25) +
  scale_x_ggtree() +
  scale_fill_jama(name="unknown_1")
p18
# rotate_tree(p14, 90)
ggsave("phage/endolysin/tree_new/Global_struct_tree7.pdf",width = 9,height = 7,limitsize = FALSE)