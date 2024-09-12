library(bio3d)
library(tibble)
library(stringr)
library(magrittr)

## prepare SRs domains
chainsaw <- fread("phage/endolysin/tree_new/chainsaw.SRs.tsv")

for (i in 1:nrow(chainsaw)) {
  # i <- 12
  pdbx <- read.pdb(paste0("phage/endolysin/tree_new/tree_pdbs/",chainsaw$chain_id[i],".pdb"))
  ndom <- chainsaw$ndom[i]
  if(ndom>0){
    doms <- str_split(chainsaw$chopping[i],fixed(","))[[1]]
    for (j in 1:ndom) {
      # j <- 1
      if(str_detect(doms[j],"_",negate = T)) {
        r.start <- str_split(doms[j],"-(?=[^-]*$)")[[1]][1]
        r.end <- str_split(doms[j],"-(?=[^-]*$)")[[1]][2]
        domain.sele <- atom.select(pdbx,# "calpha",
                                   resno=r.start:r.end)
        domain.pdb <- trim.pdb(pdbx, domain.sele)
        write.pdb(domain.pdb,file = paste0("phage/endolysin/tree_new/chainsaw.SRs/",chainsaw$chain_id[i],"_chainsaw_",j,".pdb"))
      } else {
        segs0 <- str_split(doms[j],fixed("_"))[[1]]
        segs <- c()
        for (k in 1:length(segs0)) { # k <- 1
          sk <- str_split(segs0[k],"-(?=[^-]*$)")[[1]]
          segs <- c(segs,sk[1]:sk[2])}
        
        domain.sele <- atom.select(pdbx,# "calpha",
                                   resno=segs)
        domain.pdb <- trim.pdb(pdbx, domain.sele)
        write.pdb(domain.pdb,file = paste0("phage/endolysin/tree_new/chainsaw.SRs/",chainsaw$chain_id[i],"_chainsaw_",j,".pdb"))
      }
    }
  } else { file.copy(from = paste0("phage/endolysin/tree_new/tree_pdbs/",chainsaw$chain_id[i],".pdb"),
                     to = paste0("phage/endolysin/tree_new/chainsaw.SRs/",chainsaw$chain_id[i],"_chainsaw_0.pdb"))}
  
}
