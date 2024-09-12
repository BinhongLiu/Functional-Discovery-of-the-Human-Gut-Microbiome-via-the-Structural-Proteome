library(bio3d)
library(tibble)
library(stringr)
library(magrittr)

PDBscore <- function(dir0){
  filename <- list.files(dir0,pattern = "\\.pdb$",full.names = TRUE)
  # filename <- filename[str_detect(filename,"-F1-model_v")]
  
  af2fold <- data.frame(matrix(nrow = length(filename),ncol = 2))
  for (i in 1:length(filename)) {
    # i <- 1
    id <- gsub("^.*/", "", filename[i])
    # read in
    pdb0 <- read.pdb(filename[i],ATOM.only=T)
    pdb0.atom <- pdb0$atom %>%
      group_by(resno) %>%
      dplyr::slice(1)
    
    af2fold[i,1] <- id
    af2fold[i,2] <- mean(pdb0.atom$b)
  }
  colnames(af2fold) <- c("Uniprot_acc","score")
  return(af2fold)
  
}
# xy <- PDBscore("phage/Typephages_protein")