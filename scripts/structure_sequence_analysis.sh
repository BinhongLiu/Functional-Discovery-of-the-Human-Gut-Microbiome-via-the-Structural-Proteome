#######################
# structure comparision and clustering
#########################

################
# 1. primary screening of similar structures performed use FoldSeek
foldseek easy-search \
      /home/liuhongbin/data/protein_prediction/Final_structures/proPhages/ \
      /home/liuhongbin/database/FoldSeek/swissprot_pdb_v4/ \
      swiss.m8 tmpFolder \
      --exhaustive-search 1 \
      -s 7.5 \
      --seq-id-mode 2 \
      --format-output query,target,pident,nident,qcov,qlen,tcov,tlen,evalue,qtmscore,ttmscore,alntmscore,lddt,prob

################
# 2. structural alignment performed use USalign
~/soft/USalign-master/USalign AF-Q877J3-F1-model_v3.cif.gz \
    -dir2 database/ \
    -split 0 \
    /home/liuhongbin/database/pdb100_zelixir/pdb_BC100_list -outfmt 2 >TMscroes.txt

# 3. prepare the distance matrix (1 - TM-score) for the structure-based clustering
# after pair-wised structural aligning, (together with well-annotated templates)
Rscript structure_clustering.R
# the clusters could be visualized in Cytoscape

# 4. structural tree visualizaiton
# a. prepare the input file of MEGA
Rscript MEGA_input.R
# b. import input distance matrix file into the MEGA11 to build the tree file with Unweighted Pair Group Method with Arithmetic mean (UPGMA)
# c. tree visualizaiton using ggtree in R
Rscript structure_tree_visualization.R




#######################
# structure domain segmentation
#########################
# get the domain boundary
python ~/soft/chainsaw/get_predictions.py --structure_directory chainsaw_input_pdbs --output chainsaw.tsv
# segment the mother PDB into domain pdb files
Rscript segment_into_domain_pdbs.R




#######################
# other structure analysis
#########################
# get the pLDDT value of AF2 models
Rscript get_AF2_pLDDT.R



#######################
# sequence similarity analysis
#########################
makeblastdb -dbtype prot -in template.faa -out template
blastp -query query.FA -db template -evalue 1000 -max_target_seqs 5 \
      -max_hsps 1 -num_threads 160 \
      -outfmt "6 qseqid sseqid evalue nident qlen slen pident qcovs" \
      -out query.template.blastp
# get the sequence identity normalized by the longer sequence
echo "query,template,identity" |sed "s/,/\t/g" >query.template.blastp.identity
awk '{if($5<$6) {$9=$4/$6} else {$9=$4/$5}; print $1"\t"$2"\t"$9}' query.template.blastp >>query.template.blastp.identity

