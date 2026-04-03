#/bin/bash
main=""
nowdir=`pwd`
Gffpath="$main/data/gff"
##prepared data 
mkdir -p Syntenyblock/bed
cat $main/species.list|awk '{print $1"\t"$2}'|while read fullname consname 
do
        awk -F "[\t=;]" '{if($3=="mRNA")print "'$consname_'"$1"\t"$4-1"\t"$5"\t'$consname'_"$10}' $Gffpath/$consname.gff >Syntenyblock/bed/$consname.bed
done

mkdir Syntenyblock/orthofinder
Orthologs="$main/Orthologs"
ln -s $Orthologs/pep Syntenyblock/peptide
for f in $Orthologs/Results_* ;do ln -s $f Syntenyblock/orthofinder ;done 

## grogram path
orthofinder_path="" #https://github.com/davidemms/OrthoFinder
MCScanX_path="" #https://github.com/wyp1125/MCScanX
# install https://github.com/jtlovell/GENESPACE 

echo '''library(GENESPACE)
setwd("'$nowdir'")
workdir <- "'$nowdir'/Syntenyblock"
path2mcscanx <- "'$MCScanX_path'"

existingOrthofinderDir <- "'$nowdir'/Syntenyblock/orthofinder"


genomeIDs <- c("Adia", "Cruf", "Etat", "Xuni", "Ajin", "Pmen", "Talbo", "Gshe", "Bful", "Csep", "Ofor", "Pgut") #c("Bful", "Csep", "Ofor", "Pgut")

ploidy <- rep(1, 12)

gpar <- init_genespace(
    genomeIDs = genomeIDs,
    ploidy = ploidy,
    wd = workdir,
    nCores = 1,
    path2mcscanx = path2mcscanx,
    path2orthofinder = "'$orthofinder_path'",
    rawOrthofinderDir = existingOrthofinderDir,
    blkSize = 5,
    synBuff = 100
)

out <- run_genespace(
    gsParam = gpar,
    overwrite = FALSE,
    overwriteBed = FALSE,
    overwriteSynHits = TRUE,
    overwriteInBlkOF = TRUE
)
''' >GENESPACE.r

Rscript GENESPACE.r
