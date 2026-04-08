## 7. ILS - Introgression Analyses
## A. QuIBL - ABBA analyses
This analyses is done to investigate trees that might undergo incomplete lineage sorting (ILS) or introgression during its evolution line.

    bash 7.ILS_introgression.sh
    
## B. Discordant tree analyses
This analyses is done to compute the contribution of discordance genes to each nodes

    # Create index file "Orthogroup_HOGID   Gene_Name"
        awk '
        > FNR==1 {fname=FILENAME; sub(/\.fa$/,"",fname)}
        > /^>/ && /{species}/ {
        >     gsub(/^>/,"")
        >     print fname "\t" $0
        > }
        > ' /path/to/OrthoFinder_Output/Single_Copy_Orthologue_Sequences/*.fa > index.txt 
    
    python 7b.quartet_rf_classify.py --confi 7b.node.confi --gene "iqtree_output/*tree" --out node.tsv
    python 7b.mergeCoords_to_genes.py
    python 7b.addValue.py 
