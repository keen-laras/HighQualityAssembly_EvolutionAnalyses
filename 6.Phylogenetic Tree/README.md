## 6. Phylogenetic Tree

### A. Ortholog prediction with [OrthoFinder (version 3.0.1)](https://github.com/davidemms/OrthoFinder)
One-to-one orthologs **(Single Copy Orthologs)** are the expected output for further analysis in this pipeline.

    `orthofinder \
        -f input_files \
        -M ~{input_method} \
        -S ~{input_sequenceSearchMethod} \
        -o ~{input_fileName} \
        -t ~{input_cpu}`

### B. Orthologs alignment with [MAFFT (version 7.471)](https://mafft.cbrc.jp/alignment/server/index.html)
This software performs gene alignment. We ran MAFFT for every single gene seperately and concatenate them eventually.

    `bash 6b.run_mafft.sh`

_(optional)_ Each sequence ID usually contains something like **>SpeciesA_geneID**. However, for our analysis, we dont need the gene ID because our end goal is to make a species tree. ASTRAL (see step D.) is not able to calculate each species' topology if the sequence ID is not the same for every genes tree (see step C.). That's why this script helps removing all the gene ID from each sequence, resulting in an **>SpeciesA** sequence ID.

    `bash 6b.edit_geneID.sh`

### C. Generating gene trees with [IQ-TREE (version 3.0.1)](https://iqtree.github.io/doc/Tutorial)
We used CDS alignment instead of protein alignment for further analyses, please proceed with your needs. CDS are used to obtain more evolutionary information as it captured both synonymous and nonsynonymous mutation that protein alignments couldn't.

    `python 6c.convert_MAFFT_to_CDS.py`

Run IQ-Tree 

    `iqtree -s {cds.align}.fa -m MFP -bb 1000 -bnni -nt 1 -pre {iqtree_output}`

and check for corrupted trees

    `bash 6c.check_iqtreeRun.sh`

### D. Generating species tree with [ASTRAL (version 5.7.1)](https://github.com/smirarab/ASTRAL)

    `cat {iqtree_output}/*.treefile > all.treefile
    astral -i all.treefile -o all.astral.speciestree > all.speciestree.log 2>&1 &`

### E. Visualize tree
Tree visualization can be done in a tree building website, for example [ITOL](https://itol.embl.de/upload.cgi)

### F. Bootstrap confirmation
Unconfident bootstraps (<100%) are tested with alternative topologies using [DiscoVista (version 1.0)](https://github.com/esayyari/discovista)

    `python DiscoVista/src/utils/discoVista.py -a {annotation}.txt -m 5 -p ./ -g Outgroup -o {output}`


