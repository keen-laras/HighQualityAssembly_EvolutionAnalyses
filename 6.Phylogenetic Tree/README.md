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
  
