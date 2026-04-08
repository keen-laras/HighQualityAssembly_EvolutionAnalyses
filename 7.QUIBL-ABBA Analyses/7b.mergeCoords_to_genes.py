# merge_hog_data.py

# input files
file1 = "index.txt"          # coord.txt
file2 = "node.tsv"           # tree concordance
file3 = "coords.ortho.txt"   # GeneID → Chr Start End

# load HOG → GeneIDE
hog_to_gene = {}
with open(file1) as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 2:
            hog, gene = parts[0], parts[1]
            hog_to_gene[hog] = gene

# load HOG → Status
hog_to_status = {}
with open(file2) as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 2:
            hog, status = parts[0], parts[1]
            hog_to_status[hog] = status

# load GeneID → Coordinates
gene_to_coord = {}
with open(file3) as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 4:
            gene, chr_, start, end = parts[0], parts[1], parts[2], parts[3]
            gene_to_coord[gene] = (chr_, start, end)

# output
with open("genediscordance.txt", "w") as out:
    out.write("HOG_ID\tGeneID\tStatus\tChr\tStart\tEnd\n")

    for hog in hog_to_gene:
        gene = hog_to_gene[hog]

        # only keep if all info exists
        if hog in hog_to_status and gene in gene_to_coord:
            status = hog_to_status[hog]
            chr_, start, end = gene_to_coord[gene]

            out.write(f"{hog}\t{gene}\t{status}\t{chr_}\t{start}\t{end}\n")
