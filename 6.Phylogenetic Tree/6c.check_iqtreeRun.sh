#this script is run to confirm if IQ-TREE is already performed for every gene

for aln in mafft_output/*.fa; do
    prefix=$(basename "$aln" .fa)
    outdir="iqtree_output/${prefix}"

    # Expected 9 IQ-TREE output files
    expected_files=(
        "${outdir}.bionj"
        "${outdir}.ckp.gz"
        "${outdir}.contree"
        "${outdir}.iqtree"
        "${outdir}.log"
        "${outdir}.midst"
        "${outdir}.model.gz"
        "${outdir}.splits.nex"
        "${outdir}.treefile"
    )

    ## Count how many exist
    count_existing=0
    for f in "${expected_files[@]}"; do
        [[ -f "$f" ]] && ((count_existing++))
    done

    ## If all 9 output files exist, skip
    if [[ $count_existing -eq 9 ]]; then
        echo "✅ Skipping ${prefix} — all output files found."
        continue
    fi

    echo "🚀 Running IQ-TREE for ${prefix}..."
    iqtree -s "$aln" -m MFP -nt 1 -bnni -bb 1000 -pre "$outdir"

    echo "✅ Finished ${prefix}"
done
