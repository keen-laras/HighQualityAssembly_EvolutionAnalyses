cd /workdir/
export PATH=/opt/software/minimap2-2.20/:/opt/software/hisat2/:$PATH
/opt/software/EviAnn-2.0.2/bin/gffread -W -y {species_1}_protein.faa -w {species_1}_rna_from_genomic.fna -g {species_1}_genomic.fna {species_1}_genomic.gff
/opt/software/EviAnn-2.0.2/bin/gffread -W -y {species_2}_protein.faa -w {species_2}_rna.fna -g {species_2}_genomic.fna {species_2}_genomic.gff
cat *protein* > proteins.faa
cat *rna* > transcripts.faa
/opt/software/EviAnn-2.0.2/bin/eviann.sh -t 64 -g /path/to/{input}.fa -e /path/to/transcripts.faa -p /path/to/proteins.fa -l

# clean output IDs for BRAKER input
sed's/ .*//; s/:/_/g; s/|/-/g' /path/to/proteins.faa > cleaned.proteins.faa
