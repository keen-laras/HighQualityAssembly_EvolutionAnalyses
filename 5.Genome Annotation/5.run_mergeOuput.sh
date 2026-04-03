# Operation
# need to be understood that the evi annotation is homology-based method (we used 3 snakes and 1 chicken), while braker used ab-initio to predict the genes. thus this pipeline is to filtered out which sites on braker's prediction is not found in the evi's prediction. overlapping sites are called out by coverages, and the sites from evi are further chosen as the backbone rather than the braker's (for the overlapping site only). as evi is homology-based, the accuracy of the protein prediction is higher than the braker, but still braker is needed because there are chances that some protein coding genes are specificially found in a species. while evi is used as the backbone for gene annotation, braker's is eventually concatenated to push accuracy and gene profiles.

#1. Filter out protein-coding genes from the EVI GFF annotation file (removing non-coding RNAs and pseudogenes) and convert them to BED format for genomic region analysis
awk '{if($3=="gene")print}' /path/to/{evi_output}.pseudo_label.gff | grep -v lncRNA | grep -v processed_pseudogene |awk -F "[\t;=]" '{print $1"\t"$4-1"\t"$5"\t"$10}' > evi.gene.bed

#2. Extract all gene location information from the breaker's GTF annotation file and convert it to BED format
awk '{if($3=="gene")print $1"\t"$4-1"\t"$5"\t"$NF}' /path/to/{braker_output}.gtf > braker.gene.bed

#3. bedtools coverage: Use the bedtools tool to calculate the coverage relationship between two bed files (braker's bed and evi's bed). -wa to keep the informations in file A, meaning that overlapping sites will not be excluded and it will write based on the braker's annotation
/share/app/bedtools/2.29.2/bin/bedtools coverage -a braker.gene.bed -b evi.gene.bed -wa | awk '{if($NF==0)print $4}' | while read line ;do awk '$4=="'$line'"' braker.gene.bed;done > braker.kept.bed

# 4. get the non-coding genes from evi
awk '{if($3=="gene")print}'  /path/to/{evi_output}.pseudo_label.gff  | grep -E "lncRNA|processed_pseudogene" |awk -F "[\t;=]" '{print $1"\t"$4-1"\t"$5"\t"$10}' > evi.lncRNA.psudo.bed
# Filter the evi.lncRNA.psudo.bed region that overlaps with braker.kept.bed to obtain evi.lncRNA.psudo.bed.kept (retained)
/share/app/bedtools/2.29.2/bin/bedtools coverage -a evi.lncRNA.psudo.bed -b braker.kept.bed -wa |awk '{if($NF==0)print $1"\t"$2"\t"$3"\t"$4}' > evi.lncRNA.psudo.bed.kept

# 5. remove redundants from braker's pepseq 
/share/app/samtools/1.11/bin/samtools faidx /path/to/{braker_output}.aa
awk -F "." '{print $1}'  /path/to/{braker_output}.aa.fai | sort | uniq |while read line ;do grep "$line\.t" /path/to/{braker_output}.aa.fai | sort -k2,2rn | head -1 ;done | awk '{print $1}'| while read ids ; do /share/app/samtools/1.11/bin/samtools faidx /path/to/{braker_output}.aa $ids;done > braker.filter.redundancy.aa 
# list for which IDs should be retained
awk '{print $NF}' braker.kept.bed | while read line ;do grep "${line}.t" braker.filter.redundancy.aa ;done | awk -F ">" '{print $2}' > {input}.list
/share/app/seqkit/0.14.0-dev/seqkit grep -f {input}.list braker.filter.redundancy.aa > braker.kept.prot.fa

# 6. remove redundants from evi's pepseq
/share/app/gffread/0.12.6/gffread /path/to/{evi_output}.pseudo_label.gff -g /path/to/{input}.fasta -x cds.fa -y pep.fa
/share/app/samtools/1.11/bin/samtools faidx pep.fa
awk -F "-" '{print $1}' pep.fa.fai | sort | uniq | while read line;do grep $line pep.fa.fai | sort -k2,2rn | head -1 |awk '{print $1}';done > 1
/share/app/seqkit/0.14.0-dev/seqkit grep -f 1 pep.fa > pep_removed_redundant.fa
cat braker.kept.prot.fa pep_removed_redundant.fa > all.pep.fa

# 7. create a gff file, merging both retained annotation
/share/app/bedtools/2.29.2/bin/bedtools intersect -f 0.99 -a /path/to/{braker_output}.gtf -b braker.kept.bed -wa -u > braker.kept.gff
python 5.cleanup_braker_gtf.py braker.kept.gff braker.kept_fixed.gff

# modify evi's gff
grep ">" pep_removed_redundant.fa | sed 's/>//g' > evi.keptID.txt
awk 'BEGIN{while((getline<"evi.keptID.txt")>0) ids[$1]=1}
     $0 ~ /^#/ {print; next}
     {
         match($0, /ID=([^;]+)/, a);
         match($0, /Parent=([^;]+)/, b);
         if (ids[a[1]] || ids[b[1]]) print
     }' /path/to/{braker_output}.pseudo_label.gff \
     > evi.kept.gff
cat evi.kept.gff braker.kept_fixed.gff > raw.gff
sort -k1,1 -k4,4n raw.gff > final.gff

grep ">" all.pep.fa | awk -F ">" '{print $2}' > keep_id_for_gff.txt
/share/app/gffread/0.12.6/gffread final.gff -o final.fixed.gff
python /path/to/filter_redundant_gff.py keep_id_for_gff.txt final.fixed.gff tmp.FINAL.gff
/share/app/gffread/0.12.6/gffread tmp.FINAL.gff -o {output}.FINAL.struc.gff

# transcript + cds generating
/share/app/gffread/0.12.6/gffread {output}.FINAL.struc.gff -g /path/to/{input}.fasta -w all.transcript.fa -x all.cds.fa

#8. QC
echo 'Number of Protein Coding Genes:' > qc.txt
grep ">" all.pep.fa | wc -l >> qc.txt # number of protein coding genes
awk '$3=="gene"||$3=="mRNA" {len=$5-$4+1; sum+=len; n++} END{print "Average gene length:", sum/n}' final.gff >> qc.txt  
awk '$3=="exon"{gene=$10; gsub(/[";]/,"",gene); exon_count[gene]++; exon_len[gene]+=$5-$4+1}
END{
  for (i in exon_count){
    total_exon_count+=exon_count[i];
    total_exon_len+=exon_len[i];
    n++;
  }
  print "Average exon number per gene:", total_exon_count/n;
  print "Average exon length:", total_exon_len/total_exon_count;
}' final.gff >> qc.txt

awk '
# Capture exon features from both EviAnn and AUGUSTUS
$3 == "exon" {
    id = "";
    if ($9 ~ /Parent=/) {
        match($9, /Parent=([^;]+)/, a);
        id = a[1];
    } else if ($9 ~ /transcript_id/) {
        match($9, /transcript_id "([^"]+)"/, a);
        id = a[1];
    }
    if (id != "") {
        start[id] = start[id]" "$4;
        end[id]   = end[id]" "$5;
        exon_len[id] += ($5 - $4 + 1);
        exon_count[id]++;
    }
}
END {
    for (id in start) {
        split(start[id], s, " ");
        split(end[id], e, " ");
        n_exons = length(s);
        if (n_exons > 1) {
            asort(s);
            asort(e);
            for (i = 2; i <= n_exons; i++) {
                intron_len = s[i] - e[i-1] - 1;
                if (intron_len > 0) {
                    total_intron_len += intron_len;
                    total_intron_count++;
                }
            }
        }
        total_exon_count += exon_count[id];
        total_exon_len += exon_len[id];
        total_transcripts++;
    }
    print "Average exon number per transcript:", total_exon_count / total_transcripts;
    print "Average exon length:", total_exon_len / total_exon_count;
    print "Average intron number per transcript:", total_intron_count / total_transcripts;
    print "Average intron length:", total_intron_len / total_intron_count;
}' final.gff >> qc.txt

awk '
$3 == "CDS" {
    id = "";
    # Extract transcript ID for both annotation styles
    if ($9 ~ /Parent=/) {
        match($9, /Parent=([^;]+)/, a);
        id = a[1];
    } else if ($9 ~ /transcript_id/) {
        match($9, /transcript_id "([^"]+)"/, a);
        id = a[1];
    }
    if (id != "") {
        cds_len[id] += ($5 - $4 + 1);
    }
}
END {
    for (id in cds_len) {
        total_len += cds_len[id];
        n++;
    }
    print "Total transcripts with CDS:", n;
    print "Average CDS length (bp):", total_len / n;
}' final.gff >> qc.txt 

# 9. transcript + cds
/share/app/gffread/0.12.6/gffread braker.kept_fixed.gff -g /path/to/{input}.fasta -w braker_transcripts.fa
/share/app/gffread/0.12.6/gffread braker.kept_fixed.gff -g /path/to/{input}.fasta -x braker_cds.fa
/share/app/gffread/0.12.6/gffread evi.kept.gff -g /path/to/{input}.fasta -w evi_transcripts.fa
/share/app/gffread/0.12.6/gffread evi.kept.gff -g /path/to/{input}.fasta -x evi_cds.fa
cat braker_transcripts.fa  evi_transcripts.fa > all.transcript.fa
cat braker_cds.fa  evi_cds.fa > all.cds.fa
