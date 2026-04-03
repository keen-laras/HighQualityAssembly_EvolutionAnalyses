## 4. Repeat Annotation
This repository provides a pipeline for identifying repetitive elements and generating a fully masked genome using [RepeatMasker (version 4.2.2)](https://github.com/Dfam-consortium/RepeatMasker) and [BEDTools](https://github.com/arq5x/bedtools2).

**- Run software depending on your library**

    `bash run_denovo.sh` or `bash run_lib.sh`

**- Mask output file**
        `perl extract_lowercase_bed.pl /path/to/repeatmasker/{sample}.fa.masked {output}.nmasked.bed \
        cat /denovo/{output}.nmasked.bed /dfam_lib/{output}.nmasked.bed > {output}_combined_nmasked.raw.bed sort -k1,1 -k2,2n         {output}_combined_nmasked.raw.bed > {output}_combined_nmasked.sorted.bed \
        bedtools/2.29.2/bin/mergeBed -i {output}_combined_nmasked.sorted.bed > {output}_combined_nmasked.merged.bed \ 
        bedtools/2.29.2/bin/maskFastaFromBed -mc N -fi /path/to/{sample}.fa -fo {output}.FINAL.masked.fa -bed {output}_combined_nmasked.merged.bed`
