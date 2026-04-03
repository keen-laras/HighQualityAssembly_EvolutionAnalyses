import sys

input_gff = sys.argv[1]
output_gff = sys.argv[2]

with open(input_gff) as f, open(output_gff, "w") as out:
    out.write("##gff-version 3\n")

    for line in f:
        if line.startswith("#") or line.strip() == "":
            continue

        parts = line.strip().split("\t")
        if len(parts) < 9:
            continue

        chrom, src, feature, start, end, score, strand, phase, attr = parts

        # GENE
        if feature == "gene":
            gene_id = attr.strip()
            out_attr = f"ID={gene_id}"
            out.write("\t".join([chrom, src, feature, start, end, score, strand, phase, out_attr]) + "\n")

        # TRANSCRIPT
        elif feature == "transcript":
            tid = attr.strip()
            gid = tid.split(".")[0]
            out_attr = f"ID={tid};Parent={gid}"
            out.write("\t".join([chrom, src, feature, start, end, score, strand, phase, out_attr]) + "\n")

        # EXON / CDS / others
        else:
            # extract transcript_id "g6.t1"
            if 'transcript_id "' in attr:
                tid = attr.split('transcript_id "')[1].split('"')[0]
            else:
                continue

            out_attr = f"Parent={tid}"
            out.write("\t".join([chrom, src, feature, start, end, score, strand, phase, out_attr]) + "\n")
