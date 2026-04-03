import pandas as pd

# ======================
# 1. Load GeneSpace data
# ======================
df = pd.read_csv("syntenicBlock_coordinates.csv")

# ======================
# Clean + remove selfBlk
# ======================
df["genome1"] = df["genome1"].str.strip()
df["genome2"] = df["genome2"].str.strip()

df = df[df["genome1"] != df["genome2"]]

# Ensure numeric coordinates
for col in ["startBp1","endBp1","startBp2","endBp2"]:
    df[col] = pd.to_numeric(df[col], errors="coerce")

# ======================
# 2. Build per-species table
# ======================
df1 = df[["genome1","chr1","startBp1","endBp1"]].rename(
    columns={
        "genome1":"species",
        "chr1":"chromosome",
        "startBp1":"start",
        "endBp1":"end"
    }
)

df2 = df[["genome2","chr2","startBp2","endBp2"]].rename(
    columns={
        "genome2":"species",
        "chr2":"chromosome",
        "startBp2":"start",
        "endBp2":"end"
    }
)

combined = pd.concat([df1, df2])

# Remove invalid intervals
combined = combined[combined["end"] > combined["start"]]

# ======================
# 3. Merge overlapping regions
# ======================
def merge_intervals(df):
    merged_lengths = []

    for (species, chrom), group in df.groupby(["species", "chromosome"]):
        intervals = group[["start", "end"]].sort_values("start").values

        merged = []
        for start, end in intervals:
            if not merged or start > merged[-1][1]:
                merged.append([start, end])
            else:
                merged[-1][1] = max(merged[-1][1], end)

        total_len = sum(e - s for s, e in merged)

        merged_lengths.append([species, chrom, total_len])

    return pd.DataFrame(merged_lengths, columns=["species","chromosome","syntenic_length"])

summary = merge_intervals(combined)

# ======================
# 4. Load chromosome lengths
# ======================
len_df = pd.read_csv("all_sum.len", sep="\t", header=None)
len_df.columns = ["chromosome","chr_length"]

# ======================
# 5. Merge + calculate coverage
# ======================
summary = summary.merge(len_df, on="chromosome", how="left")

# Check missing lengths
print("⚠️ Missing chr_length:", summary["chr_length"].isna().sum())

summary["coverage_pct"] = (summary["syntenic_length"] / summary["chr_length"]) * 100
summary["coverage_pct"] = summary["coverage_pct"].clip(upper=100)

# ======================
# 6. Final format
# ======================
summary = summary[[
    "species",
    "chromosome",
    "chr_length",
    "syntenic_length",
    "coverage_pct"
]]

summary.to_csv("x.final_synteny_coverage.tsv", sep="\t", index=False)

print("✅ Final table ready!")
