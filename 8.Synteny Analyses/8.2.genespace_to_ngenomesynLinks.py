import pandas as pd

# Load GeneSpace file
df = pd.read_csv("syntenicRegion_coordinates.csv")
df = df[df["genome1"] != df["genome2"]]

# Select relevant columns
links = df[[
    "chr1", "startBp1", "endBp1",
    "chr2", "startBp2", "endBp2"
]]

# Save as GenomeSyn link file (space-separated, no header)
links.to_csv("genomesyn.links", sep="\t", index=False, header=False)

#and later grep the species u want to create one-on-one links
