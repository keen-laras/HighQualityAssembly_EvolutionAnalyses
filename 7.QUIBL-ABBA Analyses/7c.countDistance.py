import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp

# =====================
# PARAMETERS
# =====================
INPUT = "geneDiscordance.value.txt"
N_PERM = 500   # increase to 1000+ for publication
SEED = 42

np.random.seed(SEED)

# =====================
# LOAD DATA
# =====================
df = pd.read_csv(INPUT, sep="\t", header=None)
df.columns = ["HOG", "GeneID", "Status", "Scaffold", "Start", "End", "Value"]

# clean numeric
df["Start"] = pd.to_numeric(df["Start"], errors="coerce")
df["End"]   = pd.to_numeric(df["End"], errors="coerce")
df = df.dropna(subset=["Start", "End"])

# midpoint
df["Mid"] = (df["Start"] + df["End"]) / 2

# sort genome order
df = df.sort_values(["Scaffold", "Mid"]).reset_index(drop=True)

# =====================
# FUNCTION: get D–D distances
# =====================
def get_dd_distances(dataframe):
    distances = []
    
    for scaffold, group in dataframe.groupby("Scaffold"):
        group = group.sort_values("Mid")
        mids = group["Mid"].values
        status = group["Status"].values
        
        for i in range(len(group) - 1):
            if status[i] == "D" and status[i+1] == "D":
                dist = mids[i+1] - mids[i]
                distances.append(dist)
                
    return distances

# =====================
# OBSERVED
# =====================
obs_distances = get_dd_distances(df)

print(f"Observed D-D pairs: {len(obs_distances)}")

# =====================
# PERMUTATION (NULL)
# =====================
perm_distances_all = []

for i in range(N_PERM):
    shuffled = df.copy()
    shuffled["Status"] = np.random.permutation(shuffled["Status"])
    
    perm_distances = get_dd_distances(shuffled)
    perm_distances_all.extend(perm_distances)

    if (i+1) % 50 == 0:
        print(f"Permutation {i+1}/{N_PERM}")

# subsample null to match observed size
perm_distances_all = np.array(perm_distances_all)

if len(perm_distances_all) > len(obs_distances):
    perm_distances = np.random.choice(perm_distances_all, size=len(obs_distances), replace=False)
else:
    perm_distances = perm_distances_all

# =====================
# STATS
# =====================
stat, pval = ks_2samp(obs_distances, perm_distances)

print("\nKS test:")
print("Statistic =", stat)
print("p-value   =", pval)

# =====================
# PLOTTING (log scale)
# =====================
plt.figure(figsize=(6,5))

# avoid log(0)
obs = np.array(obs_distances)
perm = np.array(perm_distances)

obs = obs[obs > 0]
perm = perm[perm > 0]

bin_width = 10000  # 10 kb

max_val = max(np.max(obs), np.max(perm))
bins = np.arange(0, max_val + bin_width, bin_width)

plt.hist(obs, bins=bins, alpha=0.5, label="Observed (D-D)")
plt.hist(perm, bins=bins, alpha=0.5, label="Null (shuffled)")

#plt.xscale("log")
plt.xlim(0, 1000000) 
plt.xlabel("Distance (bp)")
plt.ylabel("Frequency")
plt.title("Discordant Gene Clustering Test")
plt.legend()

plt.tight_layout()
plt.savefig("discordance_clustering.png", dpi=300)
plt.show()

# =====================
# SAVE OUTPUT
# =====================
pd.DataFrame({"distance": obs}).to_csv("observed_distances.txt", index=False)
pd.DataFrame({"distance": perm}).to_csv("null_distances.txt", index=False)
