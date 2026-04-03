## 5. Genome Annotation
This repository contains a genome annotation pipeline that integrates homology-based with [BRAKER (version 3.0.8)](https://github.com/Gaius-Augustus/BRAKER) and ab initio gene prediction with [EviAnn (version 2.0.4)](https://github.com/alekseyzimin/EviAnn_release) to produce a high-confidence gene set.

- Run both software

        `bash run_eviann.sh` and `bash run_braker.sh`

- Merge both pipeline's output

        `bash run_mergeOutput.sh`

