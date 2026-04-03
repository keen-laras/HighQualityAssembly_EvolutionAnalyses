## 8. Synteny Analyses
### A. Generating syntenic blocks with [GENESPACE (version 1.2.3)](https://github.com/jtlovell/GENESPACE)
Syntenic blocks are obtained from OrthoFinder output [(see 6.Phylogenetic Tree)](https://github.com/keen-laras/HighQualityAssembly_EvolutionAnalyses/tree/main/6.Phylogenetic%20Tree)

    bash 8a.GENESPACE.sh

### B. Visualization using [NGenomeSyn (version 1.43)](https://github.com/hewm2008/NGenomeSyn)
This visualization is for chromosome-level assemblies only

    python 8b.genespace_to_ngenomesynLinks.py
    NGenomeSyn-1.43/bin/NGenomeSyn -InConfi config.txt -OutPut output
    

    
