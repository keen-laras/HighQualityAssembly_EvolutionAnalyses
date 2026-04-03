# ANNOTATE REPEAT BASED ON PUBLIC LIBRARY

conda deactivate
source ~/.bashrc
source activate RepeatMasker
cd /workdir/
run RepeatMasker

# pick between your library options
# OPTION 1: based on NCBI's database, make sure the '-species' can be found in NCBI taxonomy key (eg. mouse, serpentes, vertebrates, Homo sapiens)
RepeatMasker -xsmall -gff -species serpentes -a -pa 88 -dir  repeatmasker/ /path/to/{input}.fa 

# OPTION 2: based on published library (dfam or repbase)
RepeatMasker -xsmall -gff -lib /path/to/{database_lib} -a -pa 88 -dir repeatmasker/ /path/to/{input}.fa
