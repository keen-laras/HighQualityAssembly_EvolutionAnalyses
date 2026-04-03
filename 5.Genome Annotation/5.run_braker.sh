export PATH="/home/stereonote/miniconda3/bin:${PATH}"
source /home/stereonote/miniconda3/etc/profile.d/conda.sh
source /home/stereonote/miniconda3/bin/activate /home/stereonote/miniconda3/envs/braker3
cd /workdir/
perl /home/stereonote/BRAKER/scripts/braker.pl --threads 88 --genome=/path/to/{masked_input}.fa \
--prot_seq=/path/to/cleaned/proteins.fa
