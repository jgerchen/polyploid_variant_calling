module add picard-2.8.1
set +eu
source /storage/brno2/home/gerchenj/mambaforge-pypy3/bin/activate
conda activate gatk4
set -eu
GATK4=/storage/brno2/home/gerchenj/software/gatk-4.3.0.0/gatk
module load bcftools
