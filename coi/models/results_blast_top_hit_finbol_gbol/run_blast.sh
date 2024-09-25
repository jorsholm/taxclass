#!/usr/bin/env bash
#SBATCH --job-name=blast_top_hit
#SBATCH --account=project_2005718
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=4800M
#SBATCH --cpus-per-task=4
#SBATCH --array=4
#SBATCH --output=blast_%a.out
#SBATCH --error=blast_%a.out
#SBATCH --mail-type=ALL

module load blast
export PATH="/appl/opt/time/1.9/bin:$PATH"

TIME="$(which time) --verbose"
$TIME makeblastdb -in ../../data/finbol-gbol/train_finbol-gbol_nt_sintax.fasta\
                  -dbtype nucl\
                  -out train_finbol-gbol_nt_$SLURM_ARRAY_TASK_ID
$TIME blastn -query ../../data/finbol-gbol/test_finbol-gbol_nt_sintax.fasta\
             -db train_finbol-gbol_nt_$SLURM_ARRAY_TASK_ID\
             -out blast_finbol_gbol_nt_$SLURM_ARRAY_TASK_ID.raw\
             -strand plus\
             -num_threads $SLURM_ARRAY_TASK_ID\
             -outfmt "6 delim=; qseqid sseqid pident"\
             -max_target_seqs 1\
             -qcov_hsp_perc 80
echo "ID	class	Prob_class	order	Prob_order	"\
     "family	Prob_family	subfamily	Prob_subfamily	"\
     "tribe	Prob_tribe	genus	Prob_genus	species	Prob_species"\
     >blast_top_hit_finbol_gbol_nt_$SLURM_ARRAY_TASK_ID.txt
cp blast_top_hit_finbol_gbol_nt_$SLURM_ARRAY_TASK_ID.txt\
   blast_thresh_finbol_gbol_nt_$SLURM_ARRAY_TASK_ID.txt

sed -r 's/([^;]+);[^;]+;tax=/\1/; s/[kpcofgs]:([^,]+),?/\t\1\t1.0/g'\
    blast_finbol_gbol_nt_$SLURM_ARRAY_TASK_ID.raw\
    >>blast_top_hit_finbol_gbol_nt_$SLURM_ARRAY_TASK_ID.txt

awk -F"[,;]"    'BEGIN {
       n=split("90,93,95,96,96.5,97,98", thresh)
     }
     NR>1 {printf "%s", $1
       sub(/tax=/, "")
       for (i = 1; i <= n; i++) {
         if ($10 >= thresh[i]) {
           sub(/[kpcofgs]:/, "", $(i + 2))
           printf "\t%s\t%s", $(i + 2), "1.0"
         } else {
           printf "\tNA\tNA"
         }}
       print ""
     }'\
     blast_finbol_gbol_nt_$SLURM_ARRAY_TASK_ID.raw\
     >>blast_thresh_finbol_gbol_nt_$SLURM_ARRAY_TASK_ID.txt

$TIME blastn -query ../../data/finbol-gbol/testshort_finbol-gbol_nt_sintax.fasta\
             -db train_finbol-gbol_nt_$SLURM_ARRAY_TASK_ID\
             -out blast_finbol_gbol_nt_short_$SLURM_ARRAY_TASK_ID.raw\
             -strand plus\
             -num_threads $SLURM_ARRAY_TASK_ID\
             -outfmt "6 delim=; qseqid sseqid pident"\
             -max_target_seqs 1\
             -qcov_hsp_perc 80
echo "ID	class	Prob_class	order	Prob_order	"\
     "family	Prob_family	subfamily	Prob_subfamily	"\
     "tribe	Prob_tribe	genus	Prob_genus	species	Prob_species"\
     >blast_top_hit_finbol_gbol_nt_short_$SLURM_ARRAY_TASK_ID.txt
cp blast_top_hit_finbol_gbol_nt_short_$SLURM_ARRAY_TASK_ID.txt\
   blast_thresh_finbol_gbol_nt_short_$SLURM_ARRAY_TASK_ID.txt
sed -r 's/([^;]+);[^;]+;tax=/\1/; s/;[0-9.]+$//; s/[kpcofgs]:([^,]+),?/\t\1\t1.0/g'\
     blast_finbol_gbol_nt_short_$SLURM_ARRAY_TASK_ID.raw\
     >>blast_top_hit_finbol_gbol_nt_short_$SLURM_ARRAY_TASK_ID.txt

awk -F"[,;]"    'BEGIN {
       n=split("90,93,95,96,96.5,97,98", thresh)
     }
     NR>1 {printf "%s", $1
       sub(/tax=/, "")
       for (i = 1; i <= n; i++) {
         if ($10 >= thresh[i]) {
           sub(/[kpcofgs]:/, "", $(i + 2))
           printf "\t%s\t%s", $(i + 2), "1.0"
         } else {
           printf "\tNA\tNA"
         }}
       print ""
     }'\
     blast_finbol_gbol_nt_short_$SLURM_ARRAY_TASK_ID.raw\
     >>blast_thresh_finbol_gbol_nt_short_$SLURM_ARRAY_TASK_ID.txt

$TIME makeblastdb -in ../../data/finbol-gbol/train_finbol-gbol_aa_sintax.fasta\
                  -dbtype prot\
                  -out train_finbol-gbol_aa_$SLURM_ARRAY_TASK_ID
$TIME blastp -query ../../data/finbol-gbol/test_finbol-gbol_aa.fasta\
             -db train_finbol-gbol_aa_$SLURM_ARRAY_TASK_ID\
             -out blast_finbol_gbol_aa_$SLURM_ARRAY_TASK_ID.raw\
             -num_threads $SLURM_ARRAY_TASK_ID\
             -outfmt "6 delim=; qseqid sseqid"\
             -max_target_seqs 1\
             -qcov_hsp_perc 80
echo "ID	class	Prob_class	order	Prob_order	"\
     "family	Prob_family	subfamily	Prob_subfamily	"\
     "tribe	Prob_tribe	genus	Prob_genus	species	Prob_species"\
     >blast_top_hit_finbol_gbol_aa_$SLURM_ARRAY_TASK_ID.txt
sed -r 's/([^;]+);[^;]+;tax=/\1/; s/[kpcofgs]:([^,]+),?/\t\1\t1.0/g'\
     blast_finbol_gbol_aa_$SLURM_ARRAY_TASK_ID.raw\
     >>blast_top_hit_finbol_gbol_aa_$SLURM_ARRAY_TASK_ID.txt
