#!/usr/bin/env bash
#SBATCH --job-name=blast_top_hit
#SBATCH --account=project_2005718
#SBATCH --partition=small
#SBATCH --mem-per-cpu=4800M
#SBATCH --cpus-per-task=40
#SBATCH --array=40
#SBATCH --output=blast_top_hit_%a.out
#SBATCH --error=blast_top_hit_%a.out
#SBATCH --mail-type=ALL

TIME="$(which time) --verbose"
$TIME makeblastdb -in ../../data/finbol-gbol/train_finbol-gbol_nt_sintax.fasta\
                  -dbtype nucl\
                  -out train_finbol-gbol_nt_$ARRAY_JOB_TASK_ID
$TIME blastn -query ../../data/finbol-gbol/test_finbol-gbol_nt_sintax.fasta\
             -db train_finbol-gbol_nt_$ARRAY_JOB_TASK_ID\
             -out blast_top_hit_finbol_gbol_nt_$ARRAY_JOB_TASK_ID.raw\
             -strand plus\
             -num_threads 8\
             -outfmt "6 delim=\t qseqid sseqid"\
             -max_target_seqs 1\
             -qcov_hsp_perc 80
echo "ID	class	Prob_class	order	Prob_order	"\
     "family	Prob_family	subfamily	Prob_subfamily	"\
     "tribe	Prob_tribe	genus	Prob_genus	species	Prob_species"\
     >blast_top_hit_finbol_gbol_nt_$ARRAY_JOB_TASK_ID.txt
sed -r 's/([^;]+);[^;]+;tax=/\1/; s/[kpcofgs]:([^,]+),?/\t\1\t1.0/g'\
    blast_top_hit_finbol_gbol_nt_$ARRAY_JOB_TASK_ID.raw\
    >>blast_top_hit_finbol_gbol_nt_$ARRAY_JOB_TASK_ID.txt

$TIME blastn -query ../../data/finbol-gbol/testshort_finbol-gbol_nt_sintax.fasta\
             -db train_finbol-gbol_nt_$ARRAY_JOB_TASK_ID\
             -out blast_top_hit_finbol_gbol_nt_short_$ARRAY_JOB_TASK_ID.raw\
             -strand plus\
             -num_threads 8\
             -outfmt "6 delim=; qseqid sseqid"\
             -max_target_seqs 1\
             -qcov_hsp_perc 80
echo "ID	class	Prob_class	order	Prob_order	"\
     "family	Prob_family	subfamily	Prob_subfamily	"\
     "tribe	Prob_tribe	genus	Prob_genus	species	Prob_species"\
     >blast_top_hit_finbol_gbol_nt_short_$ARRAY_JOB_TASK_ID.txt
sed -r 's/([^;]+);[^;]+;tax=/\1/; s/[kpcofgs]:([^,]+),?/\t\1\t1.0/g'\
     blast_top_hit_finbol_gbol_nt_short_$ARRAY_JOB_TASK_ID.raw\
     >>blast_top_hit_finbol_gbol_nt_short_$ARRAY_JOB_TASK_ID.txt

$TIME makeblastdb -in ../../data/finbol-gbol/train_finbol-gbol_aa_sintax.fasta\
                  -dbtype prot\
                  -out train_finbol-gbol_aa_$ARRAY_JOB_TASK_ID
$TIME blastp -query ../../data/finbol-gbol/test_finbol-gbol_aa_sintax.fasta
             -db train_finbol-gbol_aa_$ARRAY_JOB_TASK_ID\
             -out blast_top_hit_finbol_gbol_aa_$ARRAY_JOB_TASK_ID.raw\
             -num_threads 8\
             -outfmt "6 delim=\t qseqid sseqid"\
             -max_target_seqs 1\
             -qcov_hsp_perc 80
echo "ID	class	Prob_class	order	Prob_order	"\
     "family	Prob_family	subfamily	Prob_subfamily	"\
     "tribe	Prob_tribe	genus	Prob_genus	species	Prob_species"\
     >blast_top_hit_finbol_gbol_aa_$ARRAY_JOB_TASK_ID.txt
sed -r 's/([^;]+);[^;]+;tax=/\1/; s/[kpcofgs]:([^,]+),?/\t\1\t1.0/g'\
     blast_top_hit_finbol_gbol_aa_$ARRAY_JOB_TASK_ID.raw\
     >>blast_top_hit_finbol_gbol_aa_$ARRAY_JOB_TASK_ID.txt
