#!/usr/bin/env bash
TIME="$(which time) --verbose"
$TIME makeblastdb -in ../../data/finbol-gbol/train_finbol-gbol_nt_sintax.fasta\
                  -dbtype nucl\
                  -out train_finbol-gbol_nt
$TIME blastn -query ../../data/finbol-gbol/test_finbol-gbol_nt_sintax.fasta\
             -db train_finbol-gbol_nt -out blast_top_hit_finbol_gbol_nt.raw\
             -strand plus\
             -num_threads 8\
             -outfmt "6 delim=\t qseqid sseqid"\
             -max_target_seqs 1\
             -qcov_hsp_perc 80
echo "ID	class	Prob_class	order	Prob_order	family	Prob_family	subfamily	Prob_subfamily	tribe	Prob_tribe	genus	Prob_genus	species	Prob_species" >blast_top_hit_finbol_gbol_nt.txt
sed -r 's/([^;]+);[^;]+;tax=/\1/; s/[kpcofgs]:([^,]+),?/\t\1\t1.0/g' blast_top_hit_finbol_gbol_nt.raw >>blast_top_hit_finbol_gbol_nt.txt

$TIME blastn -query ../../data/finbol-gbol/testshort_finbol-gbol_nt_sintax.fasta\
             -db train_finbol-gbol_nt\
             -out blast_top_hit_finbol_gbol_nt_short.raw\
             -strand plus\
             -num_threads 8\
             -outfmt "6 delim=; qseqid sseqid"\
             -max_target_seqs 1\
             -qcov_hsp_perc 80
echo "ID	class	Prob_class	order	Prob_order	family	Prob_family	subfamily	Prob_subfamily	tribe	Prob_tribe	genus	Prob_genus	species	Prob_species" >blast_top_hit_finbol_gbol_nt_short.txt
sed -r 's/([^;]+);[^;]+;tax=/\1/; s/[kpcofgs]:([^,]+),?/\t\1\t1.0/g' blast_top_hit_finbol_gbol_nt_short.raw >>blast_top_hit_finbol_gbol_nt_short.txt

$TIME makeblastdb -in ../../data/finbol-gbol/train_finbol-gbol_aa_sintax.fasta\
                  -dbtype prot\
                  -out train_finbol-gbol_aa
$TIME blastn -query ../../data/finbol-gbol/test_finbol-gbol_aa_sintax.fasta
             -db train_finbol-gbol_aa\
             -out blast_top_hit_finbol_gbol_aa.raw\
             -num_threads 8\
             -outfmt "6 delim=\t qseqid sseqid"\
             -max_target_seqs 1\
             -qcov_hsp_perc 80
echo "ID	class	Prob_class	order	Prob_order	family	Prob_family	subfamily	Prob_subfamily	tribe	Prob_tribe	genus	Prob_genus	species	Prob_species" >blast_top_hit_finbol_gbol_aa.txt
sed -r 's/([^;]+);[^;]+;tax=/\1/; s/[kpcofgs]:([^,]+),?/\t\1\t1.0/g' blast_top_hit_finbol_gbol_aa.raw >>blast_top_hit_finbol_gbol_aa.txt
