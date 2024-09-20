#!/usr/bin/env bash
TIME="$(which time) --verbose"
$TIME vsearch --makeudb_usearch ../../data/finbol-gbol/train_finbol-gbol_nt_sintax.fasta --output train_finbol-gbol_nt.udb
$TIME vsearch --sintax ../../data/finbol-gbol/test_finbol-gbol_nt_sintax.fasta --db train_finbol-gbol_nt.udb --tabbedout sintax_finbol_gbol_nt.raw --threads 8
echo "ID	class	Prob_class	order	Prob_order	family	Prob_family	subfamily	Prob_subfamily	tribe	Prob_tribe	genus	Prob_genus	species	Prob_species" >sintax_finbol_gbol_nt.txt
sed -r 's/ [^\t]*\t//; s/[kpcofgs]:([^,\t]+)\(([01][.][0-9]+)\),?/\t\1\t\2/g; s/\t\+.*//' sintax_finbol_gbol_nt.raw >>sintax_finbol_gbol_nt.txt

$TIME vsearch --sintax ../../data/finbol-gbol/testshort_finbol-gbol_nt_sintax.fasta --db train_finbol-gbol_nt.udb --tabbedout sintax_finbol_gbol_nt_short.raw --threads 8
echo "ID	class	Prob_class	order	Prob_order	family	Prob_family	subfamily	Prob_subfamily	tribe	Prob_tribe	genus	Prob_genus	species	Prob_species" >sintax_finbol_gbol_nt_short.txt
sed -r 's/ [^\t]*\t//; s/[kpcofgs]:([^,\t]+)\(([01][.][0-9]+)\),?/\t\1\t\2/g; s/\t\+.*//' sintax_finbol_gbol_nt_short.raw >>sintax_finbol_gbol_nt_short.txt
