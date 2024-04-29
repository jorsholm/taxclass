#!/usr/bin/env bash
sed '/^>/!s/-//g' ../../data/finbol-gbol/train_finbol-gbol_genus_sintax.fasta >train_finbol-gbol_genus_sintax.fasta
sed '/^>/!s/-//g' ../../data/finbol-gbol/train_finbol-gbol_species_sintax.fasta >train_finbol-gbol_species_sintax.fasta
sed '/^>/!s/-//g' ../../data/finbol-gbol/test_finbol-gbol_genus.fasta >test_finbol-gbol_genus.fasta
sed '/^>/!s/-//g' ../../data/finbol-gbol/test_finbol-gbol_species.fasta >test_finbol-gbol_species.fasta
sed '/^>/!s/-//g' ../../data/finbol-gbol/test_finbol-gbol_genus.fasta >testshort_finbol-gbol_genus.fasta
sed '/^>/!s/-//g' ../../data/finbol-gbol/test_finbol-gbol_species.fasta >testshort_finbol-gbol_species.fasta
TIME="$(which time) --verbose"
#$TIME vsearch --makeudb_usearch train_finbol-gbol_genus_sintax.fasta --output train_finbol-gbol_genus.udb
#$TIME vsearch --sintax test_finbol-gbol_genus.fasta --db train_finbol-gbol_genus.udb --tabbedout sintax_finbol_gbol_genus.raw --threads 8
echo "ID	class	Prob_class	order	Prob_order	family	Prob_family	subfamily	Prob_subfamily	genus	Prob_genus" >sintax_finbol_gbol_genus.txt
sed 's/ [^\t]*\t//; s/)*,*[pcofgs]:/\t/g; s/(/\t/g; s/)\t+.*//' sintax_finbol_gbol_genus.raw >>sintax_finbol_gbol_genus.txt

#$TIME vsearch --sintax testshort_finbol-gbol_genus.fasta --db train_finbol-gbol_genus.udb --tabbedout sintax_finbol_gbol_genus_short.raw --threads 8
echo "ID	class	Prob_class	order	Prob_order	family	Prob_family	subfamily	Prob_subfamily	genus	Prob_genus" >sintax_finbol_gbol_genus_short.txt
sed 's/ [^\t]*\t//; s/)*,*[pcofgs]:/\t/g; s/(/\t/g; s/)\t+.*//' sintax_finbol_gbol_genus_short.raw >>sintax_finbol_gbol_genus_short.txt


#$TIME vsearch --makeudb_usearch train_finbol-gbol_species_sintax.fasta --output train_finbol-gbol_species.udb
#$TIME vsearch --sintax test_finbol-gbol_species.fasta --db train_finbol-gbol_species.udb --tabbedout sintax_finbol_gbol_species.raw --threads 8
echo "ID	class	Prob_class	order	Prob_order	family	Prob_family	subfamily	Prob_subfamily	genus	Prob_genus	species	Prob_species" >sintax_finbol_gbol_species.txt
sed 's/ [^\t]*\t//; s/)*,*[pcofgs]:/\t/g; s/(/\t/g; s/)\t+.*//' sintax_finbol_gbol_species.raw >>sintax_finbol_gbol_species.txt

#$TIME vsearch --sintax testshort_finbol-gbol_species.fasta --db train_finbol-gbol_species.udb --tabbedout sintax_finbol_gbol_species_short.raw --threads 8
echo "ID	class	Prob_class	order	Prob_order	family	Prob_family	subfamily	Prob_subfamily	genus	Prob_genus	species	Prob_species" >sintax_finbol_gbol_species_short.txt
sed 's/ [^\t]*\t//; s/)*,*[pcofgs]:/\t/g; s/(/\t/g; s/)\t+.*//' sintax_finbol_gbol_species_short.raw >>sintax_finbol_gbol_species_short.txt
