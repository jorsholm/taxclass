# Instructions to download raw data and generate COI datasets

## BOLD

Downloaded public BOLD snapshot `BOLD_Public.29-Mar-2024.tar.gz`
from [http://doi.org/10.5883/DP-BOLD_Public.29-Mar-2024] on 21 May 2024

## GBOL (test)

Accessed [https://data.bolgermany.de/gbol1/ergebnisse/results] on 2023-11-22
filters:
 - Sub-/Pyhulum[sic] : Arthropoda
 - Collected: 11
 - With Barcode: 11
Export textfile (tsv) -> `data_bolgermany_de_gbol1-search_results-guest-2023-11-22184603.csv`

I _believe_ this means I downloaded data for all specimens in phylum Arthropoda, but only up to 11 per species, but also those 11 should have barcodes.

Extract BOLD process-ID (column 24, removing web api address if present) and barcode sequence (column 27) from GBOL for all specimens with a process ID and >= 600bp COI or COI-5P barcode; then use the process id to look up the taxonomy from class to species (fields 17-23) from BOLD, and output everything into FASTA format:

```sh
awk -F"\t" '
  BEGIN{file="GBOL"};
  file=="GBOL" &&
  length($24)>0 &&
  length($27)>=600 &&
  $26 ~ /COI(-5P)?$/{
    sub(/http.*processid=/, "", $24);
    seq[$24] = $27
  };
  ENDFILE{file="BOLD"};
  file=="BOLD" &&
  $22 != "None" &&
  $23 != "None" &&
  $16=="Arthropoda" &&
  ($1 in seq) &&
  !($23 ~ /sp[.]|aff[.]|cf[.]|nr[.]|agg[.]|t[.]|cluster/) {
    for (i=17;i<=21;i++) {
      if ($i == "None") {
        $i = "dummy_" $(i-1);
        sub(/dummy_dummy/, "dummy", $i)
      }
    }
    sub(/ var\..*/, "", $23) 
    gsub(/ /, "_", $23)
    taxon=$17
    for (i=18; i<=23; i++) taxon = taxon "|" $i
    if (taxon ~ /[0-9]/) next
    print ">" $1, taxon
    print seq[$1]
    delete seq[$1]
  }
' <(zcat data_bolgermany_de_gbol1-search_results-guest-2023-11-22184603.csv.gz)\
 <(tar -xOf BOLD_Public.29-Mar-2024.tar.gz BOLD_Public.29-Mar-2024.tsv)\
 >test_gbol_raw.fasta
```

This yielded 27200 sequences with species ID.

```sh
awk -F"[ |]" '/^>/{
  for (i=2;i<=NF;i++) c[i-1][$i]++
};
END{
  taxon[1]="classes";
  taxon[2]="orders";
  taxon[3]="families";
  taxon[4]="subfamilies";
  taxon[5]="tribes";
  taxon[6]="genera";
  taxon[7]="species";
  for (i=1;i<=7;i++) print "- " length(c[i]) " unique " taxon[i]
}' test_gbol_raw.fasta
```

The taxonomic breakdown is:

- 11 unique classes
- 51 unique orders
- 564 unique families
- 1030 unique subfamilies
- 1504 unique tribes
- 3795 unique genera
- 9031 unique species

## FinBOL (train)

Now, similarly extract the FinBOL sequences from the BOLD public database.  This is easier because they have a BOLD project code, so we don't need to use any other reference file to extract them.

```sh
tar -xOf BOLD_Public.29-Mar-2024.tar.gz BOLD_Public.29-Mar-2024.tsv |
awk -F"\t" '
  $23 != "None" &&\
  $74 ~ /DS-FINPRO/ &&\
  !($23 ~ /sp[.]|aff[.]|cf[.]|nr[.]|agg[.]|t[.]|cluster/) &&\
  $69 ~ /COI(-5P)?$/ &&\
  length($65)>=600 {
    for (i=17;i<=21;i++) {
      if ($i == "None") {
        $i = "dummy_" $(i-1);
        sub(/dummy_dummy/, "dummy", $i)
      }
    }
    gsub(/ /, "_", $23)
    taxon=$17
    for (i=18; i<=23; i++) taxon = taxon "|" $i
    if (taxon ~ /[0-9]/) next
    print ">" $1, taxon
    gsub(/-/, "", $65)
    print $65
}' >train_finbol_raw.fasta
```

36833 sequences

```sh
awk -F"[ |]" '/^>/{
  for (i=2;i<=NF;i++) c[i-1][$i]++
};
END{
  taxon[1]="classes";
  taxon[2]="orders";
  taxon[3]="families";
  taxon[4]="subfamilies";
  taxon[5]="tribes";
  taxon[6]="genera";
  taxon[7]="species";
  for (i=1;i<=7;i++) print "- " length(c[i]) " unique " taxon[i]
}' train_finbol_raw.fasta
```

- 2 unique classes
- 20 unique orders
- 477 unique families
- 916 unique subfamilies
- 1372 unique tribes
- 3895 unique genera
- 11229 unique species

# Alignment via MACSE pipeline

Download [MACSE v2.0.7 jar file](https://www.agap-ge2pop.org/wp-content/uploads/macse/releases/macse_v2.07.jar)
and [MACSE pipelines](https://github.com/ranwez/MACSE_V2_PIPELINES/releases/tag/11.05).

Build singularity container for REPRESENTATIVE_SEQ pipeline from `representative_seqs_v01_sing_3.3.def` and save `representative_seqs_v01_sing_3.3.sif` in this directory.

## Generate representative sequences

I downloaded the repseq of the cox1 CDS from _Drosophila melanogaster_ to use as reference for the alignment.

Then the MACSE pipeline selects 100 sequences from the dataset which align to the reference, and cover the diversity present.

Setting `TMPDIR` and `APPTAINER_TMPDIR` may not be necessary; and would certainly need to be changed for a different system.

```sh
TMPDIR=/home/brfurnea/tmp \
APPTAINER_TMPDIR=/home/brfurnea/tmp \
./representative_seqs_v01_sing_3.3.sif\
  --in_refSeq Drosophila_melanogaster_cox1_refseq.fasta\
  --in_seqFile train_finbol_raw.fasta\
  --in_geneticCode 5\
  --debug
```

## MACSE alignment of representative sequences

```sh
cat Drosophila_melanogaster_cox1_refseq.fasta >>representative_seq_NT.fasta
java -jar macse_v2.07.jar\
 -prog alignSequences\
 -gc_def 5\
 -out_AA finbol_repseq_aln_aa.fasta\
 -out_NT finbol_repseq_aln_nt.fasta\
 -seq representative_seq_NT.fasta\
 -gap_ext_term 0.2\
 -gap_op_term 2.0
```

## MACSE alignment and translation of other sequences

This runs test and train in parallel, but it still takes ~8 hours.

```sh
for data in test_gbol train_finbol; do
  time java -jar macse_v2.07.jar\
    -prog enrichAlignment\
    -align finbol_repseq_aln_nt.fasta\
    -seq ${data}_raw.fasta\
    -out_AA ${data}_raw_aa.fasta\
    -out_NT ${data}_raw_nt.fasta\
    -out_tested_seq_info ${data}_macse_stats.csv\
    -maxSTOP_inSeq 0\
    -maxINS_inSeq 0\
    -fixed_alignment_ON\
    -new_seq_alterable_ON\
    -output_only_added_seq_ON\
    -gc_def 5\
    -gap_ext_term 0.2\
    -gap_op_term 2.0 &
done
```

## Remove gap columns from alignment

```sh
for data in test train; do
  source=$([ $data = test ] && echo gbol || echo finbol)
  sed -i '/^>/!y/!/-/' ${data}_${source}_raw_nt.fasta
  sed -ri -e '/^>/!s/^(-*)!/\1-/' -e '/^>/!s/!(-*)$/-\1/' ${data}_${source}_raw_aa.fasta
  awk '
    BEGIN{
      RS=">"
      FS="\n"
    }
    FNR==1{next}
    FNR==NR{
      seq=""
      for (i=2;i<=NF;i++) seq = seq $i
      seq=gensub(/(-{1,2})([ACGT]{1,2})(-*)$/, "\\2\\1\\3", 1, seq)
      seqlen = length(gensub(/[^ACGT]/, "", "g", seq))
      if (seqlen >= 600) {
        select[$1]=1
        print ">" $1 > "'../data/${data}_nt_aln_label.fasta'"
        print seq > "'../data/${data}_nt_aln_label.fasta'"
      }
      next
    }
    $1 in select {
      seq=""
      for (i=2;i<=NF;i++) seq = seq $i
      print ">" $1 > "'../data/${data}_aa_aln_label.fasta'"
      print seq > "'../data/${data}_aa_aln_label.fasta'"
    }
  ' <(esl-alimask --outformat afa -g ${data}_${source}_raw_nt.fasta)\
    <(esl-alimask --outformat afa -g ${data}_${source}_raw_aa.fasta)
done
```

## Full BOLD taxonomy (only used by PROTAX)

```sh
tar -xOf BOLD_Public.29-Mar-2024.tar.gz BOLD_Public.29-Mar-2024.tsv |
awk -F"\t" '
  $16=="Arthropoda" {
    # remove all placeholder taxa
    for (i=17; i<=23; i++) {
      if ($i ~ /[Uu]nknown|[Uu]nclassified|[Uu]nassigned|[Ii]ncertae|sp[.]|aff[.]|cf[.]|nr[.]|agg[.]|t[.]|cluster|[0-9]|^[a-z]|^[A-Z].*[A-Z]/) $i="None"
      if (i < 23 && $i ~ /[_ ]/) $i = "None"
    }
    # convert "None" to standardized placeholders at ranks above genus
    for (i=17;i<=21;i++) {
      if ($i == "None") {
        $i = "dummy_" $(i-1);
        sub(/dummy_dummy/, "dummy", $i)
      }
    }
    # remove trailing placeholders/None
    for (i=23; i>17; i--) {
      if ($i == "None" && $(i-1) ~ /^dummy_/) $(i-1) = "None"
    }
    if ($17=="None") next
    gsub(/ /, "_", $23);
    taxon=$17
    for (i=18; i<=23; i++) {
      if ($i != "None") taxon = taxon "|" $i
    }
    if (taxon in c) next
    if (taxon~/[0-9]/) next
    c[taxon]=1
    print $1, taxon
  }' >../data/full_tax.txt
```
