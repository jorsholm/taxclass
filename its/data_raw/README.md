---
editor_options: 
  markdown: 
    wrap: 72
---

# Instructions to download raw data and generate ITS datasets

## Unite

Download [Unite 10.0 general fasta release for Fungi
2](https://doi.org/10.15156/BIO/2959332), saving as default name
`sh_general_release_04.04.2024.tgz`, which includes the sequences we
will use as test sequences, along with their taxonomy.

Also download [Unite 10.0 Full Unite+INSD dataset for
Fungi](https://doi.org/10.15156/BIO/2959330) and save as default name
`UNITE_public_21.04.2024.fasta.gz`. This will be used to lookup the SH
identification for the training sequences.

Finally, download the Unite SH training data; this will be used to
lookup taxonomy for SH. Code from Unite
[sh_matching_pub](https://github.com/TU-NHM/sh_matching_pub/blob/c40762ba571ff84d69e4fc5492f36263389d2185/sh_matching.def#L60C5-L61C73)
git repository.

```         
wget https://s3.hpc.ut.ee/plutof-public/original/5dcbe93a-b50c-4ff5-8898-528160c4e593.zip
mv 5dcbe93a-b50c-4ff5-8898-528160c4e593.zip sh_matching_data_0_5.zip
```

## CBG/Westerdijk Institute yeast and filamentous datasets

Download [PopSet
1102638064](https://www.ncbi.nlm.nih.gov/popset/1102638064) as fasta,
save as `Westerdijk_yeast.fasta`.

Download [PopSet
1473241625](https://www.ncbi.nlm.nih.gov/popset/1473241625) as fasta,
save as `Westerdijk_filamentous.fasta`.

## Test data - Unite RefSeq

These are the sequences in Unite which have been specifically chosen by
a curator to represent species hypotheses (SH) so we expect them to be
more accurately identified than other sequences.

``` sh
tar -xzOf sh_general_release_s_04.04.2024.tgz sh_general_release_dynamic_s_04.04.2024.fasta |
grep -A1 --no-group-separator "refs" \
 >unite_refs.fasta
```

This yields 18895 sequences.

## Train data -- Westerdijk Institute yeast and filamentous datasets

These are from a large, well-curated culture collection, which are very
likely to be correctly identified. Many type strains are included, which
are correctly identified by definition, as long as the species that they are
the type for is not synonymized with an older species. For many species
there are several representatives. However, because the dataset is
culture-based, it does not include any fungi which cannot be cultured
using present techniques.

First we combine the two datasets and remove any sequences which are not
identified to a single, validated species.

``` sh
gawk '
  BEGIN{
    FS="|"
    file = 1
    all_seq = 0
    missing_accno = 0
    missing_sh = 0
    multi_taxon = 0
    multitax_seqs = 0
    wrote_seq = 0
  }
  file==1 && $1~/^>/ && !($1 ~ /^>UDB/) {
    sub(/^>/, "", $1)
    sh[$1]=$3
  }
  file==2 {
    if (!($2 ~ /unspecified/)) {
      taxonomy[$1]=$2
    }
  }
  ENDFILE{
    file++
    if (file == 2) {
      FS="\t"
    }
    if (file == 3) {
      FS="\n"
      OFS="|"
      RS=">"
    }
  }
  file>2 && FNR==1 { next }
  file>2 {
    all_seq++
    species=gensub(/^[A-Z_]+[0-9.]+ (.+) culture.+/, "\\1", 1, $1)
    accno=gensub(/^([A-Z_]+[0-9]+).+/, "\\1", 1, $1)
    myseq=""
    for (i=2; i<=NF; i++) myseq = myseq $i
    sub(/>/, "", myseq)
    if (! (accno in sh)) {
       print "missing accno: " accno >"/dev/stderr"
       missing_accno++
       next
    }
    if (! (sh[accno] in taxonomy)) {
      print "missing taxonomy for \"" sh[accno] "\" for accno \"" accno >"/dev/stderr"
      missing_sh++
      next
    }
    found_tax++
    map[species][taxonomy[sh[accno]]]++
    seq[accno] = myseq
    sp[accno] = species
  }
  END {
    for (species in map) {
      if (length(map[species]) > 1) {
        print "multiple taxonomies for \"" species "\": " >"/dev/stderr"
        for (taxon in map[species]) {
          print "  " taxon >"/dev/stderr"
          multitax_seqs += map[species][taxon]
        }
        multi_taxon++
      }
    }
    for (accno in seq) {
      if (length(map[sp[accno]][taxonomy[sh[accno]]]) == 1) {
        print ">" accno, taxonomy[sh[accno]]
        print seq[accno]
        wrote_seq++
      }
    }
    print "total sequences: " all_seq >"/dev/stderr"
    print "sequences with accno not found: " missing_accno >"/dev/stderr"
    print "sequences with SH taxonomy not found: " missing_sh >"/dev/stderr"
    print "sequences belonging to inconsistent species: " multitax_seqs >"/dev/stderr"
    print "sequences written: " wrote_seq >"/dev/stderr"
    print "number of mapped species: " length(map) >"/dev/stderr"
    print "species with inconsistent taxonomy: " multi_taxon >"/dev/stderr"
  }
  ' <(zcat UNITE_public_21.04.2024.fasta.gz)\
  <(unzip -p sh_matching_data_0_5.zip data/shs_out.txt)\
  Westerdijk_yeast.fasta\
  Westerdijk_filamentous.fasta\
  >train_westerdijk_raw.fasta
```

### Remove overlap between test and train

Test sequences are almost all "singletons" because they are the single
representative sequence of their SH. Thus, when a test sequence is removed, its
entire species is removed.
(An exception is when a single species name is used by more than one species
hypothesis.)
However the training set contains an mean of 2.75 sequences per species name.

Thus, our first choice when a sequence occurs in both sets is to remove it from
the training set, if there is another sequence from the same species.

On the other hand, if it is the only sequence for that species in both test and
train, then we retain it in train.  There are already plenty of sequences in
test which belong to species missing from train.


```sh
awk -F"[|]" '
  FNR==NR && /^>/ {
    last_accno=$1
    if (/_sp$/) next
    train++
    sub(/^>/, "", $1)
    last_accno=$1
    train_accno[$1]=1
    train_tax[$1]=$2
    train_count[$2]++
    next
  }
  FNR==NR {
  train_seq[last_accno] = train_seq[last_accno] $1
  }
  FNR!=NR && /^>/ {
    last_accno=$2
    if (/_sp$/) next
    test++
    test_accno[$2]=1
    test_tax[$2]=$5
    test_count[$5]++
    if ($2 in train_accno) {
      overlap++
      if (train_count[$5] == 1) {
        if (test_count[$5] == 1) {
         single_overlap++
        } else if (test_count[$5] == 2) {
         single_overlap--
        }
      }
    }
    next
  }
  FNR!=NR {
  test_seq[last_accno]=test_seq[last_accno] $1
  }
  END {
    print test " test sequences in " length(test_count) " species"
    print train " train sequences in " length(train_count) " species"
    print overlap " overlapping"
    print single_overlap " overlapping singletons"
    for (accno in train_tax) {
      if (accno in test_tax) {
        if (train_count[train_tax[accno]] == 1) {
          print ">" accno "|" train_tax[accno] "\n" train_seq[accno] >"train_nt_raw.fasta"
        } 
      } else {
        print ">" accno "|" train_tax[accno] "\n" train_seq[accno] >"train_nt_raw.fasta"
      }
    }
    for (accno in test_tax) {
      if (accno in train_tax) {
        if (train_count[test_tax[accno]] > 1) {
          print ">" accno "|" test_tax[accno] "\n" test_seq[accno] >"test_nt_raw.fasta"
        }
      } else {
        print ">" accno "|" test_tax[accno] "\n" test_seq[accno] >"test_nt_raw.fasta"
      }
    }
  }'\
  train_westerdijk_raw.fasta\
  unite_refs.fasta
```

### Extract ITS region

```sh
# micromamba create -n ITSx -c bioconda itsx
# micromamba activate ITSx

for f in test train
do
  ITSx -i ${f}_nt_raw.fasta\
       -o ${f}_nt\
       -t F\
       --nhmmer T\
       --partial 50\
       --save_regions none\
       --complement F\
       --cpu 20\
       --preserve T
done
```

In many of these cases, the end of SSU or the beginning of LSU
are not detected, because the fragment is too short. This is common in SSU
because the standard ITS1 primer is very close to the end of SSU, and the
remaining piece is too short to be confidently detected by the HMMs used by
ITSx. Commonly used reverse primers in LSU are farther from the beginning of
LSU, but this can still be a problem if the sequence is trimmed due to low
quality, and additionally the beginning of LSU is more variable than the end of
SSU and so more difficult to identify from a short fragment.

For non-detected SSU, we use `cutadapt` to trim sequences matching the end of
SSU, from just before the beginning of the ITS1 primer site to the end, with
some common variations in fungi included.

```sh
#micromamba create -n cutadapt -c bioconda cutadapt=4.2
#micromamba activate cutadapt

for f in test train
do
  gawk '
    BEGIN{
      FS="\t"
    }
    FNR == NR {
      if ($3 ~ /SSU: Not found/) {
        missing_ssu[$1] = 1
      } else {
        found_ssu[$1] = 1
      }
      next
    }
    ENDFILE {
      RS=">"
      FS="\n"
      ORS=""
    }
    FNR==1 {
      next
    }
    {
      if ($1 in missing_ssu) {
        print ">" $0 >"'"${f}_nt_no_ssu.fasta"'"
        delete missing_ssu[$1]
      } else if ($1 in found_ssu) {
        print ">" $0 >"'"${f}_nt_ssu_trimmed.fasta"'"
        delete found_ssu[$1]
      }
    }'\
    ${f}_nt.positions.txt\
    ${f}_nt.full_and_partial.fasta\
    ${f}_nt.full.fasta

  cutadapt -g GGTYTCCWAGGTGAACCWGCGGARGGATCATTN\
           -e 0.2\
           ${f}_nt_no_ssu.fasta\
           >>${f}_nt_ssu_trimmed.fasta
done
```

For non-detected LSU, we use LSUx, an approach similar to ITSx but which uses a
longer covariance model (CM) instead of a pair of HMMs, enabling it to detect
the end of 5.8S and beginning of LSU based in part on complementary base-pairing
between the two.

```sh
# micromamba create -n LSUx -c bioconda -c bfurneaux r-lsux
# micromamba activate LSUx
for f in test train
do
  gawk '
    BEGIN{
      FS="\t"
    }
    FNR == NR {
      if ($7 ~ /LSU: Not found/) {
        missing_lsu[$1] = 1
      } else {
        found_lsu[$1] = 1
      }
      next
    }
    ENDFILE {
      RS=">"
      FS="\n"
      ORS=""
    }
    FNR==1 {
      next
    }
    {
      if ($1 in missing_lsu) {
        print ">" $0 >"'"${f}_nt_no_lsu.fasta"'"
        delete missing_lsu[$1]
      } else if ($1 in found_lsu) {
        print ">" $0 >"'"${f}_nt_lsu_trimmed.fasta"'"
        delete found_lsu[$1]
      }
    }'\
    ${f}_nt.positions.txt\
    ${f}_nt_ssu_trimmed.fasta
    
  R --vanilla -e 'seq <- Biostrings::readDNAStringSet("'${f}'_nt_no_lsu.fasta")
    pos <- LSUx::lsux(seq, cpu = 20)
    pos <- pos[pos$region == "ITS2",]
    seq[pos$seq_id] <- IRanges::narrow(seq[pos$seq_id], end = pos$end)
    Biostrings::writeXStringSet(seq, "'${f}'_nt_lsu_trimmed.fasta", append = TRUE, width = 10000)'
done
```

There exist some sequences which, probably due to sequencing errors resulting in
divergent flanking sequences, are not trimmed by any of these methods. These are
unlikely to interfere with any of the assignment algorithms, and are likely also
to be present in real-world reference databases, so we retain them.

### Format "Incertae Sedis

In the fungal classification we are using only "primary" ranks, so "Incertae
sedis" represents actual taxonomic uncertainty, rather than simply groups where
"secondary" ranks are not in use. Thus it is not safe to assume that two
different cases of "Incertae sedis" within the same ancestor taxon are
equivalent.  To distinguish them, we append the highest rank enclosed taxon,
which is often the genus, to the "Incertae sedis" label.

```sh
mkdir -p ../data
for f in test train;
do
  gawk -F'[;]' '
    BEGIN{
      OFS=";"
    }
    /Incertae_sedis/ {
      named_taxon = gensub(/^g__/, "", 1, $6)
      for (rank = 5; rank > 1; rank--) {
        if ($rank ~ /Incertae_sedis$/) {
          $rank = $rank "_" named_taxon
        } else {
          named_taxon = gensub(/^[pcofg]__/, "", 1, $rank)
        }
      }
    }
    {
      print $0
    }'\
    ${f}_nt_lsu_trimmed.fasta\
    >../data/${f}_nt_unite.fasta
done
```

#### Remaining sequences:

```sh
gawk '
  BEGIN{
    FS="[|;\n]"
    RS=">"
  }
  FNR==NR {
    if (FNR==1) next
    train++
    for (rank=1; rank <=7; rank++) {
      taxon[rank][$(rank+1)] = 1
    }
    next
  }
  FNR==1{
    next
  }
  {
    test++
    for (rank=1; rank <=7; rank++) {
      if ($(rank+1) in taxon[rank]) {
        taxon[rank][$(rank+1)]=0
      } else {
        taxon[rank][$(rank+1)]=-1
      }
    }
  }
  END {
    print "sequences in train: " train " in test: " test
    for (rank=1; rank <=7; rank++) {
      train_only = 0
      test_only = 0
      shared = 0
      for (tax in taxon[rank]) {
        if (taxon[rank][tax] == 1) {
          train_only++
        } else if (taxon[rank][tax] == 0) {
          shared++
        } else {
          test_only++
        }
      }
      print "rank " rank " train only: " train_only " shared: " shared " test only: " test_only
    }
  }'\
  ../data/train_nt_unite.fasta\
  ../data/test_nt_unite.fasta
```

## Full Unite taxonomy -- for Protax only

```sh
unzip -p sh_matching_data_0_5.zip data/shs_out.txt |
gawk -F'[;\t]' '
  BEGIN {
    rank[2] = "kgd"
    rank[3] = "phy"
    rank[4] = "cl"
    rank[5] = "ord"
    rank[6] = "fam"
    rank[7] = "gen"
    rank[8] = "sp"
  }
  {
    # remove rank prefixes and placeholder taxa
    for (i=2; i <= 8; i++) {
      sub(/^[kpcofgs]__/, "", $i)
      if ($i ~ /_sp$/) $i="None"
      if ($i ~ /unspecified$/) $i = "None"
    }
    # remove non-fungi
    if ($2 != "Fungi") next
    # remove incomplete taxonomy; Protax does not deal with "dangling" taxa
    # (i.e., those who have no children at species level)
    if ($8 == "None") next
    
    # convert "None" to standardized placeholders at ranks above genus
    next_known = $7
    for (i=6; i >= 2; i--) {
      if ($i ~ /Incertae_sedis/) {
        if (next_known != "None") {
          $i = $i "_" next_known
        } else {
          $i = "None"
        }
      } else {
        next_known = $i
      }
    }
    
    # generate taxon name
    taxon=$2
    for (i=3; i <= 8; i++) {
      if ($i != "None") {
        taxon = taxon "|" $i
      }
    }
    if (taxon in c) next
    c[taxon]=1
    print $1, taxon
  }'\
  >../data/full_tax.txt
```
