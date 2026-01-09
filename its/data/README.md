# Formatting train/test data

The commands shown here generate all variants of the test and train data based
on the files `train_nt_unite.fasta` and `test_nt_unite.fasta`,
which are generated from files in `../data_raw`

## Artificial metabarcoding -- testshort

Generate short test sequences, including only ITS2.
This is accomplished by deleting everything up to the end of 5.8S, as detected
by the Rfam covariance model.

```sh
# micromamba create -n LSUx -c bioconda -c bfurneaux r-lsux
# micromamba activate LSUx
R --vanilla -e '
  seq <- Biostrings::readDNAStringSet("test_nt_unite.fasta")
  pos <- inferrnal::cmsearch(inferrnal::cm_5_8S(), seq, cpu = 20)
  pos <- pos[pos$inc == "!",]
  seq[pos$target_name] <- Biostrings::subseq(seq[pos$target_name], pos$seq_to + 1L)
  Biostrings::writeXStringSet(seq, "testshort_nt_unite.fasta", width = 10000L)
'
```
## Alternate annotations

### Unlabeled sequences

```sh
for f in *_unite.fasta
do
  sed '/^>/s/[|].*//' $f >${f%_unite.fasta}.fasta
done
```

### Sintax-style labels

Sintax requires taxonomic annotations with single-letter codes representing the
primary Linnean ranks.

```sh
sed -r '/^>/{
  s/[|]k__/;tax=k:/
  s/;([pcofgs])__/,\1:/g
}'\
 train_nt_unite.fasta\
 >train_nt_sintax.fasta
```

### "Basic" labels

```sh
tr '|' ' ' <train_nt_unite.fasta |
sed -r '/^>/ {
  s/k__//
  s/;([pcofgs])__/|/g
 }'\
 >train_nt_label.fasta
```

### "RDP" labels

```
tr '\t|' ' ;' <train_nt_label.fasta |
sed 's/ / root;/' >train_nt_rdp.fasta
```

### `raxtax` labels

`raxtax` labels are also similar to Sintax labels, but do not have rank indicators and
a semicolon at the end is required.

```sh
sed 's/[kpcofgst]://g; /^>/s/.*$/&;/' train_nt_sintax.fasta >train_nt_raxtax.fasta
```

## ID-taxonomy maps (without sequences)

Many of these formats use two different delimiters; one to delimit the ID from
the taxonomic classification, and one to delimit taxa at different ranks in the
classification.

### Basic (space+pipe)

```sh
for type in test train
do
  tr '|' ' ' <${type}_nt_unite.fasta |
  sed -nr '/^>/ {
    s/^>//
    s/k__//
    s/;([pcofgs])__/|/g
    p
   }'\
   >${type}_tax.txt
done
```

### TSV with headers

This is the format used by Dnabarcoder. It is also convenient to load into R.

```sh
for f in *_tax.txt
do
  echo -e "id\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies" > ${f%.txt}.tsv
  tr ' |' '\t\t' <$f >>${f%.txt}.tsv
done
```

### Gappa (tab+semicolon)

Only needed for the training data.

```sh
tr ' |' '\t;' <train_tax.txt >train_gappa.tax
```

### Protax (tab+comma, with depth)

Only needed for the training data. This works for the seq2tax file, not the
reference taxonomy.

```sh
tr ' |' '\t,' <train_tax.txt |
sed 's/\t/\t7\t/' >train_protax.tax
```
