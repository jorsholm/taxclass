# Formatting train/test data

The commands shown here generate all variants of the test and train data based
on the files `train_nt_aln_label.fasta`, `train_aa_aln_label.fasta`,
`test_nt_aln_label.fasta` and `test_aa_aln_label.fasta`,
which are generated from files in `../data_raw`

## Artificial metabarcoding -- testshort

Generate simulated short "metabarcoding" test sequences, starting only after
the BF3 primer site, which ends 240bp (80 aa) after the alignment start.
240 gaps are included so that the alignment is still accurate.

```sh
awk '
  BEGIN{
    blank=sprintf("%240s", "");
    gsub(/ /, "-", blank)
  };
  /^>/{print; next};
  {print blank substr($1, 241)}
' \
test_nt_aln_label.fasta \
>testshort_nt_aln_label.fasta


awk '
  BEGIN{
    blank=sprintf("%80s", "");
    gsub(/ /, "-", blank)
  };
  /^>/{print; next};
  {print blank substr($1, 81)}
' \
test_aa_aln_label.fasta \
>testshort_aa_aln_label.fasta
```

## Unaligned sequences

```sh
for f in *_aln_label.fasta
do
  sed '/^>/!s/-//g' $f >${f%aln_label.fasta}label.fasta
done
```

## Alternate annotations

### Unlabeled sequences

```sh
for f in *_label.fasta
do
  sed '/^>/s/ .*//' $f >${f%_label.fasta}.fasta
done
```

### Sintax-style labels

Sintax requires taxonomic annotations with single-letter codes representing the
primary Linnean ranks.  In this case we are using secondary ranks (subfamily
and tribe) but the total number of ranks still fits within the number Sintax
accepts.

```sh
for f in train_*_label.fasta;
do
  awk -F"[ |]" '
    /^>/ {
      print $1 ";tax=k:" $2 ",p:" $3 ",c:" $4 ",o:" $5 ",f:" $6 ",g:" $7 ",s:" $8;
      next
    }
    {
      print
    }
  ' $f >${f%_label.fasta}_sintax.fasta
done
```

### Unite-style labels

Some algorithms are configured to use taxonomic annotations as used by the Unite
database. These are similar to Sintax, but use different delimiters.

```sh
for f in train_*_sintax.fasta
do
  sed -r 's/;tax=/|/; s/([kpcofgst]):/\1__/g; y/,/;/' ${f}\
      >${f%sintax.fasta}unite.fasta
done
```

## ID-taxonomy maps (without sequences)

Many of these formats use two different delimiters; one to delimit the ID from
the taxonomic classification, and one to delimit taxa at different ranks in the
classification.

### Basic (space+pipe)

```sh
for type in test train
do
  awk -F" " '
    /^>/ {
      sub(/^>/, "", $1)
      print
    }
  ' ${type}_nt_label.fasta >${type}_tax.txt
done
```

### TSV with headers

This is the format used by Dnabarcoder. It is also convenient to load into R
(although the ranks need to be renamed).

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
