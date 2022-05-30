# Matching ostrich scaffolds to chicken chromosomes

Using LASTZ, we identify those scaffolds matching to chromosomes 4 and 5 of chicken. 

`bash lastz_run.sh`

The output has the following structure. The top two rows are shown where NC_006090.5 is a
chicken chromosome and scaffold53 is the ostrich scaffold.

Table 1. Top two rows of the LASTZ run of ostrich assembly against chicken assembly.
| #name1  | start1 | end1  |  strand1| name2  | start2+ |end2+  | strand2 |score|
| ------- | ------ | ----- |  ------ | -----  | ------- | ----- | ------- |-----|
| NC_006090.5   |   797670  |849736  |+     |  scaffold53    |  2196657 |2250872| +    |   2082544|
| NC_006090.5   |  900362 | 949856  |+     |  scaffold53   |   2311895| 2366219 |+     | 1902427|


## Obtain ostrich scaffolds matching with chicken chromosomes 4 and 5 

```
python chicken_to_ostrich_autosomes.py ../../../data/lastz/chicken_NC_chr.txt ../../../data/lastz/chicken_ostrich.lastz \
../../../data/genome/Struthio_camelus.20130116.OM.fa.fai > ../../../data/lastz/gg_ostrich_macrochr.txt
```

## Create bed files with ostrich scaffolds matching chromosomes 4 and 5 of chicken

```
for chrom in {4..5}; do
echo chromosome_${chrom}
grep  chromosome_${chrom} ../../../data/lastz/gg_ostrich_macrochr.txt | \
awk 'BEGIN{print "chrom""\t""chromStart""\t""chromEnd"}{print $4"\t""0""\t"$5}' | uniq > ../../../data/lastz/gg_chr${chrom}_ostrich.bed
awk '{if(NR > 1) sum+=$3} END{print sum}' ../../../data/lastz/gg_chr${chrom}_ostrich.bed; done
```

Table 2. Matches between ostrich scaffolds and chromosome 4 of chicken
|chrom  | chromStart    |  chromEnd|
| ----- | ------------- | -------- |
|scaffold198   |  0     |  3123540|
|scaffold554   |  0    |   1015449|
|scaffold25    |  0    |   2288681|
|superscaffold11 |0     |  41821358|
|superscaffold31 |0     |  4751311|
|superscaffold28 |0     |  9720582|
|scaffold120     |0     |  1010486|
|superscaffold57 |0     |  17714005|
|superscaffold44 |0     |  9842373|
|superscaffold66 |0    |   6810023|
|scaffold1404   | 0    |   184221|
|scaffold1321    |0   |    133612|

Table 2. Matches between ostrich scaffolds and chromosome 5 of chicken
|chrom   chromStart      chromEnd
| ----- | ------------- | -------- |
|scaffold928    | 0     |  280674|
|superscaffold8 | 0     |  36774883|
|superscaffold41| 0     |  5004022|
|superscaffold86| 0     |  776085|
|scaffold760    | 0     |  412945|
|scaffold1008   | 0     |  44810|
|scaffold1128   | 0     |  102074|
|scaffold967    | 0     |  155924|

