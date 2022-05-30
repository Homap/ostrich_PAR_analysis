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
| chrom  | chromlength |
| ----- | -------- |
|scaffold198   |  3123540|
|scaffold554   |    1015449|
|scaffold25    |   2288681|
|superscaffold11 |  41821358|
|superscaffold31 | 4751311|
|superscaffold28 | 9720582|
|scaffold120     |  1010486|
|superscaffold57 | 17714005|
|superscaffold44 |  9842373|
|superscaffold66 |  6810023|
|scaffold1404   |   184221|
|scaffold1321    |  133612|

Table 3. Matches between ostrich scaffolds and chromosome 5 of chicken
|chrom  | chromEnd
| ----- | -------- |
|scaffold928    | 280674|
|superscaffold8 | 36774883|
|superscaffold41| 5004022|
|superscaffold86| 776085|
|scaffold760    | 412945|
|scaffold1008   | 44810|
|scaffold1128   | 102074|
|scaffold967    | 155924|

