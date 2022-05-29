
## Move the final autosomal vcf into the directory called LD_vcf
```
mkdir -p ../data/vcf/LD_vcf
mv ../data/vcf/${species}.*.hwe.filtered.vcf.gz ../data/vcf/LD_vcf
```

## Concatenate PAR and non-PAR positions
```
zgrep -v "#" ../data/vcf/LD_vcf/${species}.A.hwe.filtered.vcf.gz | \
awk '{print $1"\t"$2}' > ../data/vcf/LD_vcf/${species}.A.pos

zgrep -v "#" ../data/vcf/LD_vcf/${species}.PAR.hwe.filtered.vcf.gz | \
awk '{print $1"\t"$2}' > ../data/vcf/LD_vcf/${species}.PAR.pos

zgrep -v "#" ../data/vcf/LD_vcf/${species}.nonPAR.hwe.filtered.vcf.gz | \
awk '{print $1"\t"$2}' > ../data/vcf/LD_vcf/${species}.nonPAR.pos

cat ../data/vcf/LD_vcf/${species}.PAR.pos ../data/vcf/LD_vcf/${species}.nonPAR.pos | \
sort -k1,1 -k2n > ../data/vcf/LD_vcf/${species}.Z.pos

awk '{print $1"\t"$2-1"\t"$2}' ../data/vcf/LD_vcf/${species}.Z.pos > ../data/vcf/LD_vcf/${species}.Z.pos.bed 
python scaffold_to_chr_vcf.py ../data/vcf/LD_vcf/${species}.Z.pos.bed > ../data/vcf/LD_vcf/${species}.Z.pos.coordinates.txt
```