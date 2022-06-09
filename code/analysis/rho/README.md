# Estimating population scaled recombination rate (ρ)


# Run LDhat with genotype data
# Producing input for ldhat run from bin
# In /Users/homapapoli/Documents/projects/ostrich_Z/ldhat_dir
# conda activate 
# conda activate ldhat
./vcf_to_ldhat_out.py black.PAR.hwe.filtered.vcf.gz par_scaf.bed 20 19 superscaffold36 100
export vcf=/proj/snic2020-16-269/private/homap/ostrich_z/data/vcf/black.PAR.PASS.biallelic.nomissing.hwe.allhet.fixedalt.fixedref.vcf.gz
./vcf_to_LDhat_genotype.py $vcf par_scaf.bed 20 19 superscaffold36 100
#*********************************************************************************************#
# In personal computer /Users/homapapoli/Documents/projects/ostrich_Z/4_Ostrich_polymorphism/manuscript_tables
# Create Z coordinates 
for species in black blue red
do 
echo $species
python scaffold_to_chr_vcf.py sfs_measures/${species}.PAR.sfs.txt sfs_measures/${species}.nonPAR.sfs.txt > sfs_measures/${species}.sfs.Z.txt
done

# Prepare LD output for plotting in R
# In /proj/snic2020-16-269/private/homap/ostrich_z/data/LD/scaf.split 
mkdir -p LD_chromosome_plot
# Run LD_plot_scaffold.sh in the /proj/snic2020-16-269/private/homap/ostrich_z/bin

#*********************************************************************************************#
# Get the PAR fasta sequence
reference/black.repeat.depth.masked.fa
# Convert fasta sequence into consesuns for each individual with SNPs
../data/vcf/${species}.PAR.filtered.vcf.gz


for segment in superscaffold36:3524263-9394175 superscaffold35:1-4625539 superscaffold54:1-16379243 superscaffold26:1-25310599
do 
echo $segment
samtools faidx reference/Struthio_camelus.20130116.OM.fa $segment | bcftools consensus consensus_fasta/black.PAR.filtered.vcf.gz -H 2 -s P1878_107 -o consensus_fasta/P1878_107.${segment}.fa
samtools faidx reference/Struthio_camelus.20130116.OM.fa $segment | bcftools consensus consensus_fasta/black.PAR.filtered.vcf.gz -H 2 -s P1878_108 -o consensus_fasta/P1878_108.${segment}.fa
samtools faidx reference/Struthio_camelus.20130116.OM.fa $segment | bcftools consensus consensus_fasta/black.PAR.filtered.vcf.gz -H 2 -s P1878_109 -o consensus_fasta/P1878_109.${segment}.fa
samtools faidx reference/Struthio_camelus.20130116.OM.fa $segment | bcftools consensus consensus_fasta/black.PAR.filtered.vcf.gz -H 2 -s P1878_110 -o consensus_fasta/P1878_110.${segment}.fa
samtools faidx reference/Struthio_camelus.20130116.OM.fa $segment | bcftools consensus consensus_fasta/black.PAR.filtered.vcf.gz -H 2 -s P1878_111 -o consensus_fasta/P1878_111.${segment}.fa
samtools faidx reference/Struthio_camelus.20130116.OM.fa $segment | bcftools consensus consensus_fasta/black.PAR.filtered.vcf.gz -H 2 -s P1878_112 -o consensus_fasta/P1878_112.${segment}.fa
samtools faidx reference/Struthio_camelus.20130116.OM.fa $segment | bcftools consensus consensus_fasta/black.PAR.filtered.vcf.gz -H 2 -s P1878_113 -o consensus_fasta/P1878_113.${segment}.fa
samtools faidx reference/Struthio_camelus.20130116.OM.fa $segment | bcftools consensus consensus_fasta/black.PAR.filtered.vcf.gz -H 2 -s P1878_114 -o consensus_fasta/P1878_114.${segment}.fa
samtools faidx reference/Struthio_camelus.20130116.OM.fa $segment | bcftools consensus consensus_fasta/black.PAR.filtered.vcf.gz -H 2 -s P1878_115 -o consensus_fasta/P1878_115.${segment}.fa
samtools faidx reference/Struthio_camelus.20130116.OM.fa $segment | bcftools consensus consensus_fasta/black.PAR.filtered.vcf.gz -H 2 -s P1878_116 -o consensus_fasta/P1878_116.${segment}.fa
done

# Produce 2000 SNP intervals for each scaffold
# Running in /proj/snic2020-16-269/private/homap/ostrich_z/result/ldhat/ldhat_dir/
python vcf_to_ldhat_out_black.py vcf_rho/black.PAR.hwe.filtered.vcf.gz vcf_rho/par_scaf.bed 2000 500 superscaffold36 1000
python vcf_to_ldhat_out_black.py vcf_rho/black.PAR.hwe.filtered.vcf.gz vcf_rho/par_scaf.bed 2000 500 superscaffold35 1000
python vcf_to_ldhat_out_black.py vcf_rho/black.PAR.hwe.filtered.vcf.gz vcf_rho/par_scaf.bed 2000 500 superscaffold54 1000
python vcf_to_ldhat_out_black.py vcf_rho/black.PAR.hwe.filtered.vcf.gz vcf_rho/par_scaf.bed 2000 500 superscaffold26 1000

# Create the exact likelihood table
./interval -seq superscaffold36.2000.500.1.sites.txt -loc superscaffold36.2000.500.1.locs.txt -lk ostrich_genotypenew_lk.txt -prefix superscaffold36.first.run -its 10000000 -bpen 5 -samp 2000
./interval -seq superscaffold35.2000.500.1.sites.txt -loc superscaffold35.2000.500.1.locs.txt -lk ostrich_genotypenew_lk.txt -prefix superscaffold35.first.run -its 10000000 -bpen 5 -samp 2000
./interval -seq superscaffold54.2000.500.1.sites.txt -loc superscaffold54.2000.500.1.locs.txt -lk ostrich_genotypenew_lk.txt -prefix superscaffold54.first.run -its 10000000 -bpen 5 -samp 2000
./interval -seq superscaffold26.2000.500.1.sites.txt -loc superscaffold26.2000.500.1.locs.txt -lk ostrich_genotypenew_lk.txt -prefix superscaffold26.first.run -its 10000000 -bpen 5 -samp 2000

for i in {1..24} 
do
echo $i
sbatch ldhat_slurm.sh superscaffold36.2000.500.${i}.sites.txt superscaffold36.2000.500.${i}.locs.txt superscaffold36.first.runnew_lk.txt superscaffold36.2000.500.${i}.
done
./interval -seq test/superscaffold36.2000.500.25.sites.txt -loc test/superscaffold36.2000.500.25.locs.txt -lk ostrich_genotypenew_lk.txt -prefix test/superscaffold36.2000.500.25. -its 10000000 -bpen 5 -samp 2000

for i in {1..18} 
do
echo $i
sbatch ldhat_slurm.sh superscaffold35.2000.500.${i}.sites.txt superscaffold35.2000.500.${i}.locs.txt ostrich_genotypenew_lk.txt superscaffold35.2000.500.${i}.
done
./interval -seq test/superscaffold35.2000.500.19.sites.txt -loc test/superscaffold35.2000.500.19.locs.txt -lk ostrich_genotypenew_lk.txt -prefix test/superscaffold35.2000.500.19. -its 10000000 -bpen 5 -samp 2000

for i in {7..58} 
do
echo $i
sbatch ldhat_slurm.sh superscaffold54.2000.500.${i}.sites.txt superscaffold54.2000.500.${i}.locs.txt ostrich_genotypenew_lk.txt superscaffold54.2000.500.${i}.
done
./interval -seq test/superscaffold54.2000.500.59.sites.txt -loc test/superscaffold54.2000.500.59.locs.txt -lk ostrich_genotypenew_lk.txt -prefix test/superscaffold54.2000.500.59. -its 10000000 -bpen 5 -samp 2000

for i in {1..100} 
do
echo $i
sbatch ldhat_slurm.sh superscaffold26.2000.500.${i}.sites.txt superscaffold26.2000.500.${i}.locs.txt ostrich_genotypenew_lk.txt superscaffold26.2000.500.${i}.
done
./interval -seq test/superscaffold26.2000.500.101.sites.txt -loc test/superscaffold26.2000.500.101.locs.txt -lk ostrich_genotypenew_lk.txt -prefix test/superscaffold26.2000.500.101. -its 10000000 -bpen 5 -samp 2000

for i in {1..24} 
do
echo $i
./interval -seq superscaffold36.2000.500.${i}.sites.txt -loc superscaffold36.2000.500.${i}.locs.txt -lk ostrich_genotypenew_lk.txt -prefix superscaffold36.2000.500.${i}. -its 10000000 -bpen 5 -samp 2000
done

for i in {1..18} 
do
echo $i
./interval -seq superscaffold35.2000.500.${i}.sites.txt -loc superscaffold35.2000.500.${i}.locs.txt -lk ostrich_genotypenew_lk.txt -prefix superscaffold35.2000.500.${i}. -its 10000000 -bpen 5 -samp 2000
done

for i in {7..58} 
do
echo $i
./interval -seq superscaffold54.2000.500.${i}.sites.txt -loc superscaffold54.2000.500.${i}.locs.txt -lk ostrich_genotypenew_lk.txt -prefix superscaffold54.2000.500.${i}. -its 10000000 -bpen 5 -samp 2000
done

for i in {1..100} 
do
echo $i
./interval -seq superscaffold26.2000.500.${i}.sites.txt -loc superscaffold26.2000.500.${i}.locs.txt -lk ostrich_genotypenew_lk.txt -prefix superscaffold26.2000.500.${i}. -its 10000000 -bpen 5 -samp 2000
done

python vcf_to_ldhat_out_black.py vcf_rho/black.PAR.hwe.filtered.vcf.gz vcf_rho/par_scaf.bed 20 5 superscaffold36 10


for i in {1..10} 
do
echo $i
./interval -seq superscaffold36.20.5.${i}.sites.txt -loc superscaffold36.20.5.${i}.locs.txt -lk ostrich_genotypenew_lk.txt -prefix superscaffold36.20.5.${i}. -its 10000000 -bpen 5 -samp 2000
done
