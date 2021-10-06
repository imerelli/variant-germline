#!/bin/bash

# da fare partire in output/
echo Run replace RD

# Replace Read Groups:

picard AddOrReplaceReadGroups I=/home/smanessi/scRNAvcall/input/bam/$1\_possorted_genome_bam.bam O=/home/smanessi/scRNAvcall/output/$1\_out/$1\_rg.bam RGID=$1 RGLB=scrnaseq RGPL=ILLUMINA RGSM=PT01 RGPU=Nextseq CREATE_INDEX=true CREATE_MD5_FILE=true SO=coordinate

echo Run MarkDuplicates
# Marking duplicates:

picard MarkDuplicates I=$1\_out/$1\_rg.bam O=$1\_out/$1\_rg_md.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=$1\_out/$1\_rg_md.metrics

echo Run SplitNCigarReads
# SplitnCigar: split reads that span exons and remove N chars (junction)

gatk SplitNCigarReads -R /home/smanessi/scRNAvcall/metadata/genome.fa -I $1\_out/$1\_rg_md.bam -O $1\_out/$1\_rg_md_spl.bam

echo Run Base Quality Recalibration
# Base Quality Recalibration: detect and correct systematic errors in base quality scores.

gatk BaseRecalibrator -I $1\_out/$1\_rg_md_spl.bam -R /home/smanessi/scRNAvcall/metadata/genome.fa -O $1\_out/$1\_rg_md_spl_bqsr.table --known-sites /home/smanessi/scRNAvcall/metadata/dbSNP.v152.vcf.gz

gatk ApplyBQSR -I $1\_out/$1\_rg_md_spl.bam -O $1\_out/$1\_rg_md_spl_bqsr.bam --bqsr-recal-file $1\_out/$1\_rg_md_spl_bqsr.table --emit-original-quals true -R /home/smanessi/scRNAvcall/metadata/genome.fa

mv $1\_out/$1\_rg_md_spl_bqsr.bai $1\_out/$1\_rg_md_spl_bqsr.bam.bai #per far funzionare poi l'analisi vartrix 

# Evaluate and compare base quality score recalibration tables

gatk AnalyzeCovariates -bqsr $1\_out/$1\_rg_md_spl_bqsr.table -plots $1\_out/$1\_rg_md_spl_bqsr.pdf -csv $1\_out/$1\_rg_md_spl_bqsr.csv

echo Run Variant call 
# Variant call with HaplotypeCaller

gatk HaplotypeCaller -R /home/smanessi/scRNAvcall/metadata/genome.fa -I $1\_out/$1\_rg_md_spl_bqsr.bam -stand-call-conf 20.0 -O $1\_out/$1\_rg_md_spl_bqsr.vcf -bamout $1\_out/$1\_rg_md_spl_bqsr_hc.bam

# vcf backup
cp $1\_out/$1\_rg_md_spl_bqsr.vcf $1\_out/$1\_rg_md_spl_bqsr.vcf.copy

echo Run Variant Filtration
# Variant filtration: place tags according to Broad institute guidelines

gatk VariantFiltration -R /home/smanessi/scRNAvcall/metadata/genome.fa -V $1\_out/$1\_rg_md_spl_bqsr.vcf -window 35 -cluster 3 --filter-name FS --filter-expression "FS > 30.0" --filter-name QD --filter-expression "QD < 2.0" --filter-name LowDP --filter-expression "DP < 7"  --genotype-filter-name HQ --genotype-filter-expression "GQ > 90" -O $1\_out/$1\_rg_md_spl_bqsr_vft.vcf -OVI true -OVM true

echo Run Annotation with SnpEff
# Annotation using snpeffect
snpEff -v -csvStats $1\_out/$1\_rawstats_canon.csv -s $1\_out/$1\_rawstats_canon.html -canon GRCh38.86 $1\_out/$1\_rg_md_spl_bqsr_vft.vcf > $1\_out/$1\_raw_annot.vcf

# Adding info by using SnpSift:

echo Cosmic Annotatio Coding
# Cosmic annotation (coding)
SnpSift annotate /home/smanessi/scRNAvcall/metadata/CosmicCodingMuts.vcf.gz $1\_out/$1\_raw_annot.vcf > $1\_out/$1\_cosmicC_annot.vcf
gatk IndexFeatureFile -I $1\_out/$1\_cosmicC_annot.vcf

echo Cosmic Annotation non coding
# Cosmic annotation (non coding)

SnpSift annotate /home/smanessi/scRNAvcall/metadata/CosmicNonCodingVariants.vcf.gz $1\_out/$1\_cosmicC_annot.vcf > $1\_out/$1\_cosmic_annot.vcf
gatk IndexFeatureFile -I $1\_out/$1\_cosmic_annot.vcf

echo dbSNP annotatio
# dbSNP annotation

SnpSift annotate /home/smanessi/scRNAvcall/metadata/dbSNP.v152.vcf.gz $1\_out/$1\_cosmic_annot.vcf > $1\_out/$1\_cosmic_dbsnp_annot.vcf
vcftools --vcf $1\_out/$1\_cosmic_dbsnp_annot.vcf --remove-filtered-all --recode --recode-INFO-all --out $1\_out/$1\_cosmic_dbsnp_annot_PASS.vcf
gatk IndexFeatureFile -I $1\_out/$1\_cosmic_dbsnp_annot.vcf
gatk IndexFeatureFile -I $1\_out/$1\_cosmic_dbsnp_annot_PASS.vcf.recode.vcf

echo subset missens variants
# subsetting only missense variants

SnpSift filter "ANN[*].EFFECT has 'missense_variant'" $1\_out/$1\_cosmic_dbsnp_annot.vcf > $1\_out/$1\_missense.vcf
vcftools --vcf $1\_out/$1\_missense.vcf --remove-filtered-all --recode --recode-INFO-all --out $1\_out/$1\_missense_PASS.vcf
gatk IndexFeatureFile -I $1\_out/$1\_missense.vcf
gatk IndexFeatureFile -I $1\_out/$1\_missense_PASS.vcf.recode.vcf

echo subset inframe INDEL
# subsetting only inframe deletions/insertions

SnpSift filter "(ANN[*].EFFECT has 'disruptive_inframe_deletion') || (ANN[*].EFFECT has 'disruptive_inframe_insertion')" $1\_out/$1\_cosmic_dbsnp_annot.vcf > $1\_out/$1\_inframe_indels.vcf
vcftools --vcf $1\_out/$1\_inframe_indels.vcf --remove-filtered-all --recode --recode-INFO-all --out $1\_out/$1\_inframe_indels_PASS.vcf
gatk IndexFeatureFile -I $1\_out/$1\_inframe_indels.vcf
gatk IndexFeatureFile -I $1\_out/$1\_inframe_indels_PASS.vcf.recode.vcf

echo subset stop lost/gain
# subsetting stop lost/gain

SnpSift filter "(ANN[*].EFFECT has 'stop_lost') || (ANN[*].EFFECT has 'stop_gained') || (ANN[*].EFFECT has 'start_lost')" $1\_out/$1\_cosmic_dbsnp_annot.vcf > $1\_out/$1\_stop_lossgain.vcf
vcftools --vcf $1\_out/$1\_stop_lossgain.vcf --remove-filtered-all --recode --recode-INFO-all --out $1\_out/$1\_stop_lossgain_PASS.vcf
gatk IndexFeatureFile -I $1\_out/$1\_stop_lossgain.vcf
gatk IndexFeatureFile -I $1\_out/$1\_stop_lossgain_PASS.vcf.recode.vcf

echo subset splice donor/acceptor
# subsetting splice donor/acceptor

SnpSift filter "(ANN[*].EFFECT has 'splice_donor') || (ANN[*].EFFECT has 'splice_acceptor') " $1\_out/$1\_cosmic_dbsnp_annot.vcf > $1\_out/$1\_splice_donoracc.vcf
vcftools --vcf $1\_out/$1\_splice_donoracc.vcf --remove-filtered-all --recode --recode-INFO-all --out $1\_out/$1\_splice_donoracc_PASS.vcf
gatk IndexFeatureFile -I $1\_out/$1\_splice_donoracc.vcf
gatk IndexFeatureFile -I $1\_out/$1\_splice_donoracc_PASS.vcf.recode.vcf

echo subset frameshift
# subsetting frameshift

SnpSift filter "(ANN[*].EFFECT has 'frameshift_variant')" $1\_out/$1\_cosmic_dbsnp_annot.vcf > $1\_out/$1\_frameshift.vcf
vcftools --vcf $1\_out/$1\_frameshift.vcf --remove-filtered-all --recode --recode-INFO-all --out $1\_out/$1\_frameshift_PASS.vcf
gatk IndexFeatureFile -I $1\_out/$1\_frameshift.vcf
gatk IndexFeatureFile -I $1\_out/$1\_frameshift_PASS.vcf.recode.vcf


