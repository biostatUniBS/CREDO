# CREDO
CREDO: Highly confident disease-relevant A-to-I RNA-editing discovery in breast cancer

## Prerequisites

Java 1.6 higher
R 3.0 higher
R library:VariantAnnotation, data.table, gplots, VennDiagram, parallel, bit64
Samtools 1.4.1
Picard 1.12 higher
VarScan 2.3.7
GATKAnalysisTK 3.3.0

## Pipeline Running

**CREDO** requires sequencing data of RNA and Exome of matched sample (normal and disease sample) for a case. Processes from 1 to 3 should be carried out for each sample, i.e., normal or disease sample, independently in matched sample case. Then, confident RNA editing will be discovered based on comparison between matched samples in process 5. However, the pipeline for single sample case is designed too. In single sample case, user executes processes 1 to 3 and performs process 4 for a sample.

## RNA-seq processing to obtain read counts of each loci
Aligned bam files from different aligners can be obtained from any database or user’s own aligning jobs. User might want to use single bam file from one aligner, then you can skip 1.1. 

### 1 Merge bam files from multiple alignments

#### 1.1 Merge bam files
samtools merge -f {output_file} {input_bamfile_aligner1} {input_bamfile_aligner2}

#### 1.2 Assign a read group ID
java –jar picard.jar AddOrReplaceReadGroups I={input_file} O={output_file} SO=coordinate VALIDATION_STRINGENCY=SILENT RGID={read_group_id} RGLB=RNA RGPL=Illumina RGPU=illumina  RGSM={read_group_id}

#### 1.3 Remove duplicate reads 
samtools rmdup –S {input_file} {output_file}

#### 1.4 Reorder a bam file
java -jar picard.jar ReorderSam VALIDATION_STRINGENCY=SILENT  I={input_bamfile} O={output_bamfile} REFERENCE={reference_genome_fastafile}

#### 1.5 index a bam file
samtools index {bam_file}

#### 1.6 Pile up reads onto a genome
samtools mpileup -f {reference_genome_fasta} {input_file} | awk '{if($4 != 0) print $0}'  > {output_prefix}.mpileup

#### 1.7 Count the number of reads piled up on each loci
java -jar VarScan.v2.3.7.jar readcounts {input_mpileup_file} --min-base-qual {minimum_base_quality}  --output-file  {output_vcf_file}

#### 1.8 A-to-G loci Acquisition
java -jar FilterRNAAtoG.jar {input_readcnt_vcffile}  {output_filename}

### 2 DNA variant calling (aligned with GRCh37-lite.fa)

#### 2.1 Mark duplicate reads
Java -jar picard.jar MarkDuplicates I={input_file} O={output_file} REMOVE_DUPLICATES=true CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M={ output.metrics_filename }

#### 2.2 InDel realignment
java -jar pathgatk/GenomeAnalysisTK.jar -T RealignerTargetCreator -R litegenome_dir/GRCh37-lite.fa -I {input_file} -o {output_prefix}.list -known {known_indel_file} 
java -jar pathgatk/GenomeAnalysisTK.jar -T IndelRealigner -R litegenome_dir/GRCh37-lite.fa -I {input_file} -targetIntervals {output_prefix}.list -o {output_filename} -known {known_indel_file}
java -jar picard.jar FixMateInformation INPUT={input_file} OUTPUT={output_file_name} SO=coordinate VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true

#### 2.3 Base quality recalibrator
java -jar pathgatk/GenomeAnalysisTK.jar -T BaseRecalibrator -R litegenome_dir/GRCh37-lite.fa -I {input_file} -o {recal.table_name} -knownSites {known_dbsnp_indel_files} 
java -jar pathgatk/GenomeAnalysisTK.jar -T PrintReads -R litegenome_dir/GRCh37-lite.fa -I {input_file} -o {output_file} -BQSR { recal.table_name }

#### 2.4 Variant call
java -jar pathgatk/GenomeAnalysisTK.jar -T HaplotypeCaller -R litegenome_dir/GRCh37-lite.fa -I {input_file} -o {output_file}  -stand_call_conf 30.0 -stand_emit_conf 10.0

#### 2.5 Variant recalibration
java -jar pathgatk/GenomeAnalysisTK.jar -T VariantRecalibrator -R litegenome_dir/GRCh37-lite.fa -input {input_file} -resource:hapmap,known=false,training=true,truth=true,prior=15.0 genomedir/hapmap_3.3.b37.vcf  -resource:omni,known=false,training=true,truth=true,prior=12.0 genomedir/1000G_omni2.5.b37.vcf -resource:1000G,known=false,training=true,truth=false,prior=10.0 genomedir/1000G_phase1.snps.high_confidence.b37.vcf  -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 genomedir/dbsnp_138.b37.vcf  -an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile recalibrate_SNP.recal -tranchesFile recalibrate_SNP.tranches -rscriptFile recalibrate_SNP_plots.R                
java -jar pathgatk/GenomeAnalysisTK.jar -T ApplyRecalibration -R  litegenome_dir/GRCh37-lite.fa -input {input_file} -mode SNP --ts_filter_level 99.0 -recalFile recalibrate_SNP.recal -tranchesFile recalibrate_SNP.tranches -o {output_file}                
java -jar pathgatk/GenomeAnalysisTK.jar -T VariantRecalibrator -R litegenome_dir/GRCh37-lite.fa -input {input_file} -resource:mills,known=true,training=true,truth=true,prior=12.0 genome_dir/Mills_and_1000G_gold_standard.indels.b37.vcf -an QD -an DP -an FS -an SOR -an MQRankSum -an ReadPosRankSum  -mode INDEL -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 --maxGaussians 4 -recalFile recalibrate_INDEL.recal -tranchesFile recalibrate_INDEL.tranches -rscriptFile recalibrate_INDEL_plots.R                
java -jar pathgatk/GenomeAnalysisTK.jar -T ApplyRecalibration -R litegenome_dir/GRCh37-lite.fa -input {input_file} -mode INDEL --ts_filter_level 99.0 -recalFile recalibrate_INDEL.recal -tranchesFile recalibrate_INDEL.tranches -o {output_file}

#### 2.6 DNA A-to-G variant Acquisition
	bcftools filter -i 'AC==1 && REF=="A" && ALT=="G" && DP>4 && QUAL>=50' {input_vcf_file} -o {output_vcf_file}

### 3 DNA-seq processing to obtain read counts of each loci

#### 3.1 Pile up onto a genome
samtools mpileup -f {reference_genome_fasta} {output_bamfile_from2.3} | awk '{if($4 != 0) print $0}'  > {prefix}.mpileup

#### 3.2 Count the number of reads piled up on each loci
java -jar VarScan.v2.3.7.jar readcounts {input_mpileup_file} --min-base-qual {minimum_base-quality}  --output-file  {output_vcf_file}

#### 3.3 A-to-G loci Acquisition
java -jar FilterDNAAtoG.jar {input_readcnt_vcffile}  {output_filename}

### 4 RNA Editing sites discovery Process for single sample case

#### 4.1 Compute confidence zero score for DNA A-to-G loci
Rscript rnaediting_merged_computeaz_singlesample.R {output_file_from1.8} {output_file_from3.3} {output_filename}

#### 4.2 Discover confident RNA editing loci
Rscript rnaediting_merged_singlesample.R {output_file_from1.8} {output_file_from4.1} {output_file_from2.6} {output_filename}

### 5 RNA Editing sites discovery Process for matched sample case

#### 5.1 Compute confidence zero score for DNA A-to-G loci for both samples
Rscript rnaediting_merged_computeaz_matchedsample.R {sample1_output_file_from1.8} {sample1_output_file_from3.3} {sample2_output_file_from1.8} {sample2_output_file_from3.3} {sample1_output_filename} {sample2_output_filename}

#### 5.2 Discover confident RNA editing loci from matched samples
Rscript rnaediting_merged_matchedsample.R {sample1_output_file_from1.8} {sample1_output_file_from5.1} {sample1_output_file_from2.6} {sample2_output_file_from1.8} {sample2_output_file_from5.1} {sample2_output_file_from2.6}   {output_filename}
