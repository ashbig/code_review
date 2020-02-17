#!/bin/bash
# Ashkan Bigdeli
#
# 5/1/2018

################## Data Bases ########
## Databases
hg38='/project/cpdlab/dev/ash/capture/cap_databases/v0/Homo_sapiens_assembly38.fasta'

db_dict='/project/cpdlab/dev/ash/capture/cap_databases/v0/Homo_sapiens_assembly38.dict'
mills_gold='/project/cpdlab/dev/ash/capture/cap_databases/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'
dbsnp_38='/project/cpdlab/dev/ash/capture/cap_databases/broad_ftp_bundle/dbsnp_138.hg38.vcf.gz'
gnomad_af='/project/cpdlab/dev/ash/capture/cap_databases/af-only-gnomad.hg38.vcf.gz'
gnomad_annotation='/project/cpdlab/dev/ash/capture/cap_databases/af-only-gnomad.hg38.annotation.vcf.gz'
cosmic_annotation='/project/cpdlab/dev/ash/capture/cap_databases/CosmicCodingMuts.with_chr.vcf.gz'
capture_PON='/project/cpdlab/dev/ash/capture/files/GBP_pon.vcf.gz'
MCL='/project/cpdlab/dev/ash/capture/files/MCL_v1.0.bed'
################################### Programs ###############

stats='/project/cpdlab/dev/ash/capture/methods/stat_caller.py'
simpletation='/project/cpdlab/dev/ash/capture/methods/simpletation.py'

# temporary dev settings

probe_bed='/project/cpdlab/dev/ash/capture/files/gbp.probes.hg38.bed'
probe_intervals='/project/cpdlab/dev/ash/capture/files/gbp.probes.hg38.intervals'
target_bed='/project/cpdlab/dev/ash/capture/files/gbp_rd.coverage.hg38.bed'
target_intervals='/project/cpdlab/dev/ash/capture/files/gbp_rd.coverage.hg38.intervals'
cpd_bed='/project/cpdlab/dev/ash/capture/files/CPD_mustcall.bed'
backbone_snp_bed='/project/cpdlab/dev/ash/capture/files/xgen-human-id-research-panel-probes.hg38.bed'
backbone_cnv_bed='/project/cpdlab/dev/ash/capture/files/xgen_cnv_backbone_panel_probes.hg38.bed'
backbone_cnv_target_bed='/project/cpdlab/dev/ash/capture/files/xgen_cnv_backbone_panel_target_snps.hg38.bed'

#load environment
source activate cap_dev

# RUN AS:
# sh /project/cpdlab/dev/ash/capture/capture_dev.sh
#       /project/cpdlab/dev/ash/capture/{runfolder}
#       /project/cpdlab/dev/ash/capture/{runfolder}/{sample}_R1_001.fastq.gz
#       /project/cpdlab/dev/ash/capture/{runfolder}/{sample}_R2_001.fastq.gz
#       {sample}

out_dir=$1
read1=$2
read2=$3
name=$4
cpu=$5
mem=$6

out_dir=${out_dir}/${name}
tmp=${out_dir}/tmp
mkdir -p ${out_dir}
cd ${out_dir}

#get fastq naming properties
FW=$(basename $read1 ".fastq.gz")
BW=$(basename $read2 ".fastq.gz")

##################### Trim, Align, BAM  #####################

trim_fq1=${FW}_val_1.fq.gz
trim_fq2=${BW}_val_2.fq.gz
out_sam=${name}.align.sam
align_stat=${name}.align.stat
out_bam=${name}.align.bam

trim_galore --illumina --fastqc --paired ${read1} ${read2} -o ${out_dir} --path_to_cutadapt cutadapt

bwa mem -t ${cpu} -K 100000000 ${hg38} ${trim_fq1} ${trim_fq2}  > ${out_sam}

samtools flagstat ${out_sam} > ${align_stat}

samtools view -Su ${out_sam} | samtools sort -m 30000000000 -o ${out_bam}

#################### Deduplicate, Fix, Sort, Clean, Validate ####################

dedup_out=${name}.nodup.bam
dedup_fix_out=${name}.nodup.fix.bam
dedup_fix_sort=${name}.nodup.fix.srt.bam
dup_qc=${name}.dup.qc
clean=${name}.clean.bam
clean_vld=${name}.clean.vr
dedup_stat=${name}.dedup.stat

picard -Xmx${mem}m MarkDuplicates TMP_DIR=${tmp} M=${dup_qc} I=${out_bam} O=${dedup_out}  REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT

picard -Xmx${mem}m AddOrReplaceReadGroups TMP_DIR=${tmp} I=${dedup_out} O=${dedup_fix_out} RGID=${name} RGLB="gbp_test" RGPL=illumina RGPU=unit1 RGSM=${name}

picard -Xmx${mem}m SortSam TMP_DIR=${tmp} I=${dedup_fix_out} O=${dedup_fix_sort} SORT_ORDER=coordinate CREATE_INDEX=true

picard -Xmx${mem}m CleanSam TMP_DIR=${tmp} I=${dedup_fix_sort} O=${clean} CREATE_INDEX=true

picard -Xmx${mem}m ValidateSamFile I=${clean} O=${clean_vld}

samtools flagstat ${clean} > ${dedup_stat}

################################## Base Recallibration #################################

recal_table=${name}.recal_data.table
recal_bam=${name}.final.bam
recal_vld=${name}.recal.vr
after_recal_table=${name}.after_recal.table
recal_plot=${name}.recal.pdf
recal_stat=${name}.recal.stat
gatk --java-options -Xmx${mem}m BaseRecalibrator -R ${hg38} --known-sites ${mills_gold} --known-sites ${dbsnp_38} -I ${clean} -O ${recal_table}
gatk --java-options -Xmx${mem}m ApplyBQSR -R ${hg38} -I ${clean} --bqsr-recal-file ${recal_table} -O ${recal_bam}
gatk --java-options -Xmx${mem}m BaseRecalibrator -R ${hg38} --known-sites ${mills_gold} --known-sites ${dbsnp_38} -I ${clean} -O ${after_recal_table}
gatk --java-options -Xmx${mem}m AnalyzeCovariates -before ${recal_table} -after ${after_recal_table} -plots ${recal_plot}
picard -Xmx${mem}m ValidateSamFile I=${recal_bam} O=${recal_vld}
samtools flagstat ${recal_bam} > ${recal_stat}

#################################### Calculate Sample Statistics ####################################

#get metrics
hs_metrics=${name}.hs_metrics.tsv
target_coverage=${name}.target_coverage.tsv
picard -Xmx${mem}m CollectHsMetrics I=${recal_bam} O=${hs_metrics} R=${hg38} BAIT_INTERVALS=${probe_intervals} TARGET_INTERVALS=${target_intervals} PER_TARGET_COVERAGE=${target_coverage} # this should be added when it permits MQ=10 Q=10 COVMAX=null

#get gene based mean coverages to balance probe set
#add naming to target coverage files
gene_mean_cvg=${name}.mean_gene_coverage.tsv
echo -e "NAME\tGENE\tMEAN_COVERAGE" > ${gene_mean_cvg}
awk -v var=${name} '{seen[$5]+=$7; count[$5]++} END{for (x in seen)print var,x, seen[x]/count[x]}' OFS='\t' ${target_coverage} >> ${gene_mean_cvg}

named_tc=${name}.named.target_coverage.tsv
awk -v var=${name} '{print var, $0}' OFS='\t' ${target_coverage} > ${named_tc}

#get interval metrics
insert_size_metrics_txt=${name}.insert_size_metrics.txt
insert_size_metrics_pdf=${name}.insert_size_metrics.pdf
picard -Xmx${mem}m CollectInsertSizeMetrics I=${recal_bam} O=${insert_size_metrics_txt} H=${insert_size_metrics_pdf}

sample_stats=${name}.stats.tsv
dropped_exons=${name}.dropped_exons.tsv
python ${stats} -n ${name} -a ${align_stat} -l ${target_coverage} -hs ${hs_metrics} -dqc ${dup_qc} -im ${insert_size_metrics_txt} -de ${dropped_exons} -o ${sample_stats}


probe_bedcov=${name}.probe_bedcov.no_name.tsv
probe_bedcov_named=${name}.probe_bedcov.tsv
samtools bedcov ${probe_bed} ${recal_bam} > ${probe_bedcov}
awk -v var=${name} '{print var, $0, $6/($3-$2)}' OFS='\t' ${probe_bedcov} > ${probe_bedcov_named}


target_bedcov=${name}.target_bedcov.no_name.tsv
target_bedcov_named=${name}.target_bedcov.tsv
samtools bedcov ${target_bed} ${recal_bam} > ${target_bedcov}
awk -v var=${name} '{print var, $0, $6/($3-$2)}' OFS='\t' ${target_bedcov} > ${target_bedcov_named}

#gene_mean_cvg=${name}.mean_gene_coverage.tsv
#echo -e "NAME\tGENE\tMEAN_COVERAGE" > ${gene_mean_cvg}
#awk -v var=${name} '{seen[$5]+=$7; count[$5]++} END{for (x in seen)print var,x, seen[x]/count[x]}' OFS='\t' ${target_coverage} >> ${gene_mean_cvg}

#Mustcalllist#
cpd_bedcov=${name}.mustcall.target_bedcov.no_name.tsv
cpd_bed_named=${name}.cpdmustcall.named.bed
#samtools bedcov ${cpd_bed} ${recal_bam} > ${cpd_bedcov}
#awk -v var=${name} '{print var, $0, $5/($3-$2)}' OFS='\t' ${cpd_bedcov} > ${cpd_bed_named}

#Backbone#
backbone_snp=${name}.backbone.snp.noname.tsv
backbone_cnv=${name}.backbone.cnv.noname.tsv
backbone_cnv_target=${name}.backbone.cnv.tsnp.noname.tsv
backbone_snp_named=${name}.backbone.snp.named.tsv
backbone_cnv_named=${name}.backbone.cnv.named.tsv
backbone_cnv_target_named=${name}.backbone.cnv.tsnp.named.tsv

#samtools bedcov ${backbone_snp_bed} ${recal_bam} > ${backbone_snp}
#awk -v var=${name} '{print var, $0, $7/($3-$2)}' OFS='\t' ${backbone_snp} > ${backbone_snp_named}

#samtools bedcov ${backbone_cnv_bed} ${recal_bam} > ${backbone_cnv}
#awk -v var=${name} '{print var, $0, $5/($3-$2)}' OFS='\t' ${backbone_cnv} > ${backbone_cnv_named}

#samtools bedcov ${backbone_cnv_target_bed} ${recal_bam} > ${backbone_cnv_target}
#awk -v var=${name} '{print var, $0, $7/($3-$2)}' OFS='\t' ${backbone_cnv_target} > ${backbone_cnv_target_named}
#################################### Variant Calls ####################################

mutect_vcf=${name}.mutect.vcf.gz
somatic_bamout=${name}.mutect2.bamout.bam
gatk --java-options -Xmx${mem}m Mutect2 -R ${hg38} -I ${recal_bam} -tumor ${name} --germline-resource ${gnomad_af} -pon ${capture_PON} -O ${mutect_vcf} -bamout ${somatic_bamout}

# Run Haplotype Caller for Germline Variants
haplotypecaller_vcf=${name}.haplotypecaller.vcf.gz
germline_bamout=${name}.haplotypecaller.bamout.bam
gatk --java-options -Xmx${mem}m HaplotypeCaller -R ${hg38} -I ${recal_bam} -O ${haplotypecaller_vcf} -bamout ${germline_bamout}

#Run Haplotypecaller for gvcf
#haplo_gvcf=${name}.g.vcf
#sample_activity=${name}.activityprofile.IGV
#sample_assembly=${name}.assembly.IGV
#Tool=HaplotypeCaller
#error_model=GVCF
#erc_vcf=${name}.erc.vcf

#gatk --java-options -Xmx${mem}m ${Tool} -R ${hg38} -I ${recal_bam} -ERC ${error_model} --dbsnp ${dbsnp_38} -L ${target_intervals} -O ${haplo_gvcf} --activity-profile-out ${sample_activity} --assembly-region-out ${sample_assembly} --geno
typing-mode GENOTYPE_GIVEN_ALLELES --output-mode EMIT_ALL_SITES
#gatk --java-options -Xmx${mem}m GenotypeGVCFs  -R ${hg38}  --variant ${haplo_gvcf} -O ${erc_vcf}
##pileup=${name}.gatkpile.txt
###gatk --java-options -Xmx${mem}m Pileup -R ${hg38} -I ${recal_bam} -L ${probe_intervals} -O ${pileup}

###calcualte bases under 500##
#Persitebed=/project/cpdlab/dev/ash/capture/files/PerBaseCaller.bed
#perbaseout=${name}.genomecov
#PerbaseReport=${name}.PerBaseReport.txt
#Final_base_report=${name}.BaseReport
#bedtools genomecov -bga  -i ${Persitebed} -ibam ${recal_bam} > ${perbaseout}
#bedtools intersect -wb -a ${Persitebed} -b ${perbaseout} | awk 'OFS="\t" {print $4,$5,$6,$7}'  > ${PerBaseReport}
#awk ' $4<=500 {print}' ${PerbaseReport} | wc  -l > BasesUnder500.txt
#awk '{print}' ${PerbaseReport} | wc -l > TotalBases.txt
#cat BasesUnder500.txt TotalBases.txt > CombinedBases.txt
#awk ' NR == 1 { a = $1 }
#              NR == 2 { b = $1 }
#                     END {  print a/b } ' CombinedBases.txt > BaseReport.txt
#rm BasesUnder500.txt TotalBases.txt CombinedBases.txt

MC_out=${name}.bcf
MC_VCF=${name}.mc.vcf
###MUST CALL###
bcftools mpileup --threads 8 -Ou -R ${MCL} -f ${hg38} ${recal_bam}| bcftools call -cMA -Ob  -o ${MC_out}
bcftools view ${MC_out} > ${MC_VCF}

#################################### Variant Annotation ####################################

gs_merged_vcf=${name}.gatk_merged.vcf.gz
bcftools concat -a -d none -O z ${mutect_vcf} ${haplotypecaller_vcf} -o ${gs_merged_vcf}
tabix -p vcf ${gs_merged_vcf}

add_cosmic=${name}.gatk_merged.cosmic.vcf.gz
bcftools annotate -O z -a ${cosmic_annotation} -c ID,INFO ${gs_merged_vcf} -o ${add_cosmic}
tabix -p vcf ${add_cosmic}

add_gnomad=${name}.gatk_merged.cosmic.gnomad.vcf.gz
bcftools annotate -O z -a ${gnomad_annotation} -c ID,INFO ${add_cosmic} -o ${add_gnomad}
tabix -p vcf ${add_gnomad}


add_snpeff=${name}.anno.vcf.gz
snpEff -Xmx${mem}m -v -canon hg38 ${add_gnomad} | bgzip > ${add_snpeff}
tabix -p vcf ${add_snpeff}

final_vcf=${name}.final.vcf.gz
bcftools norm -m -any -O z -f ${hg38} ${add_snpeff} -o ${final_vcf}
tabix -p vcf ${final_vcf}

vcf2tsv=${name}.selected_fields.tsv
gatk --java-options -Xmx${mem}m VariantsToTable -R ${hg38} -V ${final_vcf} -O ${vcf2tsv} -L ${target_bed}  -F CHROM -F POS -F ID -F REF -F ALT -F QUAL \
-F AF -F DP -F MQ -F ANN -F CNT -F GNOMAD_AC -F GNOMAD_AF -GF AD -GF AF -GF MMQ -GF MBQ

# Annotations

perl_exec='/project/cpdlab/anaconda2/bin/perl'

snpEff='/project/cpdlab/anaconda2/envs/cap_dev/bin/snpEff'
vep='/project/cpdlab/Tools/ensembl-vep'
ref='/project/cpdlab/dev/ash/capture/cap_databases/v0/Homo_sapiens_assembly38.fasta'
vcf2maf='/project/cpdlab/Tools/.vep/mskcc-vcf2maf-5453f80'
vep_cache='/project/cpdlab/Tools/.vep'

zcat ${name}.mutect.vcf.gz > ${name}.mutect.vcf
zcat ${name}.haplotypecaller.vcf.gz > ${name}.haplotypecaller.vcf
snpEff_vcf=${name}.snpEff.vcf
haplo_vcf=${name}.snpEff.haplo.vcf
${snpEff} -Xmx4000m -v -canon hg38 ${name}.mutect.vcf > ${snpEff_vcf}
${snpEff} -Xmx4000m -v -canon hg38 ${name}.haplotypecaller.vcf > ${haplo_vcf}

${perl_exec} ${vcf2maf}/vcf2maf.pl --input-vcf ${snpEff_vcf} --output-maf ${name}.anno.TCGAMAF --ref-fasta ${ref} --filter-vcf 0 --vep-path ${vep} --vep-data ${vep_cache} --ncbi-build GRCh38 --retain-info ANN --retain-fmt GT,AD,AF,DP --tumor-id ${name}
${perl_exec} ${vcf2maf}/vcf2maf.pl --input-vcf ${haplo_vcf} --output-maf ${name}.haplo.TCGAMAF --ref-fasta ${ref} --filter-vcf 0 --vep-path ${vep} --vep-data ${vep_cache} --ncbi-build GRCh38 --retain-info DP,ANN --retain-fmt GT,AD,DP,GQ --tumor-id ${name}

vep_parse='/project/cpdlab/Scripts/SolidV2/solidv2_v1.2/methods/vep_parse.capture.py'
haplo_parse='/project/cpdlab/Scripts/SolidV2/solidv2_v1.2/methods/haplotype_parse.capture.py'

sample_name=`cut --complement -d'_' -f10- <<< $(basename ${read1}) | sed "s|_|\-|g"`
/project/cpdlab/anaconda2/bin/python ${vep_parse} ${name}.anno.TCGAMAF ${name}.parsed.tsv ${sample_name}
/project/cpdlab/anaconda2/bin/python ${haplo_parse} ${name}.haplo.TCGAMAF ${name}.haplo.tsv ${sample_name}

cat ${name}.parsed.tsv ${name}.haplo.tsv | awk -F'\t' '!seen[$3]++' >> ../Run_masterVarFinal.txt
final_anno=${name}.annotations.tsv

# Master Variant Table Xtreme

vep_simpletation='/project/cpdlab/Scripts/SolidV2/solidv2_v1.2/methods/vep_simpletation.py'
haplo_simpletation='/project/cpdlab/Scripts/SolidV2/solidv2_v1.2/methods/haplo_simpletation.py'

/project/cpdlab/anaconda2/bin/python ${vep_simpletation} ${name}.anno.TCGAMAF ${name}.vep_Xtreme.tsv ${sample_name}
/project/cpdlab/anaconda2/bin/python ${haplo_simpletation} ${name}.haplo.TCGAMAF ${name}.haplo_Xtreme.tsv ${sample_name}

cat ${name}.vep_Xtreme.tsv ${name}.haplo_Xtreme.tsv | awk -F'\t' '!seen[$3]++' >> ../MasterVariantTableXtreme.txt

python ${simpletation} -n ${name} -i ${vcf2tsv}  -o ${final_anno} -g ${gene_list}

#graphs
#cp /project/cpdlab/dev/ash/capture/methods/PerSamplePlots.R ${out_dir}
#Rscript PerSamplePlots.R
