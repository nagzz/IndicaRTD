# Raw sequencing data curation with AdapterRemoval tool:

AdapterRemoval --file1 $INFILE1 --file2 $INFILE2 --threads 15 --adapter-list adap.fa --minquality 25 --basename ${FName}

# Downloaded the Oryza sativa ssp. indica genome (ASM465v1.43) fasta file, gff3 and GTF files from Ensembl plant database.

# Quality filtered sequencing data was used for STAR mapping :

### STAR genome index generated with the following parameters using the downloaded genome files:
################################################################################################

STAR --runThreadN 10 --runMode genomeGenerate --genomeDir /home/ns43567/CEIB_rice_project/JHI/ensembl_indica --genomeFastaFiles Oryza_indica.ASM465v1.dna.toplevel.fa --sjdbGTFfile ./Oryza_indica.ASM465v1.43.gtf --sjdbOverhang 100 --genomeChrBinNbits 8 --genomeSAindexNbases 13 --limitGenomeGenerateRAM 32000000000 --genomeSAsparseD 5

### First-pass STAR mapping - Mapping the quality filtered reads to the genome file using STAR
#############################################################################################

STAR --runMode alignReads --genomeDir /home/ns43567/CEIB_rice_project/JHI/ensembl_indica --readFilesType Fastx --readFilesIn $INFILE --sjdbOverhang 100 --outSAMprimaryFlag AllBestScore --outFilterMismatchNmax 2/0 --outSJfilterCountTotalMin 10 5 5 5 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --alignIntronMin 60 --alignIntronMax 6000 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix /home/ns43567/CEIB_rice_project/transcriptome/NGS_P5366/pe/${FName}_ --outTmpDir /home/ns43567/CEIB_rice_project/transcriptome/NGS_P5366/pe/${FName}_STAR

#STAR output file sj.out.tab file specifications:

#column 1: chromosome
#column 2: first base of the intron (1-based)
#column 3: last base of the intron (1-based)
#column 4: strand (0: undefined, 1: +, 2: -)
#column 5: intron motif: 0: non-canonical; 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5: AT/AC, 6: GT/AT
#column 6: 0: unannotated, 1: annotated (only if splice junction database is used)
#column 7: number of uniquely mapping reads crossing the junction
#column 8: number of multi-mapping reads crossing the junction
#column 9: maximum spliced alignment overhang

### Novel SJ filteration after first-pass STAR mapping (unannotated sj):
########################################################################

cat *.tab | awk '($5 > 0 && $7 >= 1 && $6==0)' | cut -f1-6 | sort | uniq > NovelSJ_paired.tab

### STAR genome index generated with the following parameters for the genome files along with the NovelSJ_paired.tab:
#####################################################################################################################

STAR --runThreadN 10 --runMode genomeGenerate --genomeDir /home/ns43567/CEIB_rice_project/JHI/ensembl_indica --genomeFastaFiles Oryza_indica.ASM465v1.dna.toplevel.fa --sjdbGTFfile ./Oryza_indica.ASM465v1.43.gtf --sjdbOverhang 100 --genomeChrBinNbits 8 --genomeSAindexNbases 13 --limitGenomeGenerateRAM 32000000000 --genomeSAsparseD 5 --sjdbFileChrStartEnd /home/ns43567/CEIB_rice_project/transcriptome/NGS_P5366/pe/StarMapping_paired/NovelSJ_paired.tab

### Second-pass STAR mapping - Mapping the quality filtered reads to the genome file using STAR:
################################################################################################

STAR --runThreadN 2 --runMode alignReads --genomeDir /home/ns43567/CEIB_rice_project/JHI/ensembl_indica --readFilesType Fastx --readFilesIn $INFILE --sjdbOverhang 100 --outSAMprimaryFlag AllBestScore --outFilterMismatchNmax 0 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outSJfilterCountTotalMin 10 5 5 5 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --alignIntronMin 60 --alignIntronMax 6000 --alignMatesGapMax 400 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix /home/ns43567/CEIB_rice_project/transcriptome/NGS_P5366/pe/StarMappingpaired2ndPass/${FName}_paired_2ndPass_

### Total number of SJ after second-pass STAR mapping (both annotated, unannotated and cannonica sj):
##################################################################################################

cat *.tab | awk '($5 > 0 && $7 >= 1)' |  cut -f1-6 | sort | uniq |wc -l

### Assembling the mapped reads with 3 different assemblers Scallop, Stringtie and Cufflinks with the default parameters
########################################################################################################################

scallop -i $INFILE -o ${FName}_Scallop.gtf
cufflinks -g Oryza_indica.ASM465v1.43.gtf -o ${FName}_cufflinks $INFILE
stringtie2 $INFILE -o ${FName}_stringtie2.gtf -G ~/Oryza_indica.ASM465v1.43.gtf -t -f 0

###Count number of genes in the GTF file:
#########################################

for i in *.gtf; do printf $i":" && awk 'BEGIN{FS="\t"}{print $9}' $i| awk 'BEGIN{FS=";"}{print $1}'| sed 's/gene_id "//g;s/"//g;/^$/d'|sort -u | wc -l; done

### Merging the individual assemblies of all samples for each assembler with cuffmerge, stringtie and taco with default parameters
#################################################################################################################################

stringtie --merge -p 10 -o StringtieMergedScallop.gtf -F 0 -T 0 -f 0 -g 0 -i ScallopGtf.txt

stringtie --merge -p 10 -o StringtieMergedCufflinks.gtf -F 0 -T 0 -f 0 -g 0 -i CufflinksGtf.txt

stringtie --merge -p 10 -o StringtieMergedStringtie.gtf -F 0 -T 0 -f 0 -g 0 -i StringtieGtf.txt

#-F <min_fpkm>    : minimum input transcript FPKM to include in the merge (default: 1.0)
#-T <min_tpm>     : minimum input transcript TPM to include in the merge (default: 1.0)
#-f <min_iso>     : minimum isoform fraction (default: 0.01)
#-g <gap_len>     : gap between transcripts to merge together (default: 250)
#-i               : keep merged transcripts with retained introns; by default these are not kept unless there is strong evidence for them

#### Reference Transcriptome data (RTD) construction, creating the splice junction database (sj.db) and filtration of RTD with high-confident SJs
#################################################################################################################################################

cat StringtieMerged*.gtf > stringtiemergedgtflist.txt

./stringtie --merge -p 10 -o stringtiemerge3assemblers.gtf -F 0 -T 0 -g 0 -i stringtiemergedgtflist.txt

Rscript --vanilla ../multifilesgread.R stringtiemerge3assemblers.gtf > stringtiemerge3assemblers_intron.txt

awk '{print $11"\t"$1"_"$5"_"$2"_"$3}' stringtiemerge3assemblers_intron.txt > stringtiemerge3assemblers_intron_modified.txt

cat *sj.out.tab | awk '($5 > 0 && $7 >= 5 && $9 >= 10)' | cut -f1-4 | sort |uniq > sj.filtered.out.tab

#column 5: intron motif: 0: non-canonical; 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5: AT/AC, 6: GT/AT
#column 7: number of uniquely mapping reads crossing the junction
#column 9: maximum spliced alignment overhang

awk '$4=="1" { print $1"_+_"$2"_"$3} $4=="2" { print $1"_-_"$2"_"$3}' sj.filtered.out.tab > sj.db

awk 'FNR==NR{a[$1];next}($2 in a){print}' sj.filtered.out.modified.txt stringtiemerge3assemblers_intron_modified.txt > tmp.txt

awk '{print $1}' stringtiemerge3assemblers_intron_modified.txt | sort | uniq -c | sed 's/^ *//g' > transcriptID_count.txt

awk '{print $1}' tmp.txt| sort | uniq -c | sed 's/^ *//g'| sed 's/ /_/g'|sort | comm -12 - transcriptID_count.txt > matchedtranscript.txt

awk 'BEGIN{FS="_"}{print $2}' matchedtranscript.txt |sort -u|fgrep -w -f - stringtiemerge3assemblers.gtf > matchedtranscript.gtf

### Extract & merge single exon transcripts:
############################################

tail -n +3 stringtiemerge3assemblers.gtf |awk '{print $12}'|sed 's/"//g;s/;//'|sort|uniq -c > stringtiemerge3assemblers_transcriptID_count.txt

sed 's/^[ ]*//g' stringtiemerge3assemblers_transcriptID_count.txt |awk 'BEGIN{FS=" "}{if ($1==2) print $2}'|fgrep -w -f - stringtiemerge3assemblers.gtf > stringtiemerge3assemblers_single_exon.gtf

stringtie --merge -p 10  -G /home/ns43567/CEIB_rice_project/JHI/ensembl_indica/Oryza_indica.ASM465v1.43.gtf -o IndicaRTD.gtf -F 0 -T 0 -g 0 -i matchedtranscript.gtf stringtiemerge3assemblers_single_exon.gtf
