### Ensembl RTD GTF file preparation for JCC analysis:
###################################################

#Preparing input gtf file:

gffread -w Oryza_indica.ASM465v1.43.fa -g /home/ns43567/CEIB_rice_project/JHI/ensembl_indica/Oryza_indica.ASM465v1.dna.toplevel.fa Oryza_indica.ASM465v1.43.gtf

python removeNfromfas.py Oryza_indica.ASM465v1.43.fa > Oryza_indica.ASM465v1.43.withoutN.fa

grep '>' Oryza_indica.ASM465v1.43.withoutN.fa|sed 's/-.*$//;s/^>//'|fgrep -w -f - Oryza_indica.ASM465v1.43.gtf > Oryza_indica.ASM465v1.43.withoutN.gtf

awk '{if($3=="transcript" || $3=="exon") print $0}' Oryza_indica.ASM465v1.43.withoutN.gtf > tmp.gtf

#strand fix

awk '{if ($7=="+" || $7=="-") print $0}' tmp.gtf > EnsemblRTD.gtf
