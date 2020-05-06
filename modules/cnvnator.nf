reference = file(params.reference)
reference_dict = file(params.reference_dict)
reference_index = file(params.reference_index)
output_folder = file(params.output_folder)



process ReadMapping {
	publishDir "$params.output_folder/${sampleID}"
	container "hmartiniano/docker_cnvnator"
	afterScript "rm *"
	
	input:
	tuple val(sampleID), file(bam_file) 

	output:
	tuple val("${sampleID}"), path("${sampleID}.root") 

	"""
	cnvnator -root ${sampleID}.root -tree ${bam_file} -chrom 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y
	"""
	
}

process Histogram {
	publishDir "$params.output_folder/${sampleID}"
        container "hmartiniano/docker_cnvnator"
	afterScript "rm	*"

	input:
	tuple val(sampleID), file(root)
	file reference

	output:
	tuple  val("${sampleID}"), path("${sampleID}.root")

	"""
	cnvnator -root ${root} -his 100 -fasta ${reference}
	"""

}


process Statistics {
	publishDir "$params.output_folder/${sampleID}"
        container "hmartiniano/docker_cnvnator"
	afterScript "rm *"

	input:
	tuple val(sampleID), file(root) 

	output:
	tuple val("${sampleID}"), path("${sampleID}.root") 
	
	"""
	cnvnator -root ${root}  -stat 100
	"""
}


process RDSignalPartition {
	publishDir "$params.output_folder/${sampleID}"
        container "hmartiniano/docker_cnvnator"
	afterScript "rm *"

	input:
	tuple val(sampleID), file(root) 

	output:
        tuple val("${sampleID}"), path("${sampleID}.root")

	"""
	cnvnator -root ${root} -partition 100
	"""
}

process CNVCall {
	publishDir "$params.output_folder/${sampleID}"
        container "hmartiniano/docker_cnvnator"
	afterScript "rm *"

	input:
	tuple val(sampleID), file(root) 

	output:
	tuple val("${sampleID}"), file("${sampleID}_cnv.txt"), file("${sampleID}.root")
	
	"""
	cnvnator -root ${root} -call 100 > ${sampleID}_cnv.txt
	"""
}

process CnvnatorAnnotater {
        publishDir "$params.output_folder/${sampleID}"
        container "fredhutch/xenofilter:1.6"
        afterScript "rm *"
	
	input:
	tuple val(sampleID), file(cnv_file), file(root)

	output:
	file("${sampleID}_annotated_cnvnator.txt")

"""
#!/usr/bin/env Rscript

library(stringr)
library(Homo.sapiens)
library(GenomicRanges)
sampleID <- '${sampleID}'
cnv_file <- '${cnv_file}'

cnv <- read.csv(cnv_file, header = F, sep = '\\t')
colnames(cnv) <- c("CNV_type", "coordinates", "CNV_size", "normalized_RD", "e-val1", "e-val2", "e-val3", "e-val4", "q0")
positions <- str_split(cnv\$coordinates, ":")
cnv\$Chromosomes <- unlist(lapply(positions, '[[', 1))
start_stop <- str_split(unlist(lapply(positions, '[[', 2)), "-")
cnv\$Start_Position <- unlist(lapply(start_stop, '[[', 1))
cnv\$End_Position <- unlist(lapply(start_stop, '[[', 2))
cnv\$Tumor_Sample_Barcode <- sampleID
bed_df <-cnv[,10:12]
colnames(bed_df) <- c("chrom", "chromStart", "chromEnd")


cnv_intervals = makeGRangesFromDataFrame(bed_df)
mcols(cnv_intervals, level="within")\$CNV_type <- cnv\$CNV_type
mcols(cnv_intervals, level="within")\$norm_RD <- cnv\$normalized_RD
mcols(cnv_intervals, level="within")\$`e-val1` <- cnv\$`e-val1`
mcols(cnv_intervals, level="within")\$`e-val2` <- cnv\$`e-val2`
mcols(cnv_intervals, level="within")\$`e-val3` <- cnv\$`e-val3`
mcols(cnv_intervals, level="within")\$`e-val4` <- cnv\$`e-val4`
mcols(cnv_intervals, level="within")\$q0 <- cnv\$q0


gene_intervals <- genes(Homo.sapiens, columns="SYMBOL")

#query            #subject
complete_gene_hits <- findOverlaps(sort(gene_intervals), cnv_intervals , type = "within")

#these are the cnvs that had genes in them (most of them probably)
covered_cnv_intervals <- cnv_intervals[subjectHits(complete_gene_hits)]

#1 these genes are fully there
covered_gene_intervals <- sort(gene_intervals)[queryHits(complete_gene_hits)]

#add cnv columns to gene intervals (this is so manual)
mcols(covered_gene_intervals, level="within")\$CNV_type <- mcols(covered_cnv_intervals)\$CNV_type
mcols(covered_gene_intervals, level="within")\$norm_RD <- mcols(covered_cnv_intervals)\$norm_RD
mcols(covered_gene_intervals, level="within")\$`e-val1` <- mcols(covered_cnv_intervals)\$`e-val1`
mcols(covered_gene_intervals, level="within")\$`e-val2` <- mcols(covered_cnv_intervals)\$`e-val2`
mcols(covered_gene_intervals, level="within")\$`e-val3` <- mcols(covered_cnv_intervals)\$`e-val3`
mcols(covered_gene_intervals, level="within")\$`e-val4` <- mcols(covered_cnv_intervals)\$`e-val4`
mcols(covered_gene_intervals, level="within")\$q0 <- mcols(covered_cnv_intervals)\$q0

exon_intervals <- exons(Homo.sapiens, columns=c("EXONID", "SYMBOL"))

#now subtract gene_hits from cnv_intervals
leftover_intervals <- GenomicRanges::setdiff(cnv_intervals, covered_gene_intervals, ignore.strand=TRUE)

#Get exons that didn't go into genes neatly
leftover_exon_hits <- findOverlaps(leftover_intervals, exon_intervals)

#get the intervals remaining that map to exons
leftover_exon_intervals <- exon_intervals[subjectHits(leftover_exon_hits)]

#now overlap with original cnv data to get cnv columns
exon_cnv_hits <- findOverlaps(leftover_exon_intervals, cnv_intervals)
final_exon_intervals <- leftover_exon_intervals[queryHits(exon_cnv_hits)]
cnv_column_intervals <- cnv_intervals[subjectHits(exon_cnv_hits)]

mcols(final_exon_intervals, level = "within")\$CNV_type <- mcols(cnv_column_intervals)\$CNV_type
mcols(final_exon_intervals, level = "within")\$norm_RD <- mcols(cnv_column_intervals)\$norm_RD
mcols(final_exon_intervals, level="within")\$`e-val1` <- mcols(cnv_column_intervals)\$`e-val1`
mcols(final_exon_intervals, level="within")\$`e-val2` <- mcols(cnv_column_intervals)\$`e-val2`
mcols(final_exon_intervals, level="within")\$`e-val3` <- mcols(cnv_column_intervals)\$`e-val3`
mcols(final_exon_intervals, level="within")\$`e-val4` <- mcols(cnv_column_intervals)\$`e-val4`
mcols(final_exon_intervals, level="within")\$q0 <- mcols(cnv_column_intervals)\$q0

#get final                                
final_leftovers <- setdiff(leftover_intervals, final_exon_intervals, ignore.strand=TRUE)

#query              #subject
final_leftovers_cnv_hits <- findOverlaps(final_leftovers, cnv_intervals)
leftovers_cnv_intervals <- cnv_intervals[subjectHits(final_leftovers_cnv_hits)]
final_leftovers_intervals <- final_leftovers[queryHits(final_leftovers_cnv_hits)]

mcols(final_leftovers_intervals, level = "within")\$CNV_type <- mcols(leftovers_cnv_intervals)\$CNV_type
mcols(final_leftovers_intervals, level = "within")\$norm_RD <- mcols(leftovers_cnv_intervals)\$norm_RD
mcols(final_leftovers_intervals, level="within")\$`e-val1` <- mcols(leftovers_cnv_intervals)\$`e-val1`
mcols(final_leftovers_intervals, level="within")\$`e-val2` <- mcols(leftovers_cnv_intervals)\$`e-val2`
mcols(final_leftovers_intervals, level="within")\$`e-val3` <- mcols(leftovers_cnv_intervals)\$`e-val3`
mcols(final_leftovers_intervals, level="within")\$`e-val4` <- mcols(leftovers_cnv_intervals)\$`e-val4`
mcols(final_leftovers_intervals, level="within")\$q0 <- mcols(leftovers_cnv_intervals)\$q0

#1 Fully covered Genes
mcols(covered_gene_intervals, level = "within")\$EXONID <- NA
mcols(covered_gene_intervals, level = "within")\$GENE_COVERAGE <- "all"
#2 Exons covered by remaining copy number data
mcols(final_exon_intervals, level = "within")\$GENE_COVERAGE <- "partial"

#3 Copy number data that doesn't map to genes or exons
mcols(final_leftovers_intervals, level = "within")\$EXONID <- NA
mcols(final_leftovers_intervals, level = "within")\$SYMBOL<- NA
mcols(final_leftovers_intervals, level = "within")\$GENE_COVERAGE <- "none"

#I didn't want to import dplyr so this is kinda nested, basically just turns them to dataframes and then combines them
final_df <- rbind(rbind(as.data.frame(covered_gene_intervals), as.data.frame(final_exon_intervals)), as.data.frame(final_leftovers_intervals))
final_df\$Tumor_Sample_Barcode <- sampleID
#remove strand data since its sort of misleading
final_df\$strand <- NULL

final_df\$SYMBOL <- as.character(final_df\$SYMBOL)
final_df\$EXONID <- as.character(final_df\$EXONID)

write.csv(final_df, paste0(sampleID, '_annotated_cnvnator.txt'))

	"""

}

workflow cnvnator_wf {
	 
	 // needs to have sampleID and bamfile
	 get: bams_ch 
	 
	 main:
	 ReadMapping(bams_ch)
	 Histogram(ReadMapping.out, reference)
	 Statistics(Histogram.out)
	 RDSignalPartition(Statistics.out)
	 CNVCall(RDSignalPartition.out)
	 CnvnatorAnnotater(CNVCall.out)

	 emit:
	 raw_cnvnator = CNVCall.out
	 annotated_cnvnator = CnvnatorAnnotater.out
	   
}