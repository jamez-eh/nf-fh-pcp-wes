
nextflow.preview.dsl=2



process readCounts {
	container "fredhutch/ichorcna:3.6.2"
	
	input:
	tuple val(sampleID), path(bam_file) 
	
	output:
	file("${sampleID}.bam")

	"""
readCounter -b --window 1000000 --quality 20 --chromosome "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y" ${bam_file} > ${sampleID}.wig

	"""
	
}


/*
process ichorCNA {
        container "fredhutch/ichorcna:3.6.2"

	input:
	tuple val(sampleID), file(root)
	file reference

//	output:
//	tuple  val("${sampleID}"), path("${sampleID}.root")

	"""
#Rscript /path/to/ichorCNA/scripts/runIchorCNA.R --id tumor_sample \
# --WIG /path/to/tumor.wig --ploidy "c(2,3)" --normal "c(0.5,0.6,0.7,0.8,0.9)" --maxCN 5 \
#  --gcWig /path/to/ichorCNA/inst/extdata/gc_hg19_1000kb.wig \
#  --mapWig /path/to/ichorCNA/inst/extdata/map_hg19_1000kb.wig \
#  --centromere /path/to/ichorCNA/inst/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt \
#  --normalPanel /path/to/ichorCNA/inst/extdata/HD_ULP_PoN_1Mb_median_normAutosome_mapScoreFiltered_median.rds \
#  --includeHOMD False --chrs "c(1:22, \"X\")" --chrTrain "c(1:22)" \
#  --estimateNormal True --estimatePloidy True --estimateScPrevalence True \
#  --scStates "c(1,3)" --txnE 0.9999 --txnStrength 10000 --outDir ./
	"""

}
*/

workflow ichorCNA_wf {
	 
	 get: bams_ch 
	 	 
	 main:
	 readCounts(bams_ch)
	   
}

workflow {
	 
	 bams_ch = Channel
            .fromPath(params.input_csv)
            .splitCsv(header:true)
            .map{ row-> tuple(row.sampleID, row.bam) }
	 
	 main:   
	 ichorCNA_wf(bams_ch)
	 
	 publish:
	 ichorCNA_wf.out to :'${params.output_folder}/ichorCNA/'

}
	