nextflow.preview.dsl=2


params.run = "none"
params.bwa_index = null

reference = file(params.reference)
reference_index = file(params.reference_index)

reference_dict = file(params.reference_dict)
rear_index = file(params.rear_index)
indels_index = file(params.indels_index)

indels = file(params.indels)
rear = file(params.rear)
pdx_reference = file(params.pdx_reference)

input_csv = file(params.input_csv)
output_folder = file(params.output_folder)
reference = file(params.reference)
contig_dict = file(params.contig_dict)




/*
include './modules/cnvnator' params (
        reference: params.reference,
        reference_index: params.reference_index,
        reference_dict: params.reference_dict,
        output_folder: params.output_folder
)


include './modules/single_varfilter.nf' params (
        reference: params.reference,
        reference_index: params.reference_index,
        reference_dict: params.reference_dict,
        output_folder: params.output_folder,
	data_source: false,
	rear: params.rear,
	rear_index: params.rear_index,
	indels: params.indels,
	indels_index: params.indels_index
)
*/

//include CNVkit_wf from './modules/CNVkit.nf'
include gatkCNV_wf from './modules/gatkCNV.nf'
include mutect2_wf from './modules/mutect2.nf'



process mouse_sed {
        container "ubuntu"
	errorStrategy 'retry'
	maxRetries 20

	input:
	file pdx_reference 

	output:
	path("edited_${pdx_reference}") 
	
	"""
	sed 's/>chr./&_m/' ${pdx_reference} > edited_${pdx_reference}
	"""
}

process combine_fastas {
        container "ubuntu"
	errorStrategy 'retry'
	maxRetries 30

        input:
        file pdx_reference
	file reference

        output:
	path("combined_reference.fa") 
	
        """
	cat ${reference} ${pdx_reference} > combined_reference.fa
	"""

}

process bwa_index {
	container "fredhutch/bwa:0.7.17"
	errorStrategy 'retry'
	maxRetries 30

	input:
        file combined_reference

     	output:
	    file "*.{amb,ann,bwt,pac,sa,alt}"

        script:
        """
        bwa index $combined_reference
        """	
}



process bwa_mem {
	container "fredhutch/bwa:0.7.17"
  	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
  	maxRetries 50

	input:
	tuple val(sampleID), val(kitID), val(type), val(patient), file(R1), file(R2)
	file ref 
	file bwa_ind

	output:
	tuple val("${sampleID}"), val("${kitID}"), val("${type}"), val("${patient}"),  file("${sampleID}.sam")

	"""
	bwa mem -R '@RG\\tID:${params.run}\\tLB:hg38\\tPL:Illumina\\tPU:barcode\\tSM:${params.run}' -t 3  ${ref} ${R1} ${R2} > ${sampleID}.sam
	"""
}

process remove_pdx {
        container "ubuntu"
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
        maxRetries 100

        input:
	tuple val(sampleID), val(kitID), val(type), val(patient), file(hybrid_sam)

        output:
	tuple val("${sampleID}"), val("${kitID}"), val("${type}"), val("${patient}"),  file("${sampleID}_hum.sam")

	"""
	#awk '{ if (\$3  !="chr1_m" && \$3 != "chr2_m" && \$3 != "chr3_m" && \$3 != "chr4_m" && \$3 != "chr5_m" && \$3 != "chr6_m" && \$3 != "chr7_m" && \$3 != "chr8_m" && \$3 != "chr9_m" && \$3 != "chr1_m0" && \$3 != "chr1_m1" && \$3 != "chr1_m2" && \$3 != "chr1_m3" && \$3 != "chr1_m4" && \$3 != " chr1_m5" && \$3 != "chr1_m6" && \$3 != "chr1_m7" && \$3 != "chr1_m8" && \$3 != "chr1_m9" && \$3 != "chr2_m0" && \$3 != "chr2_m1" && \$3 != "chr2_m2" && \$3 != "chrX_m" && \$3 != "chrY_m" && \$3 != "chrM_m" ) print \$0}' ${hybrid_sam} > ${sampleID}_hum.sam
	awk '!/\tchr._m/' ${hybrid_sam} > ${sampleID}_hum.sam
	"""
}




process sam2fastq {
	container 'broadinstitute/gatk:4.1.4.1'
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
        maxRetries 100

	input:
	tuple val(sampleID), val(kitID), val(type), val(patient), file(sam_file)

	output:
        tuple val("$sampleID"), val("${kitID}"), val("${type}"), val("${patient}"),  file("${sampleID}_R1.fq.gz"), file("${sampleID}_R2.fq.gz")

        """
	gatk SamToFastq -I ${sam_file} \
			-F ${sampleID}_R1.fq.gz\
			-F2 ${sampleID}_R2.fq.gz \
			-FU ${sampleID}_unpaired.fq.gz
        """
}


process sam_to_bam {
	container "fredhutch/bwa:0.7.17-samtools-1.10"
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
        maxRetries 100

	input:
	tuple val(sampleID), val(kitID), val(type), val(patient), file(sam_file) 
	
	output:
	tuple val("${sampleID}"), val("${kitID}"), val("${type}"), val("${patient}"),  file("${sampleID}.bam") 

	"""
	samtools view -bhS ${sam_file} > ${sampleID}.bam
        """
} 




process bam_filtering {
	container "fredhutch/bwa:0.7.17-samtools-1.10"
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
        maxRetries 100

        input:
        tuple val(sampleID), val(kitID), val(type), val(patient), file(bam_file) 

        output:
        tuple val("${sampleID}"), val("${kitID}"), val("${type}"), val("${patient}"),  file("${sampleID}_q.bam")

        """
	samtools view -bhS -q 20 ${bam_file} > ${sampleID}_q.bam
        """
}


process picard_clean {
	container 'broadinstitute/gatk:4.1.4.1'
	errorStrategy 'ignore'

	input:
	tuple val(sampleID), val(kitID), val(type), val(patient), file(bam_q_file)

	output:
	tuple val("${sampleID}"), val("${kitID}"), val("${type}"), val("${patient}"),  file("${sampleID}_q_clean.bam")

	"""
	gatk --java-options "-Xmx30G" CleanSam -I ${bam_q_file} -O ${sampleID}_q_clean.bam
	"""
}


process sam_sort {
	container "fredhutch/bwa:0.7.17-samtools-1.10"
	errorStrategy 'ignore'

        input:
        tuple val(sampleID), val(kitID), val(type), val(patient), file(bam_q_clean_file)

        output:
	tuple val("${sampleID}"), val("${kitID}"), val("${type}"), val("${patient}"),  file("${sampleID}_q_clean_sorted.bam") 

        """
	samtools sort ${bam_q_clean_file} -o ${sampleID}_q_clean_sorted.bam
	"""
}

//7
/*
//Identifies duplicate reads. This tool locates and tags duplicate reads in a BAM or SAM file, where duplicate reads
//are defined as originating from a single fragment of DNA.  Duplicates can arise during sample preparation e.g. library
construction using PCR.  See also 
href='https://broadinstitute.github.io/picard/command-line-overview.html#EstimateLibraryComplexity'>EstimateLibraryComplexity
*/

process picard_duplicates {
	container 'broadinstitute/gatk:4.1.4.1'		
	memory '50 GB'
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
        maxRetries 100

	input:
        tuple val(sampleID), val(kitID), val(type), val(patient), file(bam_q_clean_sorted_file)

        output:
	tuple val("${sampleID}"), val("${kitID}"), val("${type}"), val("${patient}"),  file("${sampleID}_q_clean_sorted_rmdp.bam")

        """
        gatk MarkDuplicates -I ${bam_q_clean_sorted_file} -O ${sampleID}_q_clean_sorted_rmdp.bam -METRICS_FILE {sampleID}_metrics.txt  -REMOVE_DUPLICATES true -VALIDATION_STRINGENCY STRICT -ASSUME_SORTED true -CREATE_INDEX true
        """
}


//Base quality recalibration log-scale scoring
process gatk_baserecalibrator {
	container 'broadinstitute/gatk:4.1.4.1'
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
        maxRetries 100


        input:
        tuple val(sampleID), val(kitID), val(type), val(patient), file(sorted_bam)
        file ref
        file rear
	file reference_index
	file reference_dict
	file rear_index
	file indels_index

        output:
        tuple val("${sampleID}"), val("${kitID}"), val("${type}"), val("${patient}"),  val(patient), file("${sampleID}_recal_data.table") 

        """
        gatk BaseRecalibrator -R ${ref} -known-sites ${rear} -I ${sorted_bam} -O ${sampleID}_recal_data.table --java-options -Xmx8g 
  
        """
}

process gatk_printreads{
        container 'broadinstitute/gatk:4.1.4.1'
        publishDir "$params.output_folder/${sampleID}"
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
        maxRetries 100

        input:
        tuple val(sampleID), val(kitID), val(type), file(precal_data), val(patient), file(bam_file)
        file reference
        file reference_index
        file reference_dict

        output:
        tuple val("${sampleID}"), val("${kitID}"), val("${type}"), val("${patient}"), file("${sampleID}_bqsr.bam")

        """
	echo ${kitID}
	echo ${type}
	echo ${sampleID}
	echo ${precal_data}
	echo ${bam_file}
        gatk --java-options "-Xmx30G" ApplyBQSR -R ${reference} -I ${bam_file} -bqsr-recal-file ${precal_data} -O ${sampleID}_bqsr.bam
        """
}



workflow pdx_align {
	 
	 get: fqs_ch

	 main:
	   mouse_sed(pdx_reference)
	   combine_fastas(mouse_sed.out, reference)
	   bwa_index(combine_fastas.out)
	   bwa_mem(fqs_ch, combine_fastas.out, bwa_index.out)
	   remove_pdx(bwa_mem.out)
	   sam2fastq(remove_pdx.out)
	   
	 emit:
	   human_fqs = sam2fastq.out
	   mixed_bams = bwa_mem.out
}

workflow vanilla_align {

	get: human_fq_ch
	
	main:
	  bwa_index(reference)
	  bwa_mem(human_fq_ch, reference, bwa_index.out)
	  sam_to_bam(bwa_mem.out)
	  bam_filtering(sam_to_bam.out)
	  picard_clean(bam_filtering.out)
	  sam_sort(bam_filtering.out)
	  picard_duplicates(sam_sort.out)
	  gatk_baserecalibrator(picard_duplicates.out, reference, rear, reference_index, reference_dict, rear_index, indels_index)
	  gatk_printreads(gatk_baserecalibrator.out.join(picard_duplicates.out, by : [0,1,2]), reference, reference_index, reference_dict)
	
	emit:
	  raw_bams = sam_to_bam.out
	  recal_bams = gatk_printreads.out
}


workflow {

	 reference = file(params.reference)
	 reference_index = file(params.reference_index)

	 reference_dict = file(params.reference_dict)
	 rear_index = file(params.rear_index)
	 indels_index = file(params.indels_index)

	 indels = file(params.indels)
	 rear = file(params.rear)
	 pdx_reference = file(params.pdx_reference)

	 input_csv = file(params.input_csv)
	 output_folder = file(params.output_folder)
	 reference = file(params.reference)
	 contig_dict = file(params.contig_dict)
	 common_variants = file(params.common_variants)
	 common_variants_index = file(params.common_variants_index)

	 

	fqs_ch = Channel
	    .fromPath(params.input_csv)
	    .splitCsv(header:true)
	    .map{ row-> tuple(row.sampleID, row.kitID, row.type, row.patient, file(row.R1), file(row.R2)) }

	if (params.cnv || params.mutect2) {
	beds_ch = Channel
       	     .fromPath(params.input_beds)
             .splitCsv(header:true)
             .map{ row-> tuple(row.kitID, file(row.capture_bed)) }
	}
        bam_bam = Channel.empty()


	main:



	fqs_ch.branch {
                Normal : it[2] == 'Normal'
                Tumor : it[2] == 'Tumor'
               }.set { fqs_branched }


        pdx = Channel.empty()
	removed = Channel.empty()

	if(params.pdx) {
		pdx_align(fqs_branched.Tumor)
		pdx = pdx_align.out
		final_fqs = pdx.human_fqs
	}
	else {
	       final_fqs = fqs_branched.Tumor
	}
	
	
	fqs_mixed = final_fqs.mix(fqs_branched.Normal)

	vanilla_align(fqs_mixed)

	bam_bam = vanilla_align.out



	//gatk_single_snp_wf(vanilla_align.out.recal_bams)


	//CNVkit_wf(vanilla_align.out.recal_bams, beds_ch, reference, reference_dict, reference_index, common_variants, common_variants_index)
	
	if(params.cnv){
        gatkCNV_wf(beds_ch, reference, reference_dict, reference_index, contig_dict, vanilla_align.out)
	}
	if(params.mutect2){
	mutect2_wf(vanilla_align.out.recal_bams, beds_ch, reference, reference_dict, reference_index, common_variants, common_variants_index)
	}
	

	publish:
	pdx.mixed_bams to : "${params.output_folder}/mixed_bams/"
	vanilla_align.out.raw_bams to: "${params.output_folder}raw/"
	vanilla_align.out.recal_bams to: "${params.output_folder}recal_bams/"
//	gatkCNV_wf.out.plots to : "./lucap_plots/"	  
//	gatkCNV_wf.out.plots to : "${params.output_folder}denoised/"	  
 //	gatkCNV_wf.out.modelsegs to : "${params.output_folder}modelsegs"	  
//	gatkCNV_wf.out.modelPlots to : "./lucap_modplots/"	  
//	gatkCNV_wf.out.func to : './func/'
//	gatkCNV_wf.out.calledSegs to : '/fh/scratch/delete90/nelson_p/james/funcotating/calledSegs/'
//	gatkCNV_wf.out.modelsegs to : '/fh/scratch/delete90/nelson_p/james/funcotating/modelSegs/'
//	gatkCNV_wf.out.download to : '/fh/scratch/delete90/nelson_p/james/funcotating/data/'
//	mutect2_wf.out.annotated to : '/fh/scratch/delete90/nelson_p/james/funcotating/mutect2/MAFS/'
//	mutect2_wf.out.vcf to : '/fh/scratch/delete90/nelson_p/james/funcotating/mutect2/VCF/'
//	mutect2_wf.out.filteredVCF to : '/fh/scratch/delete90/nelson_p/james/funcotating/mutect2/filteredVCF/'
//	gatkCNV_wf.out.segs to : '/fh/scratch/delete90/nelson_p/james/funcotating/segs/'
//	gatkCNV_wf.out.gistic to : '/fh/scratch/delete90/nelson_p/james/funcotating/gistic/'
}


