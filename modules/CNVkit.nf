nextflow.preview.dsl=2

params.hello = 'yes'
capture1 = 'v2'
capture2 = 'v3'
capture3 = 'v3+UTR'
params.refFlat = 's3://fh-pi-nelson-p/james/references/hg38/refFlat.txt.gz'

process target_CNVkit {
	container 'etal/cnvkit'

	input:
	tuple val(kitID), path(bed_file)
	output:
		

	"""
	cnvkit.py target ${bed_file} --annotate refFlat.txt --split -o ${kitID}_target.bed
	"""

}


process access_CNVkit {
	container 'etal/cnvkit'

	input:
	output:

	"""
	cnvkit.py access hg38.fa -x excludes.bed -o access-excludes.hg38.bed
	"""

}

process tag {
        container 'etal/cnvkit'
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
        maxRetries 100

        input:
	tuple val(sampleID), val(kitID), val(type),  file(bam)


        output:
	tuple val("${sampleID}"), val("${type}"), val("${kitID}"), file("${sampleID}_${type}.bam")


        """
	mv ${bam} ${sampleID}_${type}.bam
        """


}

//
// 	     --access data/access-5kb-mappable.hg19.bed \
process batch_CNVkit {
	container 'etal/cnvkit'
        errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
        maxRetries 100

	cpus 10


	input:
	tuple val(kitID), file(bed), file(bams)
	path reference
	path refFlat
	
	output:
	tuple val("${kitID}"), path(${kitID}_cnvkit)

	"""
	cnvkit.py batch *Tumor.bam \
            --normal *Normal.bam \
	    --targets ${bed} \
	    --annotate ${refFlat} \
   	     --fasta ${reference} \
    	     --output-reference ${kitID}.cnn \
	     --output-dir ${kitID}_cnvkit / \
    	     --diagram --scatter \
	     -p ${task.cpus} \
	     -y \
	     --drop-low-coverage \

	"""

}

workflow CNVkit_wf {

         take:
         bams_ch
         beds_ch
         reference
	 refFlat
		
	 
	 main:
	 
	 tagged_bams = tag(bams_ch)

	 kit_id_bams =  tagged_bams.map { r -> [r[2], r[1], r[0], r[3]] }.groupTuple().map { r -> [r[0], r[3]] }
	 
	 beds_and_bams = beds_ch.cross(kit_id_bams).map{ r -> [r[0][0], r[0][1], r[1][1]] }
	
	 batch_CNVkit(beds_and_bams, reference, refFlat)

//	 batch_CNVkit.out.view()


	emit: 
	      cnvkit_batch = batch_CNVkit.out
	      
}


workflow {
	        
	 beds_ch = Channel
            .fromPath(params.input_beds)
            .splitCsv(header:true)
            .map{ row-> tuple(row.kitID, file(row.capture_bed)) }

	 bams_ch = Channel
	    .fromPath(params.input_csv)
	    .splitCsv(header:true)
	    .map{ row-> tuple(row.sampleID, row.kitID, row.type, file(row.bam)) }

	reference = file(params.reference)
	refFlat = file(params.refFlat)

	 main:

	 CNVkit_wf(bams_ch, beds_ch, reference, refFlat)


}