nextflow.preview.dsl=2

include	bwa_index as host_index from './alignment.nf'
include bwa_index as graft_index from './alignment.nf'
include bwa_mem as host_align from './alignment.nf'
include bwa_mem as graft_align from './alignment.nf'
include sam_to_bam as host_samtobam from './alignment.nf'
include sam_to_bam as graft_samtobam from './alignment.nf'

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


process xenofilter {
        container 'fredhutch/xenofilter'
        afterScript "rm -r *"

        input:
        tuple val(sampleID), file(graft_bam), file(host_bam)

        output:
	tuple val("${sampleID}"), file("${sampleID}_xenofiltered.bam"),	file("${sampleID}_xenofilteredXenofilteR.log")


        """
	Rscript xenofiltR.R ${sampleID} ${graft_bam} ${host_bam}
        """
}

workflow xenofilter_wf {

         get: fqs_ch
	      reference
	      pdx_reference
	      ref_name
	      pdx_ref_name
	      
         main:
                 graft_index(reference)
                 host_index(pdx_reference)
                 graft_align(fqs_ch, reference, graft_index.out, ref_name)
                 host_align(fqs_ch, pdx_reference, host_index.out, pdx_ref_name)
                 graft_samtobam(graft_align.out)
                 host_samtobam(host_align.out)
                 xenofilter(graft_samtobam.out.join(host_samtobam.out))

	emit:
		host_bams = xenofilter.out
}

workflow {
	
        fqs_ch = Channel
            .fromPath(params.input_csv)
            .splitCsv(header:true)
            .map{ row-> tuple(row.sampleID, file(row.R1), file(row.R2)) }


	
	main:
		xenofilter_wf(fqs_ch, params.reference, params.pdx_reference, params.ref_name, params.pdx_ref_name)

	publish:	
		xenofilter_wf.out.host_bams to : "${params.output_folder}/test/"

}