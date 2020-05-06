nextflow.preview.dsl=2

include BedToIntervalList from './gatkCNV.nf'

include PreprocessIntervals from './gatkCNV.nf'

include samtoolsIndex from './gatkCNV.nf'

include Funcotator from './single_varfilter.nf'


process DownloadData {
container 'broadinstitute/gatk:4.1.4.1'
errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
maxRetries 100
                
	output:
        path funcotator_dataSource

                """
                gatk FuncotatorDataSourceDownloader --somatic --validate-integrity --extract-after-download
                rm funcotator_dataSource*gz
                mv funcotator_dataSource* funcotator_dataSource
                """

}


process samtoolsRemoveSecondary {
        container "fredhutch/bwa:0.7.17"
        errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
        maxRetries 100


	input:
	tuple val(sampleID), val(kitID), val(type), file(bam_file)

	output:
	tuple val("${sampleID}"), val("${kitID}"), val("${type}"), file("${sampleID}_primary.bam")


	"""
	samtools view -bh -f 0 -F 256 ${bam_file} > ${sampleID}_primary.bam
	"""
}

process mutect2_normal_only {
        container 'broadinstitute/gatk:4.1.5.0'
        errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
        maxRetries 10
	
	input:
	tuple val(sampleID), val(kitID), val(type), file(bam_file), file(bam_index)
	file reference
        file reference_dict
        file reference_index

	output:
	tuple val("${sampleID}"), val("${kitID}"), val("${type}"), file("${sampleID}.vcf.gz")
	
	"""
	gatk Mutect2 --java-options "-Xmx30G" \
	-R ${reference} \
	-I ${bam_file} \
	-O ${sampleID}.vcf.gz

	"""

}




process GenomicsDBImport {
	container 'broadinstitute/gatk:4.1.5.0'
        errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
        maxRetries 100

	input:
	tuple val(kitID), file(interval), file(normal_list)
	file reference
	file reference_dict
	file reference_index

	output:
	tuple val("${kitID}"), path("${kitID}_db")	

"""
#!/usr/bin/python

import os
os.system('for F in *.vcf.gz ; do   tabix -f -p vcf \${F}  ; done')


files = '${normal_list}'
names = files.replace('.vcf.gz', '\t')
name_list = names.split(' ')
file_list = files.split(' ')
res = [i + j for i, j in zip(name_list, file_list)]
outF = open("map.txt", "w")

for line in res:
    outF.write(line)
    outF.write("\\n")

outF.close()

files = files.replace(' ', ' -V ')
     
script1 = 'gatk GenomicsDBImport  -R ${reference} -L ${interval} --genomicsdb-workspace-path ${kitID}_db'
       
script_f = script1 + ' --sample-name-map map.txt --merge-input-intervals'

os.system(script_f)
	   
"""
}


process CreateSomaticPanelOfNormals {
        container 'broadinstitute/gatk:4.1.5.0'
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
        maxRetries 100

	input:
	tuple val(kitID), path(db)
        file reference
        file reference_dict
        file reference_index

	output:
	tuple val("${kitID}"), file("${kitID}.vcf.gz"), file("${kitID}.vcf.gz.tbi")
	
	"""
	gatk CreateSomaticPanelOfNormals -R ${reference} -V gendb://${db} -O ${kitID}.vcf.gz
	"""
} 


process	mutect2_tumor_only {
	container 'broadinstitute/gatk:4.1.5.0'
errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
maxRetries 100


	input:
	tuple val(kitID), path(normals), path(normal_index), val(sampleID), file(bam_file), file(bam_index)	
	file reference
	file reference_dict
	file reference_index
	file common_variants
	file common_variants_index


	output:
	tuple val("${sampleID}"), file("${sampleID}.vcf.gz"),  file("${sampleID}.vcf.gz.tbi"), file("${sampleID}.vcf.gz.stats")
	"""
  	gatk  --java-options "-Xmx30G" Mutect2 \
  	     -R ${reference} \
  	     -I ${bam_file} \
	     --panel-of-normals ${normals} \
	     --germline-resource ${common_variants} \
	     --min-base-quality-score 20 \
	     --pcr-indel-model AGGRESSIVE \
	     --callable-depth 14 \
	     --minimum-allele-fraction 0.2 \
	     --base-quality-score-threshold 20 \
    	     -O ${sampleID}.vcf.gz

	"""
}


process FilterMutectCalls {
        container 'broadinstitute/gatk:4.1.4.1'
        errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
maxRetries 100

        input:
        tuple val(sampleID), file(vcf),  file(vcf_index), file(stats)
        file reference
        file reference_dict
        file reference_index

        output:
	tuple val("${sampleID}"), file("${sampleID}_filtered.vcf.gz")

	"""
	gatk FilterMutectCalls \
	     -R ${reference} \
	     -V ${vcf} \
	     -O ${sampleID}_filtered.vcf.gz
	"""
}



workflow mutect2_wf {
	 take:
	 bams_ch
	 beds_ch
	 reference
	 reference_dict
	 reference_index
	 common_variants
	 common_variants_index

	 main:
	 samtoolsRemoveSecondary(bams_ch)
	 samtoolsIndex(samtoolsRemoveSecondary.out)
	 samtoolsIndex.out.branch {
                Normal : it[2] == 'Normal'
                Tumor : it[2] == 'Tumor'
               }.set { bams_branched }

	 mutect2_normal_only(bams_branched.Normal, reference, reference_dict, reference_index)
	 BedToIntervalList(beds_ch, reference, reference_dict)
	 
	 normal_vcfs = mutect2_normal_only.out.map{ r -> [r[1], r[3]] }.groupTuple()
	 grouped_normals = BedToIntervalList.out.cross(normal_vcfs).map{ r -> [r[0][0], r[0][1], r[1][1]]}
	 

	 GenomicsDBImport(grouped_normals, reference, reference_dict, reference_index)
	 CreateSomaticPanelOfNormals(GenomicsDBImport.out, reference, reference_dict, reference_index)
	 
	 
	 tumor_bams = bams_branched.Tumor.map{ r -> [r[1],r[0],r[3], r[4]] }


	 for_mutect = CreateSomaticPanelOfNormals.out.cross(tumor_bams).map{ r -> [r[0][0], r[0][1],r[0][2], r[1][1], r[1][2], r[1][3]]}
	 mutect2_tumor_only(for_mutect, reference, reference_dict, reference_index, common_variants, common_variants_index)
	 
	 FilterMutectCalls(mutect2_tumor_only.out, reference, reference_dict, reference_index)
	 
	 DownloadData()
	 Funcotator(FilterMutectCalls.out, reference, reference_index, reference_dict, DownloadData.out)

	 emit:
	 vcf = mutect2_tumor_only.out
	 filteredVCF = FilterMutectCalls.out
	 annotated = Funcotator.out 	 

}

