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
	tuple val(sampleID), val(kitID), val(type), val(patientID), file(bam_file)

	output:
	tuple val("${sampleID}"), val("${kitID}"), val("${type}"), val("${patientID}"), file("${sampleID}_primary.bam")


	"""
	samtools view -bh -f 0 -F 256 ${bam_file} > ${sampleID}_primary.bam
	"""
}
//       tuple val("${kitID}"), val("${patientID}"), val("${type}"), val("${sampleID}"), file("${bam}"), file("${bam}.bai")

process mutect2_normal_only {
        container 'broadinstitute/gatk:4.1.5.0'
        errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
        maxRetries 10
	
	input:
	tuple val(sampleID), val(kitID), val(type), val(patientID), file(bam_file), file(bam_index)
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




	input:
	tuple val(patientID), val(kitID), path(normals), path(normal_index), val(sampleID), file(bam_file), file(bam_index)	
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

process mutect2_matched_normal {
        container 'broadinstitute/gatk:4.1.5.0'
//errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
//maxRetries 100


        input:
        tuple val(patientID), val(kitID), path(normals), path(normal_index), val(sampleID), file(bam_file), file(bam_index), val(normal_sampleID), file(normal_bam_file), file(normal_index)
        file reference
        file reference_dict
        file reference_index
        file common_variants
        file common_variants_index


        output:
        tuple val("${sampleID}"), file("${sampleID}.vcf.gz"),  file("${sampleID}.vcf.gz.tbi"), file("${sampleID}.vcf.gz.stats")
        """

	gatk AddOrReplaceReadGroups \
       	-I ${bam_file} \
       	-O ${sampleID}.bam \
       	--RGID none \
       	--RGLB none \
       	--RGPL none \
       	--RGPU none \
       	--RGSM ${sampleID}

        gatk AddOrReplaceReadGroups \
        -I ${normal_bam_file} \
        -O ${normal_sampleID}.bam \
        --RGID none \
        --RGLB none \
        --RGPL none \
        --RGPU none \
        --RGSM ${normal_sampleID}
	
	echo *
       samtools index ${sampleID}.bam
       samtools index ${normal_sampleID}.bam

        gatk  --java-options "-Xmx30G" Mutect2 \
             -R ${reference} \
             -I ${sampleID}.bam \
	     -I ${normal_sampleID}.bam \
	     -normal ${normal_sampleID} \
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

process annotateVariants {
        container 'broadinstitute/gatk:4.1.4.1'


        input:
        tuple val(sampleID), file(vcf)
        path reference
        path reference_dict
        path reference_index
	path clinvar

        output:
        tuple val("${sampleID}"), file("${sampleID}_clinvar.vcf.gz")

        """
	 gatk IndexFeatureFile -I ${clinvar}
	 gatk IndexFeatureFile -I ${vcf}
	 gatk VariantAnnotator \
   	 		  -R ${reference} \
   			  -V ${vcf} \
   			  --output ${sampleID}_clinvar.vcf.gz \
			  --resource:ClinVar ${clinvar} \
			   -E ClinVar.AF_ESP \
			   -E ClinVar.AF_EXAC \
			   -E ClinVar.AF_TGP \
			   -E ClinVar.ALLELEID \
			   -E ClinVar.CLNDNINCL \
			   -E ClinVar.CLNDISDBINCL \
			   -E ClinVar.CLNHGVS \
			   -E ClinVar.CLNSIGINCL \
			   -E ClinVar.CLNVI \
			   -E ClinVar.CLNVC \
			   -E ClinVar.MC \
			   -E ClinVar.ORIGIN \
			   -E ClinVar.AF_ESP \
			   -E ClinVar.AF_EXAC \
			   -E ClinVar.AF_TGP \
			   -E ClinVar.ALLELEID \
			   -E ClinVar.CLNDNINCL \
			   -E ClinVar.CLNDISDBINCL \
			   -E ClinVar.CLNHGVS \
			   -E ClinVar.CLNSIGINCL \
			   -E ClinVar.CLNVI \
			   -E ClinVar.CLNVC \
			   -E ClinVar.MC \
			   -E ClinVar.ORIGIN \

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
	 clinvar


	 main:
	 samtoolsRemoveSecondary(bams_ch)
	 samtoolsIndex(samtoolsRemoveSecondary.out)

//       tuple val("${kitID}"), val("${patientID}"), val("${type}"), val("${sampleID}"), file("${bam}"), file("${bam}.bai")
	 samtoolsIndex.out.branch {
                Normal : it[2] == 'Normal'
                Tumor : it[2] == 'Tumor'
               }.set { bams_branched }


//        tuple val(sampleID), val(kitID), val(type), val(patientID), file(bam_file), file(bam_index)

	 mutect2_normal_only(bams_branched.Normal.map{r ->[r[3], r[0], r[2], r[1], r[4], r[5]] }, reference, reference_dict, reference_index)
	 BedToIntervalList(beds_ch, reference, reference_dict)
	 
	 normal_vcfs = mutect2_normal_only.out.map{ r -> [r[1], r[3]] }.groupTuple()
	 grouped_normals = BedToIntervalList.out.cross(normal_vcfs).map{ r -> [r[0][0], r[0][1], r[1][1]]}
	 

	 GenomicsDBImport(grouped_normals, reference, reference_dict, reference_index)

	 
	 CreateSomaticPanelOfNormals(GenomicsDBImport.out, reference, reference_dict, reference_index)
	 
	 tumor_bams = bams_branched.Tumor.map{ r -> [r[0],r[1],r[3], r[4], r[5]] } // remove type

	 
	 for_mutect = CreateSomaticPanelOfNormals.out.cross(tumor_bams).map{ r -> [r[1][1], r[0][0], r[0][1],r[0][2], r[1][2], r[1][3], r[1][4]]} //
	 
	 

	  patient_join = for_mutect.join(bams_branched.Normal.map{r -> [r[1], r[3], r[4], r[5]]}, remainder: true).filter { it[1] != null }
	  
	  patient_join.view()
	  patient_join.branch{
		Not : it [7] == null
		matched : it[7] != null
	  }.set { patient_branch }
	  
//	  patient_branch.Not.view()
	 unmatched_mutect2 = mutect2_tumor_only(patient_branch.Not.map{r->[r[0], r[1],r[2],r[3],r[4],r[5],r[6]]}, reference, reference_dict, reference_index, common_variants, common_variants_index) 
	 

      	 matched_mutect2 = mutect2_matched_normal(patient_branch.matched, reference, reference_dict, reference_index, common_variants, common_variants_index) 
	 
	 mutect2_calls = unmatched_mutect2.mix(matched_mutect2)
	 FilterMutectCalls(mutect2_calls, reference, reference_dict, reference_index)
	 annotateVariants(FilterMutectCalls.out, reference, reference_dict, reference_index, clinvar)
	 
	 DownloadData()
	 Funcotator(FilterMutectCalls.out, reference, reference_index, reference_dict, DownloadData.out)

	 emit:
	 vcf = mutect2_tumor_only.out
	 filteredVCF = FilterMutectCalls.out
	 annotated = Funcotator.out 	 

}

