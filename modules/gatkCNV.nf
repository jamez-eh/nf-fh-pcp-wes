//process CreateSequenceDictionaryTool (make a utils file)!!!!!
nextflow.preview.dsl=2
params.hello = 'yes'
capture1 = 'v2'
capture2 = 'v3'
capture3 = 'v3+UTR'


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

process BedToIntervalList {
        container 'broadinstitute/gatk:4.1.4.1'
        errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
        maxRetries 100

        input:
        tuple val(kitID), file(capture_bed)
        file reference
        file reference_dict

        output:
	
        tuple val("${kitID}"), file("${kitID}.interval_list")

        """
        gatk BedToIntervalList \
            -I ${capture_bed} \
	    -O ${kitID}.interval_list \
	    -SD ${reference}
        """
}



process PreprocessIntervals {
        container 'broadinstitute/gatk:4.1.4.1'
        errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
        maxRetries 100

        input:
        tuple val(kitID), file(capture_intervals)
        file reference
	file reference_dict
	file reference_index

        output:
        tuple val("${kitID}"), file("${kitID}.interval_list")

        """
	gatk PreprocessIntervals \
	    -L ${capture_intervals} \
    	    -R ${reference} \
    	    --bin-length 0 \
    	    --interval-merging-rule OVERLAPPING_ONLY \
    	    -O ${kitID}_processed.interval_list
        """
}

process samtoolsIndex {
	container "fredhutch/bwa:0.7.17"	


	input:
	tuple val(sampleID), val(kitID), val(type), val(patientID),  file(bam)	
	
	output:
	tuple val("${kitID}"), val("${patientID}"), val("${type}"), val("${sampleID}"), file("${bam}"), file("${bam}.bai")


	"""
	samtools index ${bam}
	"""
}


process CollectReadCounts {
        container 'broadinstitute/gatk:4.1.4.1'
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
        maxRetries 100

        input:
        tuple val(patientID), val(type), val(sampleID), file(bam), file(bam_index), val(kitID), file(capture_intervals)
        file reference
        //file reference_index
        //file reference_dict

        output:
        tuple val("${sampleID}"), val("${type}"), val("${kitID}"), val("$patientID}"), file("${sampleID}.counts.hdf5")

        """
	gatk CollectReadCounts \
    	     -I ${bam} \
    	     -L  ${capture_intervals} \
    	     --interval-merging-rule OVERLAPPING_ONLY \
    	     -O ${sampleID}.counts.hdf5
        """
}


process CreateReadCountPanelOfNormals{
        container 'broadinstitute/gatk:4.1.4.1'
        errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
        maxRetries 100

        input:
	tuple val(kitID), file(normal_list), file(annotated)
		
        output:
	tuple val("${kitID}"), file("${kitID}.cnv.pon.hdf5")

	"""
	#!/usr/bin/python
	
	import os
	
	files = ' ${normal_list}'
	files = files.replace(' ', ' -I ')
	print(files)
	script1 = 'gatk CreateReadCountPanelOfNormals ' 
	script2 = ' --annotated-intervals ${annotated} -O ${kitID}.cnv.pon.hdf5'

	script_f = script1 + files + script2

	os.system(script_f)
        """
}

process AnnotateIntervals {
        container 'broadinstitute/gatk:4.1.4.1'
        errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
        maxRetries 100

	input:
	tuple val(kitID), file(interval)
	file reference
	file reference_dict
	file reference_index


	output:
	tuple val("${kitID}"), file("${kitID}_annotated.tsv")

          """
     	  gatk AnnotateIntervals \
          -R ${reference} \
          -L ${interval} \
          --interval-merging-rule OVERLAPPING_ONLY \
          -O ${kitID}_annotated.tsv
	  """
}


process DenoiseReadCounts{
        container 'broadinstitute/gatk:4.1.4.1'
        errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
        maxRetries 100

        input:
	tuple val(kitID), file(pon), val(sampleID), val(patientID), file(counts)
	
        output:
	tuple val("${sampleID}"), val("${patientID}"), file("${sampleID}.standardizedCR.tsv"), file("${sampleID}.denoisedCR.tsv")

	"""
	gatk --java-options "-Xmx12g" DenoiseReadCounts \
	    -I ${counts} \
    	    --count-panel-of-normals ${pon} \
    	    --standardized-copy-ratios ${sampleID}.standardizedCR.tsv  --denoised-copy-ratios ${sampleID}.denoisedCR.tsv
        """
}


process CreateSequenceDictionary {
        container 'broadinstitute/gatk:4.1.4.1'
        errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
        maxRetries 100

        input:
        file reference


        output:
	file("${reference}.dict")

        """
	gatk CreateSequenceDictionary -R ${reference} -O ${reference}.dict
        """
}


process PlotDenoisedCopyRatios{
        container 'broadinstitute/gatk:4.1.4.1'
        errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
        maxRetries 100

	input:
	tuple val(sampleID), file(standard), file(denoised)
	file dict
	
	output:
	tuple val("${sampleID}"), path("plots")


	"""
     mkdir plots
     gatk PlotDenoisedCopyRatios \
          --standardized-copy-ratios ${standard} \
          --denoised-copy-ratios ${denoised} \
          --sequence-dictionary ${dict} \
          --minimum-contig-length 46709983 \
    	  --output ./plots \
    	  --output-prefix ${sampleID}
	"""

}


process CollectAllelicCounts {
	container 'broadinstitute/gatk:4.1.4.1'
        errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
        maxRetries 100	
	
	input:
        tuple val(patientID), val(sampleID), val(type), file(bam), file(bam_index), val(kitID), file(capture_intervals)
	file reference
	file reference_index
	file reference_dict

	output:
	tuple val("${patientID}"), val("${kitID}"), val("${sampleID}"), file("${sampleID}.allelicCounts.tsv")

	"""

	gatk CollectAllelicCounts \
          -I ${bam} \
          -R ${reference} \
          -L ${capture_intervals} \
          -O ${sampleID}.allelicCounts.tsv

	 """
}


//ADD CollectAllelicCounts.out to input 
// takes tuple in shape of DenoiseReadCounts.out
process ModelSegments{
        container 'broadinstitute/gatk:4.1.4.1'
        errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
        maxRetries 100

        input:
        tuple val(sampleID), file(standard), file(denoised), file(allelic)
  
        output:
        tuple val("${sampleID}"), path("${sampleID}_modelSeg")


        """
	mkdir ${sampleID}_modelSeg

	gatk --java-options "-Xmx30G" ModelSegments \
          --denoised-copy-ratios ${denoised} \
          --output-prefix ${sampleID} \
	  --allelic-counts ${allelic} \
          -O ${sampleID}_modelSeg

        """

}

process ModelSegmentsMatched {
	container 'broadinstitute/gatk:4.1.4.1'
        errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
        maxRetries 100

        input:
        tuple val(sampleID), file(standard), file(denoised), file(allelic), file(normal_allelic)


        output:
        tuple val("${sampleID}"), path("${sampleID}_modelSeg")


        """
        mkdir ${sampleID}_modelSeg

        gatk --java-options "-Xmx30G" ModelSegments \
          --denoised-copy-ratios ${denoised} \
          --output-prefix ${sampleID} \
	  --normal-allelic-counts ${normal_allelic}
          --allelic-counts ${allelic} \
          -O ${sampleID}_modelSeg

        """
}

//takes Copy-ratio-segments .cr.seg from modelsegments
process CallCopyRatioSegments {
	container 'broadinstitute/gatk:4.1.4.1'
        errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
        maxRetries 100

	input:
	tuple val(sampleID), path(modelSeg)

	output:
	tuple val("${sampleID}"), file("${sampleID}.called.seg")
	
	"""
	echo ${modelSeg}/*
	gatk CallCopyRatioSegments \
          -I ${modelSeg}/${sampleID}.cr.seg \
          -O ${sampleID}.called.seg \
	  --neutral-segment-copy-ratio-upper-bound 1.4 \
	  --neutral-segment-copy-ratio-lower-bound 0.6

 	"""

}


process PlotModeledSegments {
        container 'broadinstitute/gatk:4.1.4.1'
        errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
        maxRetries 100

        input:
        tuple val(sampleID), path(modelSeg), path(denoised)
        path contig_dict

        output:
        tuple val("${sampleID}"), path("${sampleID}_PlotModeledSegments")

        """
        gatk PlotModeledSegments \
          --denoised-copy-ratios ${denoised} \
          --allelic-counts ${modelSeg}/${sampleID}.hets.tsv \
          --segments ${modelSeg}/${sampleID}.modelFinal.seg \
          --sequence-dictionary ${contig_dict} \
          --output-prefix ${sampleID} \
          -O ${sampleID}_PlotModeledSegments

        """

}

process FuncotateSegments{
        container 'broadinstitute/gatk:4.1.4.1'
        errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
        maxRetries 100

        input:
	tuple val(sampleID), path(modelSeg)
        path reference
        path reference_dict
        path reference_index
        path data_source

        output:
        tuple val("${sampleID}"), path("${sampleID}_func.seg")

        """

        gatk FuncotateSegments \
        --data-sources-path ${data_source} \
        --ref-version hg38 \
        --output-file-format SEG \
        -R ${reference} \
        --segments ${modelSeg}/${sampleID}.cr.seg \
        -O ${sampleID}_func.seg

     	"""
}



process getSeg {
        container "ubuntu"

        input:
        tuple val(sampleID), path(modelSeg)

        output:
        path("${sampleID}.seg")

        """
        grep -v '^\\@' ./${modelSeg}/${sampleID}.cr.seg > segfile.seg
        sed '1d' segfile.seg > segfile1.seg
        awk -F "\t" '{\$1="${sampleID}" FS \$1;}1' OFS="\t" segfile1.seg > segfile2.seg
	sed '/chr.\\+_/d' segfile2.seg > ${sampleID}.seg
	"""

}


process run_gistic {
	container 'shixiangwang/gistic:latest'

	input:
	path(segs)

	output:
	path("gistic_out")


	"""
	CWD=\$(pwd)
	echo \$CWD
	cat *.seg > combined.seg
        mv combined.seg /opt/GISTIC/
        cd /opt/GISTIC
        mkdir gistic_out

#	mv /opt/GISTIC/gp_gistic2_from_seg .
	/opt/GISTIC/gistic2 \
	-b ./gistic_out \
	-seg ./combined.seg \
	-refgene /opt/GISTIC/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat \
	-rx 0 \
	-genegistic 1 \
	--armpeel 1 \
	-conf 0.95 \
	-ta 0.2 \
	-td 0.2 \
	
	mv gistic_out \$CWD 	
	cd \$CWD
	
	echo gistic_out/*
	"""

}




workflow gatkCNV_wf {

	 take:
		 beds_ch
		 reference 
         	 reference_dict 
		 reference_index 
		 contig_dict
		 sams_ch
		 bams_ch 
		 matched
		 
	 main:
	 CreateSequenceDictionary(reference)

	 BedToIntervalList(beds_ch, reference, reference_dict)
	 
 	 PreprocessIntervals(BedToIntervalList.out, reference, reference_dict, reference_index)

	 AnnotateIntervals(PreprocessIntervals.out, reference, reference_dict, reference_index)
	 
	 samtoolsIndex(bams_ch)
	 
	 //samtoolsIndex.out has 6 tuple with kitID, patientID ....
	 //want to remove r[2], it is a second kitID
	 
	 crossed = PreprocessIntervals.out.cross(samtoolsIndex.out).map { r -> r.flatten() }.map { r -> [r[3], r[4], r[5], r[6], r[7], r[0], r[1]] }
	 //comes out val("${kitID}"), val("${patientID}"), val("${type}"), val("${sampleID}"), file("${bam}"), file("${bam}.bai")     
	crossed.view()


//################   Get Read counts for normal and tumor separately, group the normals by kitID to make PON    #########################

	CollectReadCounts(crossed, reference)
        normies = CollectReadCounts.out.filter { it[1] == 'Normal'}.map { r -> [r[2], r[3]] }.groupTuple()
        tumors = CollectReadCounts.out.filter { it[1] == 'Tumor'}.map { r -> [r[2], r[0], r[3]] }
	normal_and_annotations = normies.cross(AnnotateIntervals.out).map{ r -> [r[0][0], r[0][1], r[1][1]] }
	CreateReadCountPanelOfNormals(normal_and_annotations)
	tumors_norm_ch = CreateReadCountPanelOfNormals.out.cross(tumors).map {r -> r.flatten()}
                                                             .map { r -> [r[0], r[1], r[3], r[4]] }

	DenoiseReadCounts(tumors_norm_ch)
	PlotDenoisedCopyRatios(DenoiseReadCounts.out, CreateSequenceDictionary.out)

//#######################################################################################################################################



//###########################   Get Allelic counts for normals and Tumors ################################################################

	if(matched) {

		
		normalAllelicCounts = CollectAllelicCounts(crossed.filter { it[1] == 'Normal' }, reference, reference_index, reference_dict).
		tumorAllelicCounts = CollectAllelicCounts(crossed.filter { it[1] == 'Tumor' }, reference, reference_index, reference_dict)

		tumorJoinNormal = tumorAllelicCounts.out.join(normalAlleliccounts.out, remainder: false, by : [0,1])
				  						       		  .map { r -> [r[2], r[3], r[5]] }
		
		denoisedJoinAllelic = DenoiseReadCounts.out.join(tumorJoinNormal)

		modeled = ModelSegmentsMatched(denoisedJoinAllelic)
	}
	else {	
	        AllelicCounts = CollectAllelicCounts(crossed, reference, reference_index, reference_dict)
       	        modeled = ModelSegments(DenoiseReadCounts.out.join(AllelicCounts))
	}

//#######################################################################################################################################	 	 



//####################################   Call and Plot modelled segments ###############################################################
	
	CallCopyRatioSegments(modeled)
	modPlotData = modeled.join(DenoiseReadCounts.out.map { r -> [r[0], r[2]] })
        PlotModeledSegments(modPlotData, CreateSequenceDictionary.out)

//#######################################################################################################################################	 	 


//	getSeg(modeled)
//	run_gistic(getSeg.out.collect())
//      DownloadData()
//      FuncotateSegments(modeled, reference, reference_dict, reference_index, DownloadData.out)

	
	 emit:
	        plots = PlotDenoisedCopyRatios.out
                modelsegs = ModelSegments.out
                denoised = DenoiseReadCounts.out
                modelPlots = PlotModeledSegments.out
		calledSegs = CallCopyRatioSegments.out
		download = DownloadData.out
		segs = getSeg.out
		gistic = run_gistic.out
//              func = FuncotateSegments.out

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

	 matched = params.matched

        bams_ch = Channel
            .fromPath(params.input_csv)
            .splitCsv(header:true)
            .map{ row-> tuple(row.sampleID, row.kitID, row.type, row.patientID, file(row.bam)) }

        beds_ch = Channel
            .fromPath(params.input_beds)
            .splitCsv(header:true)
            .map{ row-> tuple(row.kitID, file(row.capture_bed)) }
	
	
        sams_ch = Channel.empty()

        main:

          gatkCNV_wf(beds_ch, reference, reference_dict, reference_index, contig_dict, sams_ch, bams_ch, matched)


        publish:

         gatkCNV_wf.out.plots to : "./plots4/"
         gatkCNV_wf.out.modelPlots to : "./plots4/"

}