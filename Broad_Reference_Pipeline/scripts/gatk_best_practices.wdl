task bwa_mem {
    File reference_genome
    File fasta_amb
    File fasta_ann
    File fasta_bwt
    File fasta_pac
    File fasta_rbwt
    File fasta_rpac
    File fasta_rsa
    File fasta_sa
    File global_indels
    File global_dbsnp
    File bwa_mem_read1_input_file
    File bwa_mem_read2_input_file
    File g_indels_vcf_idx
    File dbsnp_all_vcf_idx
	
    command {
        bwa mem -t 4 -M -a -R '@RG\tID:r\tLB:l\tSM:s\tPL:ILLUMINA' ${reference_genome} ${bwa_mem_read1_input_file} ${bwa_mem_read2_input_file}
    }
    output {    
        File response = stdout()
    }
    runtime {
        tool: "bwa-mem"
        strategy: "routed_file"
    }
}

task picard_sort_sam {
    File input_sam
    String output_dir
	
    command {
         picard SortSam I=${input_sam} O="${output_dir}sorted_output.bam" SO=coordinate CREATE_INDEX=true
    }
    output {
        File response = "sorted_output.bam"
    }
        
    runtime {
        tool: "picard-sortsam"
        strategy: "routed_file"
    }
}

task picard_mark_duplicates {
    File picard_markduplicates_input_bam
    String output_dir
	
    command {
         picard MarkDuplicates I=${picard_markduplicates_input_bam} O="${output_dir}picard_markduplicates_output.bam" M="${output_dir}picard_markduplicates_output.metrics" TMP_DIR=${output_dir} CREATE_INDEX=true
    }
    output {
        File picard_markduplicates_bam_file = "picard_markduplicates_output.bam"
        File picard_markduplicates_bai_file = "picard_markduplicates_output.bai"
    }
        
    runtime {
        tool: "picard-markduplicates"
        strategy: "routed_file"
    }
}

task gatk_realigner_target_creator {
    File reference_genome
    File reference_genome_fai
    File dict
    File global_indels
    File global_dbsnp
    File realigner_target_creator_input_bam
    File realigner_target_creator_input_bai
    String output_dir
	
    command {
        GenomeAnalysisTK -T RealignerTargetCreator -nt 4 -R ${reference_genome} -o "${output_dir}gatk_realignertargetcreator.intervals" -known:indels,vcf ${global_indels} -I ${realigner_target_creator_input_bam} 
    }
    output {
        File response = "gatk_realignertargetcreator.intervals"
    }
        
    runtime {
        tool: "gatk-realignertargetcreator"
        strategy: "routed_file"
    }
}

task gatk_indel_realigner {
    File reference_genome
    File reference_genome_fai
    File dict
    File global_indels
    File global_dbsnp
    File indel_realigner_input_bam
    File indel_realigner_input_bai
    File indel_realigner_input_intervals
    String output_dir
	
    command {
        GenomeAnalysisTK -T IndelRealigner --filter_bases_not_stored -R ${reference_genome} -targetIntervals ${indel_realigner_input_intervals} -known:indels,vcf ${global_indels} -I ${indel_realigner_input_bam} -o "${output_dir}gatk_indelrealigner_output.bam"
    }
    output {
        File response_bam = "gatk_indelrealigner_output.bam"
		File response_bai = "gatk_indelrealigner_output.bai"
    }
        
    runtime {
        tool: "gatk-indelrealigner"
        strategy: "routed_file"
    }
}

task samtools_index_indel_realigner {

    File samtools_index_input_file
    String output_dir
	
    command {
        samtools index ${samtools_index_input_file} "${output_dir}gatk_indelrealigner_output.bai"
    }

    output {
        File response = "gatk_indelrealigner_output.bai"
    }
    
    runtime {
        tool: "samtools-index"
        strategy: "routed_file"
    }
}

task gatk_base_recalibrator {
    
    File reference_genome
    File reference_genome_fai
    File dict
    File global_indels
    File global_dbsnp
    File base_recalibrator_input_bam
    File base_recalibrator_input_bai
    String output_dir
	
    command {
	GenomeAnalysisTK -T BaseRecalibrator -I ${base_recalibrator_input_bam}  -R ${reference_genome} -knownSites:mask,vcf ${global_dbsnp} -o "${output_dir}gatk_baserecalibrator.grp"
    }
    output {
        File response = "gatk_baserecalibrator.grp"
    }
        
    runtime {
        tool: "gatk-baserecalibrator"
        strategy: "routed_file"
    }
}

task gatk_print_reads {
    
    File reference_genome
    File reference_genome_fai
    File dict
    File global_indels
    File global_dbsnp
    File print_reads_input_bam
    File print_reads_input_bai
    File print_reads_input_grp
    String output_dir
	
    command {
        GenomeAnalysisTK -T PrintReads -R ${reference_genome} -I ${print_reads_input_bam} -BQSR ${print_reads_input_grp} -o "${output_dir}gatk_printreads_output.bam"
    }
    output {
        File response_bam = "gatk_printreads_output.bam"
        File response_bai = "gatk_printreads_output.bai"
    }
        
    runtime {
        tool: "gatk-printreads"
        strategy: "routed_file"
    }
}

task gatk_haplotype_caller {
    
    File reference_genome
    File reference_genome_fai
    File dict
    File global_indels
    File global_dbsnp
    File haplotype_caller_input_bam
    File haplotype_caller_input_bai
    String output_dir
	
    command {
        GenomeAnalysisTK -T HaplotypeCaller -R ${reference_genome} -I ${haplotype_caller_input_bam} -ERC GVCF --variant_index_type LINEAR --variant_index_parameter 128000 -o "${output_dir}gatk_haplotypecaller_output.vcf"
    }
    output {
        File response_vcf = "gatk_haplotypecaller_output.vcf"
    }
        
    runtime {
        tool: "gatk-haplotypecaller"
        strategy: "routed_file"
    }
}


workflow Broad_GATK_BestPractices {
    
    File reference_genome
    File global_indels
    File global_dbsnp
    File reference_genome_fai
    File g_indels_vcf_idx
    File dbsnp_all_vcf_idx
    File dict
    File fasta_amb
    File fasta_ann
    File fasta_bwt 
    File fasta_pac
    File fasta_rbwt
    File fasta_rpac
    File fasta_rsa
    File fasta_sa
    String global_output_dir
	 
    call bwa_mem
    {
	input: 
	    reference_genome = reference_genome,
	    fasta_amb = fasta_amb,
	    fasta_ann = fasta_ann,
	    fasta_bwt = fasta_bwt,
	    fasta_pac = fasta_pac,
	    fasta_rbwt = fasta_rbwt,
	    fasta_rpac = fasta_rpac,
	    fasta_rsa = fasta_rsa,
	    fasta_sa = fasta_sa,
	    global_indels = global_indels,
	    global_dbsnp = global_dbsnp,
	    g_indels_vcf_idx = g_indels_vcf_idx,
    	    dbsnp_all_vcf_idx = dbsnp_all_vcf_idx
    }
    call picard_sort_sam 
    {
        input: 
	    input_sam = bwa_mem.response,
	    output_dir = global_output_dir
    }
    call picard_mark_duplicates 
    {
        input:
	     picard_markduplicates_input_bam = picard_sort_sam.response,
	     output_dir = global_output_dir
    }
    call gatk_realigner_target_creator 
    {
        input:
	     reference_genome = reference_genome,
             reference_genome_fai = reference_genome_fai,
	     dict = dict,
	     global_indels = global_indels,
	     global_dbsnp = global_dbsnp,
	     realigner_target_creator_input_bam = picard_mark_duplicates.picard_markduplicates_bam_file, realigner_target_creator_input_bai = picard_mark_duplicates.picard_markduplicates_bai_file,
	     output_dir = global_output_dir 
    }
    call gatk_indel_realigner 
    {
	input:
            reference_genome = reference_genome,
            reference_genome_fai = reference_genome_fai,
            dict = dict,
            global_indels = global_indels,
            global_dbsnp = global_dbsnp,
            indel_realigner_input_bam = picard_mark_duplicates.picard_markduplicates_bam_file,
            indel_realigner_input_bai = picard_mark_duplicates.picard_markduplicates_bai_file,
            indel_realigner_input_intervals = gatk_realigner_target_creator.response,
            output_dir = global_output_dir
    }
    call samtools_index_indel_realigner
    {
        input: 
	    samtools_index_input_file = gatk_indel_realigner.response_bam,
	    output_dir = global_output_dir
    }
    call gatk_base_recalibrator
    {
	input:
	    reference_genome = reference_genome,
	    reference_genome_fai = reference_genome_fai,
	    dict = dict,
	    global_indels = global_indels,
	    global_dbsnp = global_dbsnp,
	    base_recalibrator_input_bam = gatk_indel_realigner.response_bam,
	    base_recalibrator_input_bai = gatk_indel_realigner.response_bai,
	    output_dir = global_output_dir
    }
    call gatk_print_reads 
    {
        input:
	    reference_genome = reference_genome,
	    reference_genome_fai = reference_genome_fai,
	    dict = dict,
	    global_indels = global_indels,
	    global_dbsnp = global_dbsnp,
	    print_reads_input_bam = gatk_indel_realigner.response_bam,
	    print_reads_input_bai = gatk_indel_realigner.response_bai,
	    print_reads_input_grp = gatk_base_recalibrator.response,
	    output_dir = global_output_dir
    }
    call gatk_haplotype_caller 
    {
        input:
	    reference_genome = reference_genome,
	    reference_genome_fai = reference_genome_fai,
	    dict = dict,
	    global_indels = global_indels,
	    global_dbsnp = global_dbsnp,
	    haplotype_caller_input_bam = gatk_print_reads.response_bam,
	    haplotype_caller_input_bai = gatk_print_reads.response_bai,
	    output_dir = global_output_dir
    }
}
