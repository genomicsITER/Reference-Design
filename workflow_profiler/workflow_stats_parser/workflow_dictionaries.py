###################################
# Copyright (c) Intel Corporation    
# PLEASE DO NOT SHARE WITHOUT NDA   
###################################

from collections import OrderedDict


# The user needs to add their own dictionary here.
# The step names for the dictionary should match the stage names in the pipeline script
# For example: if there exists run..stage1.* , run..stage2.* and run..stage3.* folders, the dictionary would be constructed as in the example provided below for sample data. 

# This dictionary is used to test the existing sample data
sample_dict = OrderedDict([('Stage1', 'stage1'), 
                           ('Stage2','stage2'), 
                           ('Stage3','stage3')])
##OHSU Ordered Dictionary 
ohsu_dict = OrderedDict([('bwa mem', 'BwaMem'),
                         ('samtools view', 'SamtoolsView'), 
                         ('samtools sort', 'SamtoolsSort'), 
                         ('MarkDuplicates', 'MarkDuplicates')]) 
 #                        ('RealignerTargetCreator', 'RealignerTargetCreator'), 
 #                        ('IndelRealigner', 'IndelRealigner'), 
 #                        ('BaseRecalibrator', 'BaseRecalibrator'), 
 #                        ('PrintReads', 'PrintReads')])


##GATK Best Practices Ordered Dictionary 
#broad_gatk_dict = OrderedDict([('bwa mem', 'BwaMem')])

broad_gatk_dict = OrderedDict([('bwa mem', 'BwaMem'),
                         ('SortSam', 'SortSam'), 
                         ('MarkDuplicates', 'MarkDuplicates'), 
                         ('RealignerTargetCreator', 'RealignerTargetCreator'), 
                         ('IndelRealigner', 'IndelRealigner'), 
                         ('BaseRecalibrator', 'BaseRecalibrator'), 
                         ('PrintReads', 'PrintReads'),
                         ('HaplotypeCaller', 'HaplotypeCaller')])


#BROAD Ordered dictionary
broad_dict = OrderedDict([('SamToFastq', 'SamToFastq'),
                          ('bwa1', 'BwaAln1'),
                          ('bwa2', 'BwaAln2'),
                          ('sampe', 'Sampe'),
                          ('MergeBamAlignment', 'MergeBamAlignment'),
                          ('MarkDuplicates', 'MarkDuplicates'),
                          ('IndelRealigner', 'IndelRealigner'),
                          ('BaseRecalibrator', 'BaseRecalibrator'),
                          ('PrintReads', 'PrintReads')])

#Broad Haplotypecaller dict
broad_hc_dict = OrderedDict([('HaplotypeCaller', 'HaplotypeCaller')])

#Broad gatk Variant Discovery
broad_variant_dict = OrderedDict([('GenotypeGVCFs', 'GenotypeGVCFs'),
                          ('VariantFiltration', 'VarFilter'),
                          ('VariantRecalibratorSNP', 'VarRecalSNP'),
                          ('ApplyRecalibrationSNP', 'ApplyRecalSNP'),
                          ('VariantRecalibratorIndel', 'VarRecalIndel'),
                          ('ApplyRecalibratioIndel', 'ApplyRecalIndel')])


#Broad gatk pipline but without haplotypcaller
broad_minus_hc_dict = OrderedDict([('bwa mem', 'BwaMem'),
                         ('SortSam', 'SortSam'),
                         ('MarkDuplicates', 'MarkDuplicates'),
                         ('RealignerTargetCreator', 'RealignerTargetCreator'),
                         ('IndelRealigner', 'IndelRealigner'),
                         ('BaseRecalibrator', 'BaseRecalibrator'),
                         ('PrintReads', 'PrintReads')])
##unordered dictionary for parsing pipline input parameter to one of the above pipeline dictionaries
# Format is:
#   'pipelineName': 'name_of_dictionary'
#
# The pipelineName will be what is entered on the command line with the -P option,
#   the name_of_dictionary is the name of the OrderedDict you have defined above
#
# If you add a new dictionary above, be sure to include it in this list
workflow_parse_dict = {
    'sample': 'sample_dict',
    'ohsu': 'ohsu_dict',
    'gatk_best_practices': 'broad_gatk_dict',
    'gatk_variant_dis': 'broad_variant_dict',
    'broad_hc': 'broad_hc_dict',
    'broad': 'broad_dict',
    'broad_minus_hc': 'broad_minus_hc_dict'
}
