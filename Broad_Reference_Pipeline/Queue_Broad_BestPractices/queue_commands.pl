my $nfs_ExampleIndelRealigner = $nfsDir . "ExampleIndelRealigner.scala";
my $nfs_ExampleBaseRecalibrator = $nfsDir . "ExampleBaseRecalibrator.scala";
my $nfs_ExamplePrintReads = $nfsDir . "ExamplePrintReads.scala";
my $nfs_queueDir = $nfsDir . "/.queue";



# #run_and_log "$myLDPRELOAD $javacmd $javaopts -Djava.io.tmpDir=$tmpDir -jar $gatk -T IndelRealigner  -rf NotPrimaryAlignment -R $refgenomeFastaFile -targetIntervals $realnInterval -known:indels,vcf $dbSNPindel -I $bamDupRemFile -o $bamRealignFile";
run_and_log("$javacmd -jar dist/Queue.jar -R $refgenomeFastaFile -I $nfs_bamDupRemFile  -D $dbSNPvcf -indels $dbSNPindel -S $nfs_ExampleIndelRealigner -l DEBUG -run -jobRunner CMPShell"
run_and_log "$samtools index $nfs_bamRealignFile";




#Multi-threading ERROR MESSAGE: We have temporarily disabled the ability to run BaseRecalibrator multi-threaded for performance reasons.
#ERROR MESSAGE: GATK Lite does not support all of the features of the full version: base insertion/deletion recalibration is not supported, please use the --disable_indel_quals argument
#ERROR MESSAGE: Invalid command line: DinucCovariate has been retired.  Please use its successor covariate ContextCovariate instead, which includes the 2 bp (dinuc) substitution model of the retired DinucCovariate as well as an indel context to model the indel error rates
run_and_log("$javacmd -jar dist/Queue.jar -I $nfs_bamRealignFile -R $refgenomeFastaFile -D $dbSNPvcf -dp illumina -S $nfs_ExampleBaseRecalibrator -l DEBUG -run -jobRunner CMPShell");





run_and_log("$javacmd -jar dist/Queue.jar -R $refgenomeFastaFile -I $nfs_bamRealignFile -B $nfs_recalOut  -l DEBUG -S $nfs_ExamplePrintReads -run -jobRunner CMPShell");

run_and_log "$myLDPRELOAD $samtools index $nfs_finalBam";
